#include "BrexWitnessChecker_UTEST.hpp"
#include "BrexWitnessChecker.hpp"
#include "common/lightCircLib/lightCircPoly.hpp"

#include "PCP/Tests/ACSP_PCP/PCPWorkflowTests.h"

#include <gtest/gtest.h>

namespace PCP_UTESTS{

using PCP_Project::BREXWitnessChecker;
using PCP_Project::BREXWitness;
using PCP_Project::BREXInstance;
using PCP_Project::ConstraintSys;
using PCP_Project::Sequence;
using PCP_Project::lightCircPoly;
using Algebra::FiniteField;
using Algebra::FieldElement;
using Algebra::PolynomialInterface;
using Algebra::UnivariatePolynomialInterface;
using Algebra::UnivariatePolynomialGeneral;
using Algebra::PolynomialDegree;
using Algebra::one;
using Algebra::generateRandom;
using Algebra::elementsSet_t;
using Infrastructure::POW2;
using std::pair;
using std::rand;
using std::vector;
using std::unique_ptr;
using std::move;

/****************************************************************
 *                Auxiliary classes templates
 ****************************************************************/
template<class element_t>
class generateElement {
public:
    virtual element_t operator()() const = 0;
};

/**
 * @class specificSequence
 * @brief A mapping of \f$\mathbb{N}\f$ into the <image_t>
 * such that the first \f$len\f$ integers are mapped to some given
 * constants, and the rest are mapped to a given constant as well
 */
template <class image_t>
class expliciteSequence : public Sequence<image_t> {
public:
    /**
     * @brief   The constructor
     * @param   n boundary of indexes that can be mapped to different constant elements (namely \f$len\f$)
     */
    expliciteSequence(size_t n, const generateElement<image_t>& gen,const image_t& rest):order_(n),rest_(rest){
        for (size_t i=0; i<order_.size(); i++){
            order_[i] = gen();
        }
    }

    /**
     * @brief   The mapping of integers to <image_t>
     * @param   index some integer
     * @return  its mapping
     */
    image_t getElementByIndex(typename Sequence<image_t>::index_t index)const {
        if (index < order_.size()) return order_[index];
        else return rest_;
    }

private:

    /**
     * The mapping is represented using a vector of <image_t> elements,
     * if an integer is in the domain of the vector coordinates
     * it is mapped to the vectors value in that coordinate,
     * otherwise it is mapped to the return value of rest_.
     */
    vector< image_t> order_;
    image_t rest_;
};

/****************************************************************
 *                Coloring auxiliary classes
 ****************************************************************/

/********************
 * Random assignment
 *******************/

class randomColorGen : public generateElement<BREXWitness::color_t> {
public:
    randomColorGen(size_t len):len_(len){};
    BREXWitness::color_t operator()() const{
        BREXWitness::color_t color;
        for (int i=0; i< len_; i++){
            color.push_back(generateRandom());
        }
        return color;
    }
private:
    size_t len_;
};

class randomColoring : public expliciteSequence<BREXWitness::color_t> {
public:
    randomColoring(size_t n, size_t vec_len):
        expliciteSequence<BREXWitness::color_t>(n,randomColorGen(vec_len), randomColorGen(vec_len)()){};
};

/****************************************************************
 *               Permutations auxiliary classes
 ****************************************************************/

/********************
 * Random permutation
 *******************/

class randOrderGen : public generateElement<size_t> {
    public:
        randOrderGen(size_t numElements): available(numElements){
            for(size_t i=0; i < numElements; i++) available[i] = i;
        }
        size_t operator()() const{
            assert(!available.empty());
            size_t elementIndex = rand() % available.size();
            size_t retElement = available[elementIndex];
            available.erase( available.begin() + elementIndex );
            return retElement;
        }
    private:
        mutable vector<size_t> available;
};

class randPermutation : public expliciteSequence<size_t> {
public:
    randPermutation(size_t numElements):
        expliciteSequence<size_t>(numElements,randOrderGen(numElements),0){};
};

class plusOneSequence : public Sequence<size_t>{
    public:
    
		size_t getElementByIndex(Sequence<size_t>::index_t index)const {
        return index+1;
    }
};

/****************************************************************
 *               Random padding for permutation
 ****************************************************************/
vector<FieldElement> getRandomPadding(const size_t& vectorLen){
    vector<FieldElement> res(vectorLen);
    for(auto& x : res) x = generateRandom();

    return res;
}

/****************************************************************
 *           Constraints System auxiliary classes
 ****************************************************************/

/********************
 * Allways saticfied
 *******************/
class allwaysSatisfiedSystem : public ConstraintSys{
public:
    allwaysSatisfiedSystem(const FiniteField& field, size_t varsAmount):
        field_(field), 
        varsAmount_(varsAmount), 
        firstUsed_(varsAmount==2){};
    const FiniteField contextField()const {return field_;}
    size_t varsAmount() const {return varsAmount_;}
    const polySet_t& constraints() const {return noPolys_;}
    bool varUsed(const size_t varId) const { 
        //verify the case where first var is
        //not routed, is handled
        if((varId ==0) && !firstUsed_) return false;
        
        //non trivial function
        //this is implemented this way 
        //so the embedding itself will be non trivial
        //This should simulate the expected case
        //where the first half of the variables belongs the
        //previous configuration + auxiliary variables used when 
        //"previous configuration is next configuration"
        //and the second half is the next configuration
        //and the auxiliary variables relevant to current check.
        //so the constraint system would probably say that it uses
        //only previous configuration,next configuration, and
        //auxiliary variables relevant to current check.
        return (varId <= (varsAmount_*0.75));
    }
private:
    FiniteField field_;
    size_t varsAmount_;
    const polySet_t noPolys_; //empty set
    const bool firstUsed_;
};

class settingSaticfyingSystem : public ConstraintSys {
public:
    settingSaticfyingSystem(
        const FiniteField& field,
        size_t varsAmount,
        size_t domainSize,
        const BREXWitness& witness,
        const BREXWitness::permutation_t& permutation,
        bool makeIncomplete = false)
        : field_(field), varsAmount_(varsAmount) {
            for (int i=0; i < domainSize; i++){
                size_t perm_img = permutation.getElementByIndex(i);
                BREXWitness::color_t c1 = witness.get_color(i);
                BREXWitness::color_t c2 = witness.get_color(perm_img);
                vector<FieldElement> assignment(c1);
                for(auto elem : c2){ assignment.push_back(elem);}
                addRootToTrie(assignment);
            }
            if (makeIncomplete){
            //remove one root, so the system
            //wont be saticfied by current setting
                size_t i = rand() % (domainSize-1);
                size_t perm_img = permutation.getElementByIndex(i);
                BREXWitness::color_t c1 = witness.get_color(i);
                BREXWitness::color_t c2 = witness.get_color(perm_img);
                vector<FieldElement> assignment(c1);
                for(auto elem : c2){ assignment.push_back(elem);}
                removeRootFromTrie(assignment);
            }
           generatePolysFromTrie();
        }
    const FiniteField contextField()const {return field_;}
    size_t varsAmount() const {return varsAmount_;}
    const polySet_t& constraints() const {return polys_;}
    ~settingSaticfyingSystem(){
        for(auto node : rootsTrie_) delete node;
    }
private:
    FiniteField field_;
    size_t varsAmount_;
    polySet_t polys_;

    struct trieNode{
        FieldElement val;
        vector<struct trieNode*> next;
        trieNode(FieldElement newVal):val(newVal){;}
        ~ trieNode(){
            for(struct trieNode* node : next){
                delete node;
            }
        }
    };

    vector<struct trieNode*> rootsTrie_;

    void addRootToTrie_rec(vector<FieldElement>::const_iterator start,vector<FieldElement>::const_iterator end, vector<struct trieNode*>& currLevel){
        if(start == end) return;

        //find if prefix is already in trie
        //if so, just add the needed suffix
        for(struct trieNode* node : currLevel){
            assert(node != NULL);
            if(node->val == *start) return addRootToTrie_rec(start+1,end,node->next);
        }

        //if prefix is not in trie, add it
        struct trieNode* newNode = new struct trieNode(*start);
        currLevel.push_back(newNode);
        return addRootToTrie_rec(start+1,end,newNode->next);

    }

    void addRootToTrie(const vector<FieldElement>& assignment){
        vector<struct trieNode*>& firstList = rootsTrie_;
        addRootToTrie_rec(assignment.begin(), assignment.end(), firstList);
    }

    bool removeRootFromTrie_rec(vector<FieldElement>::const_iterator start,vector<FieldElement>::const_iterator end, vector<struct trieNode*>& currLevel){
        //if finished reading the assignment, then is is found in the trie
        if(start == end) return true;

        //find if prefix is in trie
        //if found, remove it
        for(size_t i=0; i< currLevel.size(); i++){
            assert(currLevel[i] != NULL);
            //if the suffix is found, remove it from the trie
            if((currLevel[i]->val == *start) && removeRootFromTrie_rec(start+1,end,currLevel[i]->next)){
                //if we are currently in the last node,
                //or there is only one way to continue 
                //(it must be the way only to the assignment, otherwise false is returned)
                //for sure it should be removed
                if(currLevel[i]->next.size() <= 1){
                    delete currLevel[i];
                    currLevel.erase(currLevel.begin() + i);
                    return true;
                }

                //we are on the way to more than only the assignment, so return false
                return false;
            }
        }

        //if prefix is not in trie, just return false (not found)
        return false;
    }

    void removeRootFromTrie(const vector<FieldElement>& assignment){
        vector<struct trieNode*>& firstList = rootsTrie_;
        removeRootFromTrie_rec(assignment.begin(), assignment.end(), firstList);
    }

    void generatePolysFromTrie_rec(const lightCircPoly& selector, const vector<struct trieNode*>& currLevel, const size_t currLevelIndex){
        //recursion end
        if(currLevelIndex >= varsAmount_) return;
        
        //build polynomial
        //that vanishes if and only if the value of
        //the element in the current index of an assignment correlates
        //with some root
        {
        elementsSet_t roots;
        for(const struct trieNode* node : currLevel){
            assert(node != NULL);
            roots.insert(node->val);
        }
        const UnivariatePolynomialGeneral uniPoly(roots);
        const lightCircPoly vanishesOnCurrElem(uniPoly);
        lightCircPoly selectorForCurrNode(selector);
        selectorForCurrNode.multiplyDistinct(vanishesOnCurrElem);

        //extend to fit the amount of variable
        vector<size_t> originalVarsLocations;
        for(int i=0; i<= currLevelIndex; i++){
            originalVarsLocations.push_back(i);
        }
        polys_.insert(polyPtr_t(new lightCircPoly(selectorForCurrNode,varsAmount_,originalVarsLocations)));
        }
     
        //For each posible value for the current variable,
        //build a selector polynomial that vanishes on every other
        //option prefix, but that value.
        //This polynomial would be a factor of a polynomial
        //that would vanish only on possible values for the next variable
        for(size_t currIndex = 0; currIndex < currLevel.size(); currIndex++){
            
            //gather roots
            elementsSet_t roots;
            for(const struct trieNode* node : currLevel){
                assert(node != NULL);
                if(node->val != currLevel[currIndex]->val){
                    roots.insert(node->val);
                }
            }

            //build polynomial
            const UnivariatePolynomialGeneral currLevelSelectorUniPoly(roots);
            const lightCircPoly currLevelSelector(currLevelSelectorUniPoly);
            lightCircPoly indexSelector(selector);
            indexSelector.multiplyDistinct(currLevelSelector);
            
            //call recursion
            generatePolysFromTrie_rec(indexSelector, currLevel[currIndex]->next, currLevelIndex+1);
        }
        
    }

    void generatePolysFromTrie(){
        vector<struct trieNode*>& firstList = rootsTrie_;

        //build first constrain polynomial
        //it vanishes if and only if the value of
        //the first element of an assignment correlates
        //with some root
        {
        elementsSet_t roots;
        for(const struct trieNode* node : firstList){
            assert(node != NULL);
            roots.insert(node->val);
        }
        const UnivariatePolynomialGeneral uniPoly(roots);
        const lightCircPoly vanishesOnFirstElem(uniPoly);

        //extend to fit the amount of variable
        vector<size_t> originalVarsLocations;
        originalVarsLocations.push_back(0);
        polys_.insert(polyPtr_t(new lightCircPoly(vanishesOnFirstElem,varsAmount_,originalVarsLocations)));
        }

        //For each posible value for the first variable,
        //build a polynomial that vanishes on every other
        //option, but that value.
        //This polynomial would be a factor of a polynomial
        //that would vanish only on possible values for the next variable
        for(size_t currIndex = 0; currIndex < firstList.size(); currIndex++){
            
            //gather roots
            elementsSet_t roots;
            for(const struct trieNode* node : firstList){
                assert(node != NULL);
                if(node->val != firstList[currIndex]->val){
                    roots.insert(node->val);
                }
            }

            //build polynomial
            const UnivariatePolynomialGeneral selector(roots);
            
            //call recursion
            generatePolysFromTrie_rec(lightCircPoly(selector), firstList[currIndex]->next, 1);
        }
        
    }
};

/***************************************************
 *
 * This generates a valid BREX pair,
 * with no constraints and no permutations.
 * the only relevant test here is boundary constraints.
 *
 ***************************************************/
pair<BREXInstance,BREXWitness> generate_valid_boundary(){

    /** constants **/
    const size_t extensionDegree = 64;
    const size_t vectorLen = 3;
    const short domainSizeIndicator = 3;
    const size_t domainSize = POW2(domainSizeIndicator) - 1;
    const size_t boundary_len = rand() % domainSize;

    /** define context field and make it global **/
    FiniteField contextField(extensionDegree);
    elementsSet_t fieldBasis = contextField.getBasis();
    contextField.setContext();

    /** construct witness **/
    BREXWitness::assignment_ptr assignment(new randomColoring(domainSize,vectorLen));
    BREXWitness::permutation_ptr permutation(new randPermutation(domainSize));
    BREXWitness witness(move(assignment), move(permutation));

    /** construct instance **/
    
    //empty constraints
    BREXInstance::constraintsPtr_t constraintsAssignment(new allwaysSatisfiedSystem(contextField,vectorLen*2));
    BREXInstance::constraintsPtr_t constraintsPermutation(new allwaysSatisfiedSystem(contextField,vectorLen*2));
    
    //random boundary constraints
    BREXInstance::boundaryConstraints_t boundaryConstraints;
    for(size_t i=0; i<boundary_len ; i++){
        const BREXInstance::point_t location(rand()%domainSize , rand()% vectorLen);
        const FieldElement val = witness.get_color(location.first)[location.second];
        boundaryConstraints[location] = val;
    }

    //construct the instance
    BREXInstance instance(contextField,vectorLen,domainSizeIndicator,move(constraintsAssignment), move(constraintsPermutation), boundaryConstraints, getRandomPadding(vectorLen));	

    /** Return the result **/
    return pair<BREXInstance,BREXWitness>(move(instance),move(witness));
}

/***************************************************
 *
 * This generates a invalid BREX pair,
 * with no constraints and no permutations.
 * the only relevant test here is boundary constraints.
 *
 ***************************************************/
pair<BREXInstance,BREXWitness> generate_invalid_boundary(){

    /** constants **/
    const size_t extensionDegree = 64;
    const size_t vectorLen = 3;
    const short domainSizeIndicator = 3;
    const size_t domainSize = POW2(domainSizeIndicator) - 1;
    const size_t boundary_len = 10 + rand() % domainSize;

    /** define context field and make it global **/
    FiniteField contextField(extensionDegree);
    elementsSet_t fieldBasis = contextField.getBasis();
    contextField.setContext();

    /** construct witness **/
    BREXWitness::assignment_ptr assignment(new randomColoring(domainSize,vectorLen));
    BREXWitness::permutation_ptr permutation(new randPermutation(domainSize));
    BREXWitness witness(move(assignment), move(permutation));

    /** construct instance **/
    
    //empty constraints
    BREXInstance::constraintsPtr_t constraintsAssignment(new allwaysSatisfiedSystem(contextField,vectorLen*2));
    BREXInstance::constraintsPtr_t constraintsPermutation(new allwaysSatisfiedSystem(contextField,vectorLen*2));
    
    //random boundary constraints
    BREXInstance::boundaryConstraints_t boundaryConstraints;
    for(size_t i=0; i<boundary_len ; i++){
        const BREXInstance::point_t location(rand()%domainSize , rand()% vectorLen);
        const FieldElement val = witness.get_color(location.first)[location.second];
        boundaryConstraints[location] = val + one();
    }

    //construct the instance
    BREXInstance instance(contextField,vectorLen,domainSizeIndicator,move(constraintsAssignment), move(constraintsPermutation), boundaryConstraints, getRandomPadding(vectorLen));	

    /** Return the result **/
    return pair<BREXInstance,BREXWitness>(move(instance),move(witness));
}


/***************************************************
 *
 * This generates a valid BREX pair,
 * with some empty constraints and valid permutations.
 * The boundary is empty.
 * The assignment is random.
 * 
 ***************************************************/
pair<BREXInstance,BREXWitness> generate_valid_permutations(){

    /** constants **/
    const size_t extensionDegree = 64;
    const size_t vectorLen = 3;
    const short domainSizeIndicator = 7;
    const size_t domainSize = POW2(domainSizeIndicator) - 1;

    /** define context field and make it global **/
    FiniteField contextField(extensionDegree);
    elementsSet_t fieldBasis = contextField.getBasis();
    contextField.setContext();

    /** construct witness **/
    BREXWitness::assignment_ptr assignment(new randomColoring(domainSize,vectorLen));
    
    //construct permutations
    BREXWitness::permutation_ptr permutation(new randPermutation(domainSize));
    
    //construct witness
    BREXWitness witness(move(assignment), move(permutation));

    /** construct instance **/
    
    //empty constraints vector
    BREXInstance::constraintsPtr_t constraintsAssignment(new allwaysSatisfiedSystem(contextField,vectorLen*2));
    BREXInstance::constraintsPtr_t constraintsPermutation(new allwaysSatisfiedSystem(contextField,vectorLen*2));
    
    //construct the instance
    BREXInstance instance(contextField,vectorLen,domainSizeIndicator,move(constraintsAssignment),move(constraintsPermutation), BREXInstance::boundaryConstraints_t(), getRandomPadding(vectorLen));	

    /** Return the result **/
    return pair<BREXInstance,BREXWitness>(move(instance),move(witness));
}

/***************************************************
 *
 * This generates some valid BREX pair
 * 
 ***************************************************/
pair<BREXInstance,BREXWitness> generate_valid_pair(){

    /** constants **/
    const size_t extensionDegree = 64;
    const size_t vectorLen = 2+rand()%10;
    const short domainSizeIndicator = 3;
    const size_t domainSize = POW2(domainSizeIndicator) - 1;

    /** define context field and make it global **/
    FiniteField contextField(extensionDegree);
    elementsSet_t fieldBasis = contextField.getBasis();
    contextField.setContext();

    /** construct witness **/
    BREXWitness::assignment_ptr assignment(new randomColoring(domainSize,vectorLen));
    
    //construct permutation
    BREXWitness::permutation_ptr permutation(new randPermutation(domainSize));
    
    //construct witness
    BREXWitness witness(move(assignment), move(permutation));

    /** construct instance **/
    
    //empty constraints
    BREXInstance::constraintsPtr_t constraintsAssignment(new allwaysSatisfiedSystem(contextField,vectorLen*2));
    BREXInstance::constraintsPtr_t constraintsPermutation(new allwaysSatisfiedSystem(contextField,vectorLen*2));
    
    //construct the instance
    BREXInstance instance(contextField,vectorLen,domainSizeIndicator,move(constraintsAssignment),move(constraintsPermutation), BREXInstance::boundaryConstraints_t(), getRandomPadding(vectorLen));	

    /** Return the result **/
    return pair<BREXInstance,BREXWitness>(move(instance),move(witness));
}

/***************************************************
 *
 * This generates an invalid BREX pair,
 * with some empty constraints and valid permutations,
 * except of one, which is not a permutation.
 * The boundary is empty.
 * The assignment is random.
 * 
 ***************************************************/
pair<BREXInstance,BREXWitness> generate_invalid_permutations(){

    /** constants **/
    const size_t extensionDegree = 64;
    const size_t vectorLen = 3;
    const short domainSizeIndicator = 7;
    const size_t domainSize = POW2(domainSizeIndicator) - 1;

    /** define context field and make it global **/
    FiniteField contextField(extensionDegree);
    elementsSet_t fieldBasis = contextField.getBasis();
    contextField.setContext();

    /** construct witness **/
    BREXWitness::assignment_ptr assignment(new randomColoring(domainSize,vectorLen));
    
    //construct (non)permutation
    //with some non-bijective mapping
    BREXWitness::permutation_ptr permutation(new randPermutation(domainSize-1));
    
    //construct witness
    BREXWitness witness(move(assignment), move(permutation));

    /** construct instance **/
    
    //empty constraints vector
    BREXInstance::constraintsPtr_t constraintsAssignment(new allwaysSatisfiedSystem(contextField,vectorLen*2));
    BREXInstance::constraintsPtr_t constraintsPermutation(new allwaysSatisfiedSystem(contextField,vectorLen*2));

    //construct the instance
    BREXInstance instance(contextField,vectorLen,domainSizeIndicator,move(constraintsAssignment), move(constraintsPermutation), BREXInstance::boundaryConstraints_t(), getRandomPadding(vectorLen));	

    /** Return the result **/
    return pair<BREXInstance,BREXWitness>(move(instance),move(witness));
}

/***************************************************
 *
 * This generates a valid BREX pair,
 * with some constraints and valid permutations.
 * The boundary is empty.
 * The assignment is random.
 * 
 * Generation method:
 * generates all parameters but the
 * constraint systems.
 * And defines each constraint system to be saticfied
 * exactly from the parameters given.
 * 
 ***************************************************/
pair<BREXInstance,BREXWitness> generate_valid_constraints(){

    /** constants **/
    const size_t extensionDegree = 64;
    const size_t vectorLen = 3;
    const short domainSizeIndicator = 3;
    const size_t domainSize = POW2(domainSizeIndicator) - 1;

    /** define context field and make it global **/
    FiniteField contextField(extensionDegree);
    elementsSet_t fieldBasis = contextField.getBasis();
    contextField.setContext();

    /** construct witness **/
    BREXWitness::assignment_ptr assignment(new randomColoring(domainSize,vectorLen));
    
    //construct permutations
    BREXWitness::permutation_ptr permutation(new randPermutation(domainSize));
    
    //construct witness
    BREXWitness witness(move(assignment), move(permutation));

    /** construct instance **/

    //construct constraints
    BREXInstance::constraintsPtr_t constraintsAssignment(new settingSaticfyingSystem(contextField,vectorLen*2,domainSize-1,witness, plusOneSequence()));
    BREXInstance::constraintsPtr_t constraintsPermutation(new settingSaticfyingSystem(contextField,vectorLen*2,domainSize,witness,witness.permutation()));
    
    //construct the instance
    BREXInstance instance(contextField,vectorLen,domainSizeIndicator,move(constraintsAssignment), move(constraintsPermutation), BREXInstance::boundaryConstraints_t(), getRandomPadding(vectorLen));	

    /** Return the result **/
    return pair<BREXInstance,BREXWitness>(move(instance),move(witness));
}

/***************************************************
 *
 * This generates an invalid BREX pair,
 * with constraints and valid permutations.
 * The boundary is empty.
 * The assignment is random.
 * 
 * Generation method:
 * generates all parameters but the
 * constraint systems.
 * And defines each constraint system to be satisfied
 * exactly from the parameters given,
 * except of one victim constraints system that is
 * chosen randomely, one victim index 'i',
 * such that the test of that index 
 * (color(i),color(perm(i))) will fail.
 *
 * In this case the victim is constraintsAssignment
 * 
 ***************************************************/
pair<BREXInstance,BREXWitness> generate_invalid_constraints_Assignment(){

    /** constants **/
    const size_t extensionDegree = 64;
    const size_t vectorLen = 3;
    const short domainSizeIndicator = 4;
    const size_t domainSize = POW2(domainSizeIndicator) - 1;

    /** define context field and make it global **/
    FiniteField contextField(extensionDegree);
    elementsSet_t fieldBasis = contextField.getBasis();
    contextField.setContext();

    /** construct witness **/
    BREXWitness::assignment_ptr assignment(new randomColoring(domainSize,vectorLen));
    
    //construct permutations
    BREXWitness::permutation_ptr permutation(new randPermutation(domainSize));
    
    //construct witness
    BREXWitness witness(move(assignment), move(permutation));

    /** construct instance **/

    //construct constraints
    BREXInstance::constraintsPtr_t constraintsAssignment(new settingSaticfyingSystem(contextField,vectorLen*2,domainSize-1,witness, plusOneSequence(),true));
    BREXInstance::constraintsPtr_t constraintsPermutation(new settingSaticfyingSystem(contextField,vectorLen*2,domainSize,witness,witness.permutation()));
    
    //construct the instance
    BREXInstance instance(contextField,vectorLen,domainSizeIndicator,move(constraintsAssignment), move(constraintsPermutation), BREXInstance::boundaryConstraints_t(), getRandomPadding(vectorLen));	

    /** Return the result **/
    return pair<BREXInstance,BREXWitness>(move(instance),move(witness));
}

/***************************************************
 *
 * This generates an invalid BREX pair,
 * with constraints and valid permutations.
 * The boundary is empty.
 * The assignment is random.
 * 
 * Generation method:
 * generates all parameters but the
 * constraint systems.
 * And defines each constraint system to be satisfied
 * exactly from the parameters given,
 * except of one victim constraints system that is
 * chosen randomely, one victim index 'i',
 * such that the test of that index 
 * (color(i),color(perm(i))) will fail.
 *
 * In this case the victim is constraintsPermutation
 * 
 ***************************************************/
pair<BREXInstance,BREXWitness> generate_invalid_constraints_Permutation(){

    /** constants **/
    const size_t extensionDegree = 64;
    const size_t vectorLen = 3;
    const short domainSizeIndicator = 4;
    const size_t domainSize = POW2(domainSizeIndicator) - 1;

    /** define context field and make it global **/
    FiniteField contextField(extensionDegree);
    elementsSet_t fieldBasis = contextField.getBasis();
    contextField.setContext();

    /** construct witness **/
    BREXWitness::assignment_ptr assignment(new randomColoring(domainSize,vectorLen));
    
    //construct permutations
    BREXWitness::permutation_ptr permutation(new randPermutation(domainSize));
    
    //construct witness
    BREXWitness witness(move(assignment), move(permutation));

    /** construct instance **/

    //construct constraints
    BREXInstance::constraintsPtr_t constraintsAssignment(new settingSaticfyingSystem(contextField,vectorLen*2,domainSize-1,witness, plusOneSequence()));
    BREXInstance::constraintsPtr_t constraintsPermutation(new settingSaticfyingSystem(contextField,vectorLen*2,domainSize,witness,witness.permutation(),true));
    
    //construct the instance
    BREXInstance instance(contextField,vectorLen,domainSizeIndicator,move(constraintsAssignment), move(constraintsPermutation), BREXInstance::boundaryConstraints_t(), getRandomPadding(vectorLen));	

    /** Return the result **/
    return pair<BREXInstance,BREXWitness>(move(instance),move(witness));
}

/***************************************************
 *
 * This generates an invalid BREX pair,
 * with constraints and valid permutations.
 * The boundary is empty.
 * The assignment is random.
 * 
 * Generation method:
 * generates all parameters but the
 * constraint systems.
 * And defines each constraint system to be satisfied
 * exactly from the parameters given,
 * except of one victim constraints system that is
 * chosen randomely, one victim index 'i',
 * such that the test of that index 
 * (color(i),color(perm(i))) will fail.
 *
 * In this case the victims are both constraintsAssignment & constraintsPermutation
 * 
 ***************************************************/
pair<BREXInstance,BREXWitness> generate_invalid_constraints_both(){

    /** constants **/
    const size_t extensionDegree = 64;
    const size_t vectorLen = 3;
    const short domainSizeIndicator = 4;
    const size_t domainSize = POW2(domainSizeIndicator) - 1;

    /** define context field and make it global **/
    FiniteField contextField(extensionDegree);
    elementsSet_t fieldBasis = contextField.getBasis();
    contextField.setContext();

    /** construct witness **/
    BREXWitness::assignment_ptr assignment(new randomColoring(domainSize,vectorLen));
    
    //construct permutations
    BREXWitness::permutation_ptr permutation(new randPermutation(domainSize));
    
    //construct witness
    BREXWitness witness(move(assignment), move(permutation));

    /** construct instance **/

    //construct constraints
    BREXInstance::constraintsPtr_t constraintsAssignment(new settingSaticfyingSystem(contextField,vectorLen*2,domainSize-1,witness, plusOneSequence(),true));
    BREXInstance::constraintsPtr_t constraintsPermutation(new settingSaticfyingSystem(contextField,vectorLen*2,domainSize,witness,witness.permutation(),true));
    
    //construct the instance
    BREXInstance instance(contextField,vectorLen,domainSizeIndicator,move(constraintsAssignment), move(constraintsPermutation), BREXInstance::boundaryConstraints_t(), getRandomPadding(vectorLen));	
    
    /** Return the result **/
    return pair<BREXInstance,BREXWitness>(move(instance),move(witness));
}

/***************************************************************************
 *
 *                              GTEST tests
 *
 ***************************************************************************/

/**
 * @brief   GTEST function to test completeness of BREXWitnessChecker::verify_constraints
 */
TEST(BREXWitnessChecker,verify_constraints_completeness){
    pair<BREXInstance,BREXWitness> validPair = generate_valid_constraints();
    EXPECT_TRUE(BREXWitnessChecker::verify_constraintsAssignment(validPair.first,validPair.second));
    EXPECT_TRUE(BREXWitnessChecker::verify_constraintsPermutation(validPair.first,validPair.second));
}

/**
 * @brief   GTEST function to test soundness of BREXWitnessChecker::verify_constraints
 * the Assignment constraints should not feet
 */
TEST(BREXWitnessChecker,verify_constraints_Assignment_soundness){
    pair<BREXInstance,BREXWitness> validPair = generate_invalid_constraints_Assignment();
    EXPECT_FALSE(BREXWitnessChecker::verify_constraintsAssignment(validPair.first,validPair.second));
}

/**
 * @brief   GTEST function to test soundness of BREXWitnessChecker::verify_constraints
 * the Permutation constraints should not feet
 */
TEST(BREXWitnessChecker,verify_constraints_Permutation_soundness){
    pair<BREXInstance,BREXWitness> validPair = generate_invalid_constraints_Permutation();
    EXPECT_FALSE(BREXWitnessChecker::verify_constraintsPermutation(validPair.first,validPair.second));
}

/**
 * @brief   GTEST function to test soundness of BREXWitnessChecker::verify_constraints
 * the both Assignment & Permutation constraints should not feet
 */
TEST(BREXWitnessChecker,verify_constraints_Assignment_N_Permutation_soundness){
    pair<BREXInstance,BREXWitness> validPair = generate_invalid_constraints_both();
    EXPECT_FALSE(BREXWitnessChecker::verify_constraintsAssignment(validPair.first,validPair.second));
    EXPECT_FALSE(BREXWitnessChecker::verify_constraintsPermutation(validPair.first,validPair.second));
}

/**
 * @brief   GTEST function to test completeness of BREXWitnessChecker::verify_permutations
 */
TEST(BREXWitnessChecker,verify_permutations_completeness){
    pair<BREXInstance,BREXWitness> validPair = generate_valid_permutations();
    EXPECT_TRUE(BREXWitnessChecker::verify_permutation(validPair.first,validPair.second));
}

/**
 * @brief   GTEST function to test soundness of BREXWitnessChecker::verify_permutations
 */
TEST(BREXWitnessChecker,verify_permutations_soundneness){
    pair<BREXInstance,BREXWitness> validPair = generate_invalid_permutations();
    EXPECT_FALSE(BREXWitnessChecker::verify_permutation(validPair.first,validPair.second));
}

/**
 * @brief   GTEST function to test completeness of BREXWitnessChecker::verify_boundary
 */
TEST(BREXWitnessChecker,verify_boundary_completeness){
    pair<BREXInstance,BREXWitness> validPair = generate_valid_boundary();
    EXPECT_TRUE(BREXWitnessChecker::verify_boundary(validPair.first,validPair.second));
}

/**
 * @brief   GTEST function to test soundness of BREXWitnessChecker::verify_boundary
 */
TEST(BREXWitnessChecker,verify_boundary_soundness){
    pair<BREXInstance,BREXWitness> invalidPair = generate_invalid_boundary();
    EXPECT_FALSE(BREXWitnessChecker::verify_boundary(invalidPair.first,invalidPair.second));
}

} //PCP_UTEST namespace
