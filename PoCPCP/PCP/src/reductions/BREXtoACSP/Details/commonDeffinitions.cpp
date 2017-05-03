#include "commonDeffinitions.hpp"
#include <algorithm>

using std::vector;
using std::ceil;
using std::max;
using std::pair;
using std::count;
using std::move;
using std::unique_ptr;
using Infrastructure::Log2;
using Infrastructure::POW2;
using Algebra::FieldElement;
using Algebra::PolynomialDegree;
using Algebra::PolynomialInterface;
using Algebra::UnivariatePolynomialGeneral;
using Algebra::getStandartBasis;
using Algebra::zero;

namespace PCP_Project{
namespace BREXtoACSP{

	commonDeffinitions::commonDeffinitions(const BREXInstance& instance) : contextField_(instance.contextField()),witnessDegreeBound_(0){
    
    //getting variables partition
    {
    variablesPerm_ = getRoutedIndexes(instance);
    auto vars_n_constraints = getUnroutedIndexes_and_ConstraintsChi(instance);
    variablesNonPerm_ = vars_n_constraints.first;
    constraintsChi_ = move(vars_n_constraints.second);
    
    //initialize constraints Pi
    for(const auto& c : instance.constraintsPermutation().constraints()){
        constraintsPi_.push_back(c->clone());
    }
    }

    //get the variables location translation
    {
    brexVarLocation_.resize(instance.vectorsLen());

    for(size_t i=0; i<variablesPerm_.size() ; i++){
        brexVarLocation_[variablesPerm_[i]] = {true,i};
    }

    for(size_t i=0; i<variablesNonPerm_.size() ; i++){
        brexVarLocation_[variablesNonPerm_[i]] = {false,i};
    }
    }

    //getting dimensions and sizes
    heightSpaceDimension_ = instance.domainSizeIndicator();
    
    widthSpaceDimension_ = 0;
    if(variablesPerm_.size() > 0){
        widthSpaceDimension_ = ceil(Log2(instance.domainSizeIndicator()+2));
    }
    
    amountOfPermutationLayers_ = 0;
    if(variablesPerm_.size() > 0){
        amountOfPermutationLayers_ = 2*(1+variablesPerm_.size());
    }
    
    {
    const double spaceWidth = POW2(widthSpaceDimension_);
    const double totalConstraintsAmount = instance.constraintsAssignment().constraints().size() + instance.constraintsPermutation().constraints().size();
    const size_t vanishingLayersConstraints = ceil(totalConstraintsAmount/spaceWidth);
    vanishingLayersSpaceDimension_ = ceil(Log2(amountOfPermutationLayers_ + vanishingLayersConstraints));
    }

    {
    const size_t spaceWidth = POW2(widthSpaceDimension_);
    
    size_t acspWitnessSizeUsedForPermutation=0;
    if(variablesPerm_.size() > 0){
        acspWitnessSizeUsedForPermutation = 2*(instance.domainSize()+1)*(variablesPerm_.size()+1)*spaceWidth;
    }

    const size_t acspWitnessSizeUsedForAssignment = POW2(instance.domainSizeIndicator())*variablesNonPerm_.size();
    const size_t acspWitnessDegreeBound = acspWitnessSizeUsedForPermutation + acspWitnessSizeUsedForAssignment;
    const size_t acspWitnessDegBound_roundedToPowOf2 = POW2(ceil(Log2(acspWitnessDegreeBound)));
    witnessDegreeBound_ = PolynomialDegree(acspWitnessDegBound_roundedToPowOf2-1);
    
    const size_t vanishingSpaceSize = POW2(heightSpaceDimension_ + widthSpaceDimension_ + vanishingLayersSpaceDimension_);//ARIEL: Should this not include the witness? Michael: No, the VANISHING space is for instance only, and is independent of the witness.
	const PolynomialDegree maxTestDeg = getMaxTestDegree(instance, variablesPerm_.size() > 0, witnessDegreeBound_);

	//verify the field is big enough for reduction
	{
		const auto compositionPolyDegreeBound = PolynomialDegree::integral_t(degreeOfProduct(maxTestDeg,PolynomialDegree(vanishingSpaceSize)));
		const short soundnessFactorParameter = 4; //This is a magic number for Michael.
		const size_t minimalDegreeFromReduction = ceil(Log2(soundnessFactorParameter) + Log2(compositionPolyDegreeBound));
		const size_t brexContextFieldDegree = instance.contextField().degree();
		_COMMON_ASSERT(brexContextFieldDegree >= minimalDegreeFromReduction, "BREX context field degree is too small, can't reduce to ACSP");
	}
    }
}

PolynomialDegree commonDeffinitions::getMaxTestDegree(const BREXInstance& instance, const bool hasRoutingNetwork, const PolynomialDegree& ascpWitnessDegreeBound){
    
    // initialize a vector of zeros as input degrees,
    // for checking the maximal degree of a polynomial
    // in each constraint system
    const size_t assignmetnSize = 2*instance.vectorsLen();
    vector<PolynomialDegree> inputDegrees(assignmetnSize,ascpWitnessDegreeBound);
   
    // find maximal degree in constraint system tests
    PolynomialDegree maxDeg = PolynomialDegree::getZeroPolyDegree();

    for(const auto& p : instance.constraintsAssignment().constraints()){
        maxDeg = max(maxDeg, p->getDegreeBound(inputDegrees));
    }
    
    for(const auto& p : instance.constraintsPermutation().constraints()){
        maxDeg = max(maxDeg, p->getDegreeBound(inputDegrees));
    }

    // find maximal in case there is a routing network
    // the maximal degree test for routing is 2 (specifically booleanity testing)
    if(hasRoutingNetwork){
        maxDeg = max(maxDeg,PolynomialDegree(2));
    }

    // return maximal degree
    return maxDeg;
}


// Routing
vector<size_t> commonDeffinitions::getRoutedIndexes(const BREXInstance& instance){
    const auto vecLen = instance.vectorsLen();
    const auto& CSystemPermutation = instance.constraintsPermutation();

    vector<size_t> usedIndexes;
    //Push only the indexes that must be routed
    //meaning, only those that the constraints system uses
    //from the second vector
    for(size_t i=vecLen; i< 2*vecLen; i++){
        if(CSystemPermutation.varUsed(i)){
            usedIndexes.push_back(i-vecLen);
        }
    }

    return usedIndexes;
}

pair<vector<size_t>,vector<unique_ptr<PolynomialInterface>>> commonDeffinitions::getUnroutedIndexes_and_ConstraintsChi(const BREXInstance& instance){
    const auto vecLen = instance.vectorsLen();
    const auto& CSystemPermutation = instance.constraintsPermutation();

    vector<size_t> unusedIndexes;
    //Push only the indexes that are not routed
    //meaning, only those that the constraints system not uses
    //from the second vector
    for(size_t i=vecLen; i< 2*vecLen; i++){
        if(!CSystemPermutation.varUsed(i)){
            unusedIndexes.push_back(i-vecLen);
        }
    }

    //if there are any univariate constraints,
    //reorder the variable, so those used by the most
    //common univariate constraints are all at the beginning
    vector<size_t> reorderedUnusedIndexes;
    vector<unique_ptr<PolynomialInterface>> reorderedConstraints;
    try{
        const auto uniPolyInfo = getCommonUnivariateConstraint(instance);
        
        //push first the indexes that are unrouted
        //and used by most common univariate constraint
        vector<PolynomialInterface*> commonUnivariateConstraints;
        for(const auto p_inf : uniPolyInfo.second){
            const size_t varId = p_inf.second;
            if((count(unusedIndexes.begin(),unusedIndexes.end(),varId) > 0) && (count(reorderedUnusedIndexes.begin(),reorderedUnusedIndexes.end(),varId) == 0)){
                reorderedUnusedIndexes.push_back(varId);
                reorderedConstraints.push_back(p_inf.first->clone());
                commonUnivariateConstraints.push_back(p_inf.first);
            }
        }

        //now push the rest of unused variables
        for(const auto& v : unusedIndexes){
            if(count(reorderedUnusedIndexes.begin(),reorderedUnusedIndexes.end(),v) == 0){
                reorderedUnusedIndexes.push_back(v);
            }
        }

        //now push all the rest of the constraints
        for(const auto& c : instance.constraintsAssignment().constraints()){
            if(count(commonUnivariateConstraints.begin(),commonUnivariateConstraints.end(),c.get()) == 0){
                reorderedConstraints.push_back(c->clone());
            }
        }
    }
    catch(...){
        reorderedUnusedIndexes = unusedIndexes;
        reorderedConstraints.resize(0);
        for(const auto& c : instance.constraintsAssignment().constraints()){
            reorderedConstraints.push_back(c->clone());
        }
    }

    return pair<vector<size_t>,vector<unique_ptr<PolynomialInterface>>>(reorderedUnusedIndexes,move(reorderedConstraints));
}

//return the most common univariate constraint
//and a vector where each node is a constraint pointer that is equivalent to it,
//and the variable it uses
//If no univariate constraint is found and exception is raised
pair<UnivariatePolynomialGeneral,vector<pair<PolynomialInterface*,size_t> > > commonDeffinitions::getCommonUnivariateConstraint(const BREXInstance& instance){
    vector<pair<UnivariatePolynomialGeneral,vector<pair<PolynomialInterface*,size_t>>>> univariatePolys; 
    const size_t totalVars = instance.vectorsLen();
    const vector<PolynomialDegree> inputDegrees(totalVars*2, PolynomialDegree(1));
    
    //collect all univariate constraints and their info
    for(const auto& c : instance.constraintsAssignment().constraints()){
        size_t numVarsUsed = 0;
        size_t lastVarUsed = 0;
        for( size_t i=0; i< totalVars*2; i++){
            if(c->isEffectiveInput(i)){
                numVarsUsed++;
                lastVarUsed = i;
            }
        }
        if((numVarsUsed == 1) && (lastVarUsed < totalVars)){
            const auto degreeBound = c->getDegreeBound(inputDegrees);
            const size_t interpolationBasisSize = ceil(Log2(PolynomialDegree::integral_t(degreeBound)+1));
            const auto interpolationBasis = getStandartBasis(interpolationBasisSize);
            const vector<FieldElement> orderedBasis(interpolationBasis.begin(), interpolationBasis.end());

            //construct the evaluation
            vector<FieldElement> evaluation(POW2(interpolationBasisSize));
            {
                vector<FieldElement> assignment(totalVars*2);
                for(size_t i=0; i< evaluation.size() ; i++){
                    assignment[lastVarUsed] = getSpaceElementByIndex(orderedBasis,zero(),i);
                    evaluation[i] = c->eval(assignment);
                }
                const UnivariatePolynomialGeneral poly(evaluation,orderedBasis,zero());
                //add poly to the menaged set
                {
                    const pair<PolynomialInterface*,size_t> currPair(c.get(),lastVarUsed);
                    bool found = false;
                    for(auto& p : univariatePolys){
                        if(p.first == poly){
                            p.second.push_back(currPair);
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        const vector<pair<PolynomialInterface*,size_t>> usageVec(1,currPair);
                        const pair<UnivariatePolynomialGeneral,vector<pair<PolynomialInterface*,size_t>>> newPolyInfo(poly,usageVec);
                        univariatePolys.push_back(newPolyInfo);
                    }
                }
            }
        }
    }

    //if there are no univariate polys at all
    //rais an exception, the caller should handle it
    if(univariatePolys.size() == 0){
        throw "no univariate constraints found";
    }

    //find the most common used univariate constraint
    size_t mostUsed = 0;
    {
        size_t maxUse = 0;
        for (size_t i=0; i<univariatePolys.size(); i++){
            if (maxUse < univariatePolys[i].second.size()){
                maxUse = univariatePolys[i].second.size();
                mostUsed = i;
            }
        }
    }

    //return info about most common poly
    return univariatePolys[mostUsed];
}

} //namespace BREXtoACSP
} //namespace PCP_Project
