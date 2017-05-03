#include "languages/TinyRAM/TinyRAMDefinitions.hpp"
#include "languages/BREX/BrexWitnessChecker_UTEST.hpp"
#include "languages/ACSP/ACSPWitnessChecker.hpp"
#include "witnessReduction.hpp"
#include "reductions/BREXtoACSP/BREXtoACSP.hpp"
#include <gtest/gtest.h>
namespace{

using PCP_Project::BREXtoACSP::common;
using PCP_Project::BREXtoACSP::commonMappings;
using PCP_Project::BREXtoACSP::witnessMappings;
using PCP_Project::BREXtoACSP::witnessReduction;
using PCP_Project::CBREXtoACSP;
using PCP_Project::BREXWitness;
using PCP_Project::BREXInstance;
using PCP_Project::ACSPWitness;
using PCP_Project::ACSPInstance;
using PCP_Project::ACSPWitnessChecker;
using PCP_Project::ConstraintSys;
using PCP_Project::Sequence;
using PCP_Project::initTinyRAMParams;
using Algebra::FiniteField;
using Algebra::FieldElement;
using Algebra::PolynomialInterface;
using Algebra::UnivariatePolynomialInterface;
using Algebra::UnivariatePolynomialGeneral;
using Algebra::PolynomialDegree;
using Algebra::zero;
using Algebra::one;
using Algebra::generateRandom;
using Algebra::elementsSet_t;
using Algebra::mapFieldElementToInteger;
using Infrastructure::POW2;
using Infrastructure::Log2;
using Infrastructure::CEIL;
using std::pair;
using std::rand;
using std::vector;
using std::unique_ptr;
using std::move;
using std::max;
using std::vector;

typedef pair<BREXInstance,BREXWitness> BrexPair;

/*************************************************************
 *                   TESTS IMPLEMENTATION
 ************************************************************/

class BREXtoACSP_tester : private CBREXtoACSP, private witnessReduction{
    public:

    /*******************************************************
     *           ACSP Witness modifiers
     *     (for soundness tests of the constraints poly)
     *******************************************************/
     static ACSPWitness ruinRoutingBit(const BrexPair& brex_pair){
   
        //get common information
        common commonDef(brex_pair.first);
        witnessMappings witnessMapping(commonDef);
    
        //ordered basis for witness space
        const auto& basis = witnessMapping.getImageSpaceOrderedBasis();
        const vector<FieldElement> orderedBasis(basis.begin(),basis.end());
  
        const auto& partialInstance = brex_pair.first;
        
        //get the mapping
        auto mapping = witnessReduction::getEmbeddingMapping(partialInstance,brex_pair.second, commonDef,witnessMapping);

        //find some field element not in GF2
        FieldElement newVal = generateRandom();
        while((newVal == zero()) || (newVal == one())){
            newVal = generateRandom();
        }

        //find some victim index
        const size_t column = rand()%(POW2(commonDef.widthSpaceDimension())-2);
        const size_t row = rand()%POW2(commonDef.heightSpaceDimension());
        const size_t layer = rand()%2;

        //get the victim indicator
        const FieldElement victim = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkRoutingBit_spaceIndex(row,column,layer));
        const size_t victim_index = Algebra::getSpaceIndexOfElement(orderedBasis,Algebra::zero(),victim);

        //change the value
        mapping[victim_index] = newVal;

        //interpolate and create the witness
        unique_ptr<const ACSPWitness::polynomial> Assignment_ptr(new UnivariatePolynomialGeneral(mapping,orderedBasis,zero()));
        return ACSPWitness(move(Assignment_ptr));
     }
    
     static ACSPWitness ruinDeBruijnLastColumn(const BrexPair& brex_pair){
   
        //get common information
        common commonDef(brex_pair.first);
        witnessMappings witnessMapping(commonDef);
    
        //ordered basis for witness space
        const auto& basis = witnessMapping.getImageSpaceOrderedBasis();
        const vector<FieldElement> orderedBasis(basis.begin(),basis.end());
  
        const auto& partialInstance = brex_pair.first;
        
        //get the mapping
        auto mapping = witnessReduction::getEmbeddingMapping(partialInstance,brex_pair.second, commonDef,witnessMapping);

        //find some victim index
        const size_t lastCol = commonDef.imageWidth() - 2;
        const size_t row = rand()%POW2(commonDef.heightSpaceDimension());
        const size_t layer = rand()%(2*commonDef.variablesPerm().size());

        //get the victim indicator
        const FieldElement victim = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkElement_spaceIndex(row,lastCol,layer));
        const size_t victim_index = Algebra::getSpaceIndexOfElement(orderedBasis,Algebra::zero(),victim);

        //change the value
        //adding 1 makes it different for sure
        mapping[victim_index] += one();

        //interpolate and create the witness
        unique_ptr<const ACSPWitness::polynomial> Assignment_ptr(new UnivariatePolynomialGeneral(mapping,orderedBasis,zero()));
        return ACSPWitness(move(Assignment_ptr));
     }

     static ACSPWitness ruinAdditionalElementRouting(const BrexPair& brex_pair){
   
        //get common information
        common commonDef(brex_pair.first);
        witnessMappings witnessMapping(commonDef);
  
        //ordered basis for witness space
        const auto& basis = witnessMapping.getImageSpaceOrderedBasis();
        const vector<FieldElement> orderedBasis(basis.begin(),basis.end());
  
        const auto& partialInstance = brex_pair.first;
        
        //get the mapping
        auto mapping = witnessReduction::getEmbeddingMapping(partialInstance,brex_pair.second, commonDef,witnessMapping);

        //find some victim index
        const size_t firstCol = 0;
        const size_t layer = rand()%(2*commonDef.variablesPerm().size());

        //get the victim indicator
        const FieldElement victim = witnessMapping.mapPermutationColumnId_spaceElement(firstCol) + witnessMapping.mapPermutationLayerId_spaceElement(layer);
        const size_t victim_index = Algebra::getSpaceIndexOfElement(orderedBasis,Algebra::zero(),victim);

        //change the value
        //adding 1 makes it different for sure
        mapping[victim_index] += one();

        //interpolate and create the witness
        unique_ptr<const ACSPWitness::polynomial> Assignment_ptr(new UnivariatePolynomialGeneral(mapping,orderedBasis,zero()));
        return ACSPWitness(move(Assignment_ptr));
     }

     static ACSPWitness ruinDeBruijnVertexData(const BrexPair& brex_pair){
   
        //get common information
        common commonDef(brex_pair.first);
        witnessMappings witnessMapping(commonDef);
  
        //ordered basis for witness space
        const auto& basis = witnessMapping.getImageSpaceOrderedBasis();
        const vector<FieldElement> orderedBasis(basis.begin(),basis.end());
  
        const auto& partialInstance = brex_pair.first;
        
        //get the mapping
        auto mapping = witnessReduction::getEmbeddingMapping(partialInstance,brex_pair.second, commonDef,witnessMapping);

        //find some victim index
        const size_t lastCol = commonDef.imageWidth() - 2;
        const size_t column = rand()%lastCol;
        const size_t row = rand()%POW2(commonDef.heightSpaceDimension());
        const size_t layer = rand()%(2*commonDef.variablesPerm().size());

        //get the victim indicator
        const FieldElement victim = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkElement_spaceIndex(row,column,layer));
        const size_t victim_index = Algebra::getSpaceIndexOfElement(orderedBasis,Algebra::zero(),victim);

        //change the value
        //adding 1 makes it different for sure
        mapping[victim_index] += one();

        //interpolate and create the witness
        unique_ptr<const ACSPWitness::polynomial> Assignment_ptr(new UnivariatePolynomialGeneral(mapping,orderedBasis,zero()));
        return ACSPWitness(move(Assignment_ptr));
     }

};


/******************************************
 *     Constraints polynomial tests
 ******************************************/

TEST(BREXtoACSP_Constraints,reductionCompletness){
	const BrexPair src = PCP_UTESTS::generate_valid_pair();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = CBREXtoACSP::reduceWitness(src.first,src.second);
    
    EXPECT_TRUE(ACSPWitnessChecker::verify_vanishing(*instance,*witness));
}

TEST(BREXtoACSP_Constraints,reductionSoundness_routingBits){
	const BrexPair src = PCP_UTESTS::generate_valid_pair();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = BREXtoACSP_tester::ruinRoutingBit(src);
    
    EXPECT_FALSE(ACSPWitnessChecker::verify_vanishing(*instance,witness));
}

TEST(BREXtoACSP_Constraints,reductionSoundness_DeBruijnLastColumn){
	const BrexPair src = PCP_UTESTS::generate_valid_pair();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = BREXtoACSP_tester::ruinDeBruijnLastColumn(src);
    
    EXPECT_FALSE(ACSPWitnessChecker::verify_vanishing(*instance,witness));
}

TEST(BREXtoACSP_Constraints,reductionSoundness_DeBruijnLastAdditionalElementRouting){
	const BrexPair src = PCP_UTESTS::generate_valid_pair();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = BREXtoACSP_tester::ruinAdditionalElementRouting(src);
    
    EXPECT_FALSE(ACSPWitnessChecker::verify_vanishing(*instance,witness));
}

TEST(BREXtoACSP_Constraints,reductionSoundness_DeBruijnDataRouting){
	const BrexPair src = PCP_UTESTS::generate_valid_pair();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = BREXtoACSP_tester::ruinDeBruijnVertexData(src);
    
    EXPECT_FALSE(ACSPWitnessChecker::verify_vanishing(*instance,witness));
}

TEST(BREXtoACSP_Constraints,reductionComplitness_Constraints){
	const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = CBREXtoACSP::reduceWitness(src.first,src.second);
    
    EXPECT_TRUE(ACSPWitnessChecker::verify_vanishing(*instance,*witness));
}

TEST(BREXtoACSP_Constraints,reductionSoundness_Permutation_Constraints){
	const BrexPair src = PCP_UTESTS::generate_invalid_constraints_Permutation();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = CBREXtoACSP::reduceWitness(src.first,src.second);
    
    EXPECT_FALSE(ACSPWitnessChecker::verify_vanishing(*instance,*witness));
}

TEST(BREXtoACSP_Constraints,reductionSoundness_Assignment_Constraints){
	const BrexPair src = PCP_UTESTS::generate_invalid_constraints_Assignment();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = CBREXtoACSP::reduceWitness(src.first,src.second);
    
    EXPECT_FALSE(ACSPWitnessChecker::verify_vanishing(*instance,*witness));
}

}
