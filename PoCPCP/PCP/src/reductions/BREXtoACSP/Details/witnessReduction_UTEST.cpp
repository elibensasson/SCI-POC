#include "witnessReduction.hpp"
#include "../Routing/LongSymmetricDeBruijnNetwork.hpp"

#include "languages/BREX/BrexWitnessChecker_UTEST.hpp"

#include <gtest/gtest.h>

using PCP_Project::BREXtoACSP::witnessReduction;
using PCP_Project::BREXtoACSP::common;
using PCP_Project::BREXtoACSP::commonMappings;
using PCP_Project::BREXtoACSP::witnessMappings;
using PCP_Project::BREXtoACSP::LongSymmetricDeBruijnNetwork;
using PCP_Project::BREXWitness;
using PCP_Project::BREXInstance;
using PCP_Project::ACSPWitness;

using Infrastructure::POW2;

using Algebra::FieldElement;
using Algebra::zero;
using Algebra::one;

using std::pair;
using std::vector;

namespace{    

typedef pair<BREXInstance,BREXWitness> BrexPair;

FieldElement advanceToNextRow(const FieldElement src, const common& commonInfo){
    using NTL::GF2X;
    using NTL::GF2E;
    using NTL::power;

    //initialize the polynomial containing current location
    const GF2X srcLocationPoly = ((GF2E)src).LoopHole();
    
    //initialize the identity polynomial (aka "x")
    const GF2X x(1,1);

    //initialize a polynomial for the row location
    const GF2X srcRow = srcLocationPoly % power(x,commonInfo.heightSpaceDimension());

    //initialize the next row location polynomial
    const GF2X dstRow = (srcRow*x) % commonInfo.rowsModulus();

    //the destination location is exactly swapping the row location
    //polynomials in the row location area
    const GF2X dstLocationPoly = (srcLocationPoly - srcRow) + dstRow;

    //return result
    return FieldElement(to_GF2E(dstLocationPoly));
}

void verifyIntegrity(const BrexPair& brex_pair, const ACSPWitness& acsp_witness){

    //get common information
    const common commonDef(brex_pair.first);
    const witnessMappings witnessMapping(commonDef);
    const commonMappings& commonMapping = witnessMapping;

    //get the ACSP witness
    const auto& acspPoly = acsp_witness.assignmentPoly(); 
    
    //get the BREX pair
    const BREXInstance& instance = brex_pair.first;
    const BREXWitness& witness = brex_pair.second;
    
    // BREX permutation
    const auto& permutation = witness.permutation();

    //get amount of rows
    const size_t cyclicDomainSize = instance.domainSize();
    
    //get variables partition
    const vector<size_t> unroutedVars = commonDef.variablesNonPerm();
    const vector<size_t> routedVars = commonDef.variablesPerm();

    {
        //verify elements for circle
        size_t vecId = 0;
        FieldElement indicator = witnessMapping.mapNonPermutationElement(vecId,0) - commonMapping.mapNonPermutationElement(0);
        for(int i=0; i<cyclicDomainSize; i++){
            const auto& coloring = witness.get_color(vecId);
            const auto& permColoring = witness.get_color(permutation.getElementByIndex(vecId));

            //verify assignment on variables not in permutation
            for(size_t varLocation = 0; varLocation < unroutedVars.size(); varLocation++){
                const FieldElement currDelta = commonMapping.mapNonPermutationElement(varLocation);
                const FieldElement expectedVal = coloring[unroutedVars[varLocation]];
                const FieldElement val = acspPoly.eval(indicator + currDelta);

                EXPECT_EQ(expectedVal,val);
            }

            //verify assignment & permutation on variables in permutation
            for(size_t varLocation = 0; varLocation < routedVars.size(); varLocation++){

                //verify assignment
                {
                const FieldElement currDelta = commonMapping.mapPermutationElement_spaceElement(0,2*varLocation);
                const FieldElement expectedVal = coloring[routedVars[varLocation]];
                const FieldElement val = acspPoly.eval(indicator + currDelta);

                EXPECT_EQ(expectedVal,val);
                }

                //verify permutation
                {
                const FieldElement currDelta = commonMapping.mapPermutationElement_spaceElement(0,2*varLocation+1);
                const FieldElement expectedVal = permColoring[routedVars[varLocation]];
                const FieldElement val = acspPoly.eval(indicator + currDelta);

                EXPECT_EQ(expectedVal,val);
                }
            }
            
            
            indicator = advanceToNextRow(indicator,commonDef);
            vecId = (vecId+1) % cyclicDomainSize;
        }
    }
}

size_t cyclicShift(const size_t src, const char bitsAmount){
    
    //get mask for cyclic shift
    const size_t mask = POW2(bitsAmount) - 1;

    return ((src<<1)&mask) ^ (((src<<1)&(~mask)) >> bitsAmount);


}

void verifyRoutingNetwork(const BrexPair& brex_pair, const ACSPWitness& acsp_witness){

    //get common information
    common commonDef(brex_pair.first);
    witnessMappings witnessMapping(commonDef);

    //get the ACSP witness
    const auto& acspPoly = acsp_witness.assignmentPoly(); 

    //get the BREX pair
    const BREXInstance& instance = brex_pair.first;
    const BREXWitness& witness = brex_pair.second;

    //get amount of layers
    const size_t layersAmount = commonDef.variablesPerm().size()*2;
    
    //get amount of rows
    const size_t amountOfRows = instance.domainSize()+1;
    
    //get columns amount
    const size_t columnsAmount = LongSymmetricDeBruijnNetwork(instance.domainSizeIndicator()).getWingWidth();

    //get index of last column
    const size_t lastColumn = columnsAmount-1;

    //bits for row Id
    const size_t bitsForRowId = commonDef.heightSpaceDimension();

    for(size_t layerId=0; layerId < layersAmount; layerId++){
        for(size_t rowId=0; rowId < amountOfRows; rowId++){
            for (size_t columnId=0; columnId < columnsAmount; columnId++){
                
                //get the location of current point
                const FieldElement currLocation = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkElement_spaceIndex(rowId,columnId,layerId));

                //get the data of current point
                const FieldElement currData = acspPoly.eval(currLocation);

                //if current column is not the last column, verify DeBruijn routing
                if(columnId != lastColumn){
                    
                    //get neighbors row ids
                    const size_t DB_neighbor0_rowId = cyclicShift(rowId   , bitsForRowId);
                    //because this is a reverse DeBruijn, the xor with "1" must be done before the cyclic shift
                    const size_t DB_neighbor1_rowId = cyclicShift(rowId^1 , bitsForRowId);

                    //get the value of the routing bit
                    const FieldElement routingBitLocation = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkRoutingBit_spaceIndex(rowId,columnId,layerId));
                    const FieldElement routingBitVal = acspPoly.eval(routingBitLocation);

                    //verify routing bit is boolean
                    EXPECT_TRUE((routingBitVal == zero()) || (routingBitVal == one()));

                    //verify the routing
                    if (routingBitVal ==  zero()){
                        const FieldElement neighborLocation = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkElement_spaceIndex(DB_neighbor0_rowId,columnId+1,layerId));
                        const FieldElement neighborVal = acspPoly.eval(neighborLocation);
                        EXPECT_EQ(currData,neighborVal);
                    }
                    if (routingBitVal ==  one()){
                        const FieldElement neighborLocation = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkElement_spaceIndex(DB_neighbor1_rowId,columnId+1,layerId));
                        const FieldElement neighborVal = acspPoly.eval(neighborLocation);
                        EXPECT_EQ(currData,neighborVal);
                    }
                    
                }

                //if current column is the last column, verify identity with the twin column
                if(columnId == lastColumn){
                    
                    //get the location of twin point
                    const FieldElement twinLocation = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkElement_spaceIndex(rowId,columnId,layerId^1));

                    //get the data of twin point
                    const FieldElement twinData = acspPoly.eval(twinLocation);

                    //verify
                    EXPECT_NE(currLocation, twinLocation);
                    EXPECT_EQ(currData,twinData);
                }
            }
        }
    }
}


void verifyPermutation_aditionalElement(const BrexPair& brex_pair, const ACSPWitness& acsp_witness){

    //get common information
    common commonDef(brex_pair.first);
    witnessMappings witnessMapping(commonDef);

    //get the ACSP witness
    const auto& acspPoly = acsp_witness.assignmentPoly(); 

    //get the BREX pair
    const BREXInstance& instance = brex_pair.first;
    const BREXWitness& witness = brex_pair.second;

    //get the variables in permutation
    const vector<size_t> routedVars = commonDef.variablesPerm();
    
    for(size_t varIndex=0; varIndex < routedVars.size(); varIndex++){
                
        //get the location of current point
        const FieldElement currLocation = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkElement_spaceIndex(0,0,2*varIndex));

        //get the data of current point
        const FieldElement currData = acspPoly.eval(currLocation);

        //get the location of twin point
        const FieldElement twinLocation = witnessMapping.map_spaceIndex_to_fieldElement(witnessMapping.mapNetworkElement_spaceIndex(0,0,2*varIndex+1));

        //get the data of twin point
        const FieldElement twinData = acspPoly.eval(twinLocation);

        //expected data
        const FieldElement expectedVal = instance.paddingPi()[routedVars[varIndex]];

        //verify
        EXPECT_NE(currLocation, twinLocation);
        EXPECT_EQ(currData,expectedVal);
        EXPECT_EQ(twinData,expectedVal);

    }
}


/******************************************
 *                UTESTS
 ******************************************/

TEST(BREXtoACSPWitness,integrityWithBREXWitness){
    const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    const auto witness = witnessReduction::reduceWitness(src.first,src.second);
    verifyIntegrity(src,*witness);
}

TEST(BREXtoACSPWitness,RoutingNetwork){
    const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    const auto witness = witnessReduction::reduceWitness(src.first,src.second);
    verifyRoutingNetwork(src,*witness);
}

TEST(BREXtoACSPWitness,PermutationOnAdditionalElement){
    const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    const auto witness = witnessReduction::reduceWitness(src.first,src.second);
    verifyPermutation_aditionalElement(src,*witness);
}

} //anonymus namespace
