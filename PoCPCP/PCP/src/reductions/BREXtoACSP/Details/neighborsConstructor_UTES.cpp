#include "neighborsConstructor.hpp"
#include "witnessMappings.hpp"
#include "languages/BREX/BrexWitnessChecker_UTEST.hpp"
#include "common/Utils/TaskReporting.hpp"

#include <gtest/gtest.h>

using Infrastructure::POW2;

using Algebra::FieldElement;
using Algebra::generateRandom;
using Algebra::mapFieldElementToInteger;
using Algebra::PolynomialInterface;

using PCP_Project::BREXWitness;
using PCP_Project::BREXInstance;
using PCP_Project::ConstraintSys;
using PCP_Project::Task;

using PCP_Project::BREXtoACSP::common;
using PCP_Project::BREXtoACSP::instanceMappings;
using PCP_Project::BREXtoACSP::witnessMappings;
using PCP_Project::BREXtoACSP::AcspNeighbors;
using PCP_Project::BREXtoACSP::CS_testLocations;

using std::pair;


namespace{

typedef pair<BREXInstance,BREXWitness> BrexPair;

void verifyID(const BrexPair& brex_pair, const AcspNeighbors& neighbors){

    //get the neighbors polynomials vector
    const auto neighborsVec = neighbors.getNeighborPolynomials();

    //get ID poly index
    const auto idPolyIndex = neighbors.locationOfId();

    //verify is the ID by evaluating on a random element
    const FieldElement x = generateRandom();
    const FieldElement res = neighborsVec[idPolyIndex]->eval(x);
    EXPECT_EQ(x,res);
}

void verifyTwinLayer(const BrexPair& brex_pair, const AcspNeighbors& neighbors){
    
    const common defs(brex_pair.first);
    const witnessMappings mappings(defs);

    //get the neighbors polynomials vector
    const auto neighborsVec = neighbors.getNeighborPolynomials();

    //get TWIN LAYER poly index
    const auto twinLayerPolyIndex = neighbors.locationOfTwinLayer();

    //verify neighbor by testing on a rendom network location
    const size_t rowId = rand() % defs.imageHeight();
    const size_t columnId = rand() % defs.imageWidth();
    const size_t layerId = rand() % defs.amountOfPermutationLayers();
    const size_t layerId_twin = layerId ^ 1;
    
    const FieldElement currLocation = mappings.map_spaceIndex_to_fieldElement(mappings.mapNetworkElement_spaceIndex(rowId,columnId,layerId));
    const FieldElement res = neighborsVec[twinLayerPolyIndex]->eval(currLocation);

    const FieldElement ref = mappings.map_spaceIndex_to_fieldElement(mappings.mapNetworkElement_spaceIndex(rowId,columnId,layerId_twin));


    EXPECT_EQ(ref,res);
}

size_t cyclicShift(const size_t src, const char bitsAmount){
    
    //get mask for cyclic shift
    const size_t mask = POW2(bitsAmount) - 1;

    return ((src<<1)&mask) ^ (((src<<1)&(~mask)) >> bitsAmount);
}

void verifyDeBruijn(const BrexPair& brex_pair, const AcspNeighbors& neighbors){
    
    const common defs(brex_pair.first);
    const witnessMappings mappings(defs);

    //some constants
    const size_t bitsForRowId = defs.heightSpaceDimension();
    const size_t bitsForColumnId = defs.widthSpaceDimension();
    const size_t layersAmount = defs.variablesPerm().size()*2;
    const size_t amountOfRows = defs.imageHeight();
    const size_t columnsAmount = defs.imageWidth()-2;

    //get the neighbors polynomials vector
    const auto neighborsVec = neighbors.getNeighborPolynomials();

    //verify neighbors by testing on all network locations except
    //the last column which is not routed using DeBruijn neighbors
    for(size_t layerId=0; layerId < layersAmount; layerId++){
        for(size_t rowId=0; rowId < amountOfRows; rowId++){
            for (size_t columnId=0; columnId < columnsAmount; columnId++){

                //get the DeBruijn neighbors for curr location
                const size_t rowId_N0 = cyclicShift(rowId   , bitsForRowId);
                //because this is a reverse DeBruijn, the xor with "1" must be done before the cyclic shift
                const size_t rowId_N1 = cyclicShift(rowId^1 , bitsForRowId);

                //get the current location
                const FieldElement currLocation = mappings.map_spaceIndex_to_fieldElement(mappings.mapNetworkElement_spaceIndex(rowId,columnId,layerId));

                //calculate the coset Id, so we can
                //find the relevant DeBruijn neighbor
                const size_t rowIdLastBitLocation = bitsForRowId - 1;
                const size_t columnIdLastBitLocation = bitsForRowId + bitsForColumnId - 1;
                const size_t rowLastBit = mapFieldElementToInteger(rowIdLastBitLocation,1,currLocation);
                const size_t columnLastBit = mapFieldElementToInteger(columnIdLastBitLocation,1,currLocation);
                const short cosetId = 2*rowLastBit + columnLastBit;

                //get the neighbor polynomials
                const auto& N0 = neighborsVec[neighbors.locationOfDeBruijn(0,cosetId,layerId)];
                const auto& N1 = neighborsVec[neighbors.locationOfDeBruijn(1,cosetId,layerId)];

                //get expected locations
                const FieldElement ref0 = mappings.map_spaceIndex_to_fieldElement(mappings.mapNetworkElement_spaceIndex(rowId_N0,columnId+1,layerId));
                const FieldElement ref1 = mappings.map_spaceIndex_to_fieldElement(mappings.mapNetworkElement_spaceIndex(rowId_N1,columnId+1,layerId));

                //get results from neighbor polynomials
                const FieldElement res0 = N0->eval(currLocation);
                const FieldElement res1 = N1->eval(currLocation);

                EXPECT_EQ(ref0,res0);
                EXPECT_EQ(ref1,res1);

            }
        }
    }
}

void verifyRoutingBit(const BrexPair& brex_pair, const AcspNeighbors& neighbors){
    
    const common defs(brex_pair.first);
    const witnessMappings mappings(defs);

    //some constants
    const size_t bitsForRowId = defs.heightSpaceDimension();
    const size_t bitsForColumnId = defs.widthSpaceDimension();
    const size_t layersAmount = defs.variablesPerm().size()*2;
    const size_t amountOfRows = defs.imageHeight();
    const size_t columnsAmount = defs.imageWidth()-1;
    const size_t firstRoutingBitLayerId = layersAmount;

    //get the neighbors polynomials vector
    const auto neighborsVec = neighbors.getNeighborPolynomials();

    //verify neighbors by testing on all network locations
    for(size_t layerId=0; layerId < layersAmount; layerId++){
        
        const auto& routingBitAccessPoly = neighborsVec[neighbors.locationOfRoutingBit(layerId)];
        const size_t relevantRoutingBitsLayerId = firstRoutingBitLayerId + (layerId & 1);
        
        for(size_t rowId=0; rowId < amountOfRows; rowId++){
            for (size_t columnId=0; columnId < columnsAmount; columnId++){

                //get the current location
                const FieldElement currLocation = mappings.map_spaceIndex_to_fieldElement(mappings.mapNetworkElement_spaceIndex(rowId,columnId,layerId));

                //get expected locations
                const FieldElement ref = mappings.map_spaceIndex_to_fieldElement(mappings.mapNetworkElement_spaceIndex(rowId,columnId,relevantRoutingBitsLayerId));

                //get results from neighbor polynomials
                const FieldElement res = routingBitAccessPoly->eval(currLocation);

                EXPECT_EQ(ref,res);
            }
        }
    }
}

void verifyPermutationConstraints(const BrexPair& brex_pair, const AcspNeighbors& neighbors, const CS_testLocations& testLocations){
    
    const BREXInstance& partialInstance = brex_pair.first;
    const common defs(partialInstance);
    const witnessMappings mappings(defs);
    const instanceMappings& mappings_common(defs);

    //some constants
    const size_t rowsNum = defs.imageHeight();

    //get the neighbors polynomials vector
    const auto neighborsVec = neighbors.getNeighborPolynomials();

    //verify neighbors by testing on all permutation constraints, and all rows
    const ConstraintSys& permutationCS = partialInstance.constraintsPermutation();
    
    for(size_t currPoly_indx=0; currPoly_indx < defs.getConstraintsPi().size(); currPoly_indx++){
        const PolynomialInterface* currPoly = defs.getConstraintsPi()[currPoly_indx].get();
        const size_t testIndex = testLocations.indexOfConstraint_Permuation(currPoly_indx);

        //neighbors for current row
        for (size_t varId=0; varId < partialInstance.vectorsLen(); varId++){
            
            //verify the neighbor for current variable exists if and only if
            //it is used by the polynomial
            EXPECT_EQ(neighbors.existsPermCS(currPoly_indx,varId), currPoly->isEffectiveInput(varId));

            //if neighbor does not exist, nothing else to do here
            if(!neighbors.existsPermCS(currPoly_indx,varId))continue;

            //else, get the neighbor
            const auto& currNeighbor = neighborsVec[neighbors.locationOfPermCS(currPoly_indx,varId)];

            //verify the neighbor
            for(size_t rowId=0; rowId < rowsNum; rowId++){
                const FieldElement testLocation = mappings.mapNonPermutationElement(rowId, testIndex);
                const FieldElement rowIndicator = testLocation - mappings_common.mapNonPermutationElement(testIndex);
                const FieldElement res = currNeighbor->eval(testLocation);
                
                const auto varLocation = defs.getVarLocation(varId);
                if(varLocation.inPerm){
                    const FieldElement ref = rowIndicator + mappings.mapPermutationElement_spaceElement(0,2*varLocation.index);
                    EXPECT_EQ(ref,res);
                }
                else{
                    const FieldElement ref = mappings.mapNonPermutationElement(rowId,varLocation.index);
                    EXPECT_EQ(ref,res);
                }

                
            }
        }

        //neighbors for next row (next using witness permutation)
        const size_t shift = partialInstance.vectorsLen();
        for (size_t varId=shift; varId < shift+partialInstance.vectorsLen(); varId++){

            //verify the neighbor for current variable exists if and only if
            //it is used by the polynomial
            EXPECT_EQ(neighbors.existsPermCS(currPoly_indx,varId), currPoly->isEffectiveInput(varId));

            //if neighbor does not exist, nothing else to do here
            if(!neighbors.existsPermCS(currPoly_indx,varId))continue;

            //else, get the neighbor
            const auto& currNeighbor = neighborsVec[neighbors.locationOfPermCS(currPoly_indx,varId)];

            //verify the neighbor
            for(size_t rowId=0; rowId < rowsNum; rowId++){
                const FieldElement testLocation = mappings.mapNonPermutationElement(rowId, testIndex);
                const FieldElement rowIndicator = testLocation - mappings_common.mapNonPermutationElement(testIndex);
                const FieldElement res = currNeighbor->eval(testLocation);

                const auto varLocation = defs.getVarLocation(varId-shift);

                //A neighbor to a variable from next row in Permutation Constraints, must be to a permutation variable
                EXPECT_TRUE(varLocation.inPerm);

                const FieldElement ref = rowIndicator + mappings.mapPermutationElement_spaceElement(0,2*varLocation.index + 1);

                EXPECT_EQ(res,ref);

            }
        }
    }
}

void verifyAssignmentConstraints(const BrexPair& brex_pair, const AcspNeighbors& neighbors, const CS_testLocations& testLocations){
    
    const BREXInstance& partialInstance = brex_pair.first;
    const common defs(partialInstance);
    const witnessMappings mappings(defs);
    const instanceMappings& mappings_common(defs);

    //some constants
    const size_t rowsNum = defs.imageHeight();

    //get the neighbors polynomials vector
    const auto neighborsVec = neighbors.getNeighborPolynomials();

    //verify neighbors by testing on all permutation constraints, and all rows
    const ConstraintSys& assignmentCS = partialInstance.constraintsAssignment();
   
    const size_t rowIdMSB_index = defs.heightSpaceDimension()-1;
    
    for(size_t currPoly_indx=0; currPoly_indx < defs.getConstraintsChi().size(); currPoly_indx++){
        const PolynomialInterface* currPoly = defs.getConstraintsChi()[currPoly_indx].get();
        const size_t testIndex = testLocations.indexOfConstraint_Assignment(currPoly_indx);

        //neighbors for current row
        for (size_t varId=0; varId < partialInstance.vectorsLen(); varId++){
            
            //verify the neighbor for current variable exists if and only if
            //it is used by the polynomial
            EXPECT_EQ(neighbors.existsAssignmentCS(currPoly_indx,varId), currPoly->isEffectiveInput(varId));

            //if neighbor does not exist, nothing else to do here
            if(!neighbors.existsAssignmentCS(currPoly_indx,varId))continue;

            //else, get the neighbor
            
            //verify both versions are the same polynomial
            EXPECT_EQ(neighbors.locationOfAssignmentCS(currPoly_indx,varId,0),neighbors.locationOfAssignmentCS(currPoly_indx,varId,1));
            const auto& currNeighbor = neighborsVec[neighbors.locationOfAssignmentCS(currPoly_indx,varId,0)];

            //verify the neighbor
            for(size_t rowId=0; rowId < rowsNum; rowId++){
                const FieldElement testLocation = mappings.mapNonPermutationElement(rowId, testIndex);
                const FieldElement rowIndicator = testLocation - mappings_common.mapNonPermutationElement(testIndex);
                const FieldElement res = currNeighbor->eval(testLocation);
                
                const auto varLocation = defs.getVarLocation(varId);
                if(varLocation.inPerm){
                    const FieldElement ref = rowIndicator + mappings.mapPermutationElement_spaceElement(0,2*varLocation.index);
                    EXPECT_EQ(ref,res);
                }
                else{
                    const FieldElement ref = mappings.mapNonPermutationElement(rowId,varLocation.index);
                    EXPECT_EQ(ref,res);
                }

                
            }
        }
        //neighbors for next row (next using instance permutation)
        const size_t shift = partialInstance.vectorsLen();
        for (size_t varId=shift; varId < shift+partialInstance.vectorsLen(); varId++){

            //verify the neighbor for current variable exists if and only if
            //it is used by the polynomial
            EXPECT_EQ(neighbors.existsAssignmentCS(currPoly_indx,varId), currPoly->isEffectiveInput(varId));

            //if neighbor does not exist, nothing else to do here
            if(!neighbors.existsAssignmentCS(currPoly_indx,varId))continue;

            //else, get the neighbor

            //verify the neighbor
            for(size_t rowId=0; rowId < rowsNum; rowId++){
                const FieldElement testLocation = mappings.mapNonPermutationElement(rowId, testIndex);
                const FieldElement nextRowIndicator = mappings.mapNonPermutationElement(rowId+1,0) - mappings_common.mapNonPermutationElement(0);
                
                // get the MSB of the rowId in order to decide the right neighbor version
                const size_t neighborVersion = mapFieldElementToInteger(rowIdMSB_index,1,testLocation);
                const auto& currNeighbor = neighborsVec[neighbors.locationOfAssignmentCS(currPoly_indx,varId,neighborVersion)];
                
                 
                const FieldElement res = currNeighbor->eval(testLocation);

                const auto varLocation = defs.getVarLocation(varId-shift);
                if(varLocation.inPerm){
                    const FieldElement ref = nextRowIndicator + mappings.mapPermutationElement_spaceElement(0,2*varLocation.index);

                    switch(neighborVersion){
                        case 0: EXPECT_EQ(ref,res); break;
                        case 1: EXPECT_EQ(ref,res); break;
                        default : _COMMON_FATAL("Bad neighbor version");
                    }
                }
                else{
                    const FieldElement ref = mappings.mapNonPermutationElement(rowId+1,varLocation.index);

                    switch(neighborVersion){
                        case 0: EXPECT_EQ(ref,res); break;
                        case 1: EXPECT_EQ(ref,res); break;
                        default : _COMMON_FATAL("Bad neighbor version");
                    }
                }

            }
        }
    }
}


//
// GTESTS
//

TEST(BREXtoACSPNeighbors,verifyID){
    //initialize neighbors instance	
    const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    const common defs(src.first);
    const instanceMappings mappings(defs);
    const CS_testLocations testLocations(src.first);
    const AcspNeighbors neighbors(src.first,defs,mappings,testLocations);
    
    verifyID(src,neighbors);
}

TEST(BREXtoACSPNeighbors,verifyTwinLayer){

    //initialize neighbors instance	
    const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    const common defs(src.first);
    const instanceMappings mappings(defs);
    const CS_testLocations testLocations(src.first);
    const AcspNeighbors neighbors(src.first,defs,mappings,testLocations);

    //if there is no routing network, there is not expected a
    //TWIN LAYER naighbor
    if(!defs.hasRoutingNetwork())return;
    
    verifyTwinLayer(src,neighbors);
}

TEST(BREXtoACSPNeighbors,verifyDeBruijn){

    //initialize neighbors instance	
    const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    const common defs(src.first);
    const instanceMappings mappings(defs);
    const CS_testLocations testLocations(src.first);
    const AcspNeighbors neighbors(src.first,defs,mappings,testLocations);

    //if there is no routing network, there is not expected a
    //TWIN LAYER naighbor
    if(!defs.hasRoutingNetwork())return;
    
    verifyDeBruijn(src,neighbors);
}

TEST(BREXtoACSPNeighbors,verifyRoutingBit){

    //initialize neighbors instance	
    const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    const common defs(src.first);
    const instanceMappings mappings(defs);
    const CS_testLocations testLocations(src.first);
    const AcspNeighbors neighbors(src.first,defs,mappings,testLocations);

    //if there is no routing network, there is not expected a
    //TWIN LAYER naighbor
    if(!defs.hasRoutingNetwork())return;
    
    verifyRoutingBit(src,neighbors);
}

TEST(BREXtoACSPNeighbors,verifyPermutationConstraintsSystem){

    //initialize neighbors instance	
    const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    const common defs(src.first);
    const instanceMappings mappings(defs);
    const CS_testLocations testLocations(src.first);
    const AcspNeighbors neighbors(src.first,defs,mappings,testLocations);
    
    verifyPermutationConstraints(src,neighbors,testLocations);
}

TEST(BREXtoACSPNeighbors,verifyAssignmentConstraintsSystem){

    //initialize neighbors instance	
    const BrexPair src = PCP_UTESTS::generate_valid_constraints();
    const common defs(src.first);
    const instanceMappings mappings(defs);
    const CS_testLocations testLocations(src.first);
    const AcspNeighbors neighbors(src.first,defs,mappings,testLocations);
    
    verifyAssignmentConstraints(src,neighbors,testLocations);
}

}
