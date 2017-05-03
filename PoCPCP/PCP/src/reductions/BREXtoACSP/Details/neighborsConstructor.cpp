#include "neighborsConstructor.hpp"

#include <algorithm>

namespace PCP_Project{
namespace BREXtoACSP{

using Algebra::LinearPolynomial;
using Algebra::PolynomialInterface;
using Algebra::UnivariatePolynomialInterface;
using Algebra::zero;
using Algebra::one;
using Algebra::FieldElement;
using Algebra::mapIntegerToFieldElement;

using std::find;
using std::array;
using std::vector;
using std::map;
using std::unique_ptr;
using std::move;

typedef AcspNeighbors::polynomialIndicator_t polynomialIndicator_t;

AcspNeighbors::AcspNeighbors(const BREXInstance& brexInstance ,const common& commonDef, const instanceMappings& instanceMapping, const CS_testLocations& testLocations):
commonDef_(commonDef), 
instanceMapping_(instanceMapping),
testLocations_(testLocations),
permutationCS_usageVector_(getCsUsageVector(commonDef.getConstraintsPi(),brexInstance.constraintsPermutation().varsAmount())),
assignmentCS_usageVector_(getCsUsageVector(commonDef.getConstraintsChi(),brexInstance.constraintsAssignment().varsAmount()))
{
    //set the context field
    commonDef_.contextField().setContext(); 
    
    //Add the Identity neighbor
    locationOfId_ = addNeighbor(constructIdNeighbor());

    //Construct the Permutation CS neighbors
    constructPermCS();
    
    //Construct the Permutation CS neighbors
    constructAssignmentCS();

    //
    // Routing network related naighbors
    //
    
    if (commonDef.hasRoutingNetwork()){
        initRoutingNetworkNeighbors();
    }
}

void AcspNeighbors::initRoutingNetworkNeighbors(){
    
    //Add the TWIN LAYER neighbor
    locationOfTwinLayer_ = addNeighbor(constructTwinLayerNeighbor());

    //Add the DeBruijn network neighbors
    constructDeBruijn();

    //Add neighbors for accessing the routing bit from any layer of the network
    constructRoutingBitAccess();
}

vector<unique_ptr<const UnivariatePolynomialInterface> > AcspNeighbors::getNeighborPolynomials()const{
    vector<unique_ptr<const UnivariatePolynomialInterface> > res;
    
    for(const LinearPolynomial& neighbor : neighbors_){
        res.push_back(unique_ptr<const UnivariatePolynomialInterface>(new LinearPolynomial(neighbor)));
    }

    return move(res);
}

size_t AcspNeighbors::polynomialsNum()const{
    return neighbors_.size();
}

size_t AcspNeighbors::locationOfId()const{
    return retLocation(locationOfId_);
}

size_t AcspNeighbors::locationOfTwinLayer()const{
    return retLocation(locationOfTwinLayer_);
}

size_t AcspNeighbors::locationOfDeBruijn(const short dbNeighborId, const short affineCosetId ,const size_t layerId)const{
    return retLocation(locationOfDeBruijn_[dbNeighborId][affineCosetId][layerId]);
}

size_t AcspNeighbors::locationOfRoutingBit(const size_t layerId)const{
    return retLocation(locationOfRoutingBit_[layerId]);
}

size_t AcspNeighbors::locationOfPermCS(polynomialIndicator_t poly, const size_t varId)const{
    return retLocation(locationOfPermCS_.at(poly)[varId]);
}

bool   AcspNeighbors::existsPermCS(polynomialIndicator_t poly, const size_t varId)const{
    return locationOfPermCS_.at(poly)[varId].exists;
}

size_t AcspNeighbors::locationOfAssignmentCS(polynomialIndicator_t poly, const size_t varId, const size_t neighborVersion)const{
    _COMMON_ASSERT((neighborVersion == 0)||(neighborVersion == 1),"neighbor version must be 0 or 1");
    return retLocation(locationOfAssignmentCS_.at(poly)[varId][neighborVersion]);
}

bool   AcspNeighbors::existsAssignmentCS(polynomialIndicator_t poly, const size_t varId)const{
    return locationOfAssignmentCS_.at(poly)[varId][0].exists;
}

struct AcspNeighbors::NeighborLocation_stc AcspNeighbors::addNeighbor(const LinearPolynomial& neighbor){
    
    //check if this neighbor already exists,
    //if so, return its index
    const auto currLocation =  find(neighbors_.begin(), neighbors_.end(), neighbor);
    if(currLocation != neighbors_.end()){
        return currLocation - neighbors_.begin();
    }

    //otherwise, add it to the vector
    //and return its index
    neighbors_.push_back(neighbor);
    return neighbors_.size() - 1;
}

Algebra::FieldElement AcspNeighbors::getGenerator()const{
    using NTL::GF2X;
    using NTL::SetCoeff;
    using NTL::to_GF2E;

    GF2X generator;
    SetCoeff(generator,1);

    return (FieldElement)to_GF2E(generator % commonDef_.rowsModulus());
}

size_t AcspNeighbors::retLocation(const struct NeighborLocation_stc& loc){
    if(!loc.exists) _COMMON_FATAL("Neighbor does not exist");
    return loc.index;
}

LinearPolynomial AcspNeighbors::constructIdNeighbor(){
    return LinearPolynomial(zero(), one());
}

LinearPolynomial AcspNeighbors::constructTwinLayerNeighbor()const{
    const FieldElement STEP = instanceMapping_.mapPermutationLayerId_spaceElement(1)-instanceMapping_.mapPermutationLayerId_spaceElement(0);
    return LinearPolynomial(STEP,one());
}

LinearPolynomial AcspNeighbors::constructLinearDeBruijn_N0()const{
    const FieldElement generator = getGenerator();
    return LinearPolynomial(zero(),generator);
}

LinearPolynomial AcspNeighbors::constructLinearDeBruijn_N1()const{
    const LinearPolynomial bitXor(one(),one());
    const LinearPolynomial N0 = constructLinearDeBruijn_N0();
    
    //because we use reverse DeBruijn network, the bit xor
    //should be done before moving to "Neighbor 0"
    return N0.compose(bitXor);
}

LinearPolynomial AcspNeighbors::applyShiftedLinearOperation(const size_t layer_id, const LinearPolynomial operation)const{
    const LinearPolynomial moveToLayer(instanceMapping_.mapPermutationLayerId_spaceElement(layer_id),one());
    return moveToLayer.compose(operation.compose(moveToLayer));
}

FieldElement AcspNeighbors::DeBruijn_fixRowIdCarry()const{
    const size_t numBits = commonDef_.heightSpaceDimension();
    const size_t shift = 0;
    size_t fixAsInt = 1 + (1<<numBits);

    const FieldElement constShift = mapIntegerToFieldElement(shift,numBits+1,fixAsInt);

    return constShift;
}

FieldElement AcspNeighbors::DeBruijn_fixColumnIdCarry()const{
    const size_t shift = commonDef_.heightSpaceDimension();
    const FieldElement modulus(NTL::to_GF2E(commonDef_.columnsModulus()));

    const FieldElement constShift = modulus*mapIntegerToFieldElement(shift,1,1);

    return constShift;
}

FieldElement AcspNeighbors::DeBruijn_getFixByCoset(const short cosetId)const{
    switch(cosetId){
        case 0 : return zero();
        case 1 : return DeBruijn_fixColumnIdCarry();
        case 2 : return DeBruijn_fixRowIdCarry();
        case 3 : return DeBruijn_fixColumnIdCarry() + DeBruijn_fixRowIdCarry();
        default : _COMMON_FATAL("Coset ID in not supported");
    }
}

void AcspNeighbors::constructDeBruijn(){
    //some constants
    const size_t layersNum = commonDef_.variablesPerm().size()*2;

    const LinearPolynomial N0 = constructLinearDeBruijn_N0();
    const LinearPolynomial N1 = constructLinearDeBruijn_N1();

    for(size_t layerId=0; layerId <layersNum ; layerId++){
        
        const LinearPolynomial N0_shifted = applyShiftedLinearOperation(layerId, N0);
        const LinearPolynomial N1_shifted = applyShiftedLinearOperation(layerId, N1);

        for(short cosetId=0; cosetId<4; cosetId++){
            const FieldElement currFix = DeBruijn_getFixByCoset(cosetId);
            
            locationOfDeBruijn_[0][cosetId].push_back(addNeighbor(N0_shifted + currFix));
            locationOfDeBruijn_[1][cosetId].push_back(addNeighbor(N1_shifted + currFix));
        }
    }
}

void AcspNeighbors::constructRoutingBitAccess(){
    //some constants
    const size_t doubleLayersNum = commonDef_.variablesPerm().size();
    const size_t rBitsLayerId = doubleLayersNum;
    const FieldElement rBitsLocation = instanceMapping_.mapPermutationLayerId_spaceElement(2*doubleLayersNum);

    //initialization
    locationOfRoutingBit_.resize(2*doubleLayersNum);

    //fill neighbors
    for (size_t dLayerId=0; dLayerId < doubleLayersNum; dLayerId++){
        const FieldElement currLocation = instanceMapping_.mapPermutationLayerId_spaceElement(2*dLayerId);
        const FieldElement shift = rBitsLocation - currLocation;
        locationOfRoutingBit_[2*dLayerId] = addNeighbor(LinearPolynomial(shift,one()));

        //the neighbor for the twin layer is wxactly the same one,
        //because both the twin layers ids and the routing bit layers ids differ only
        //in the least segnificant bit
        locationOfRoutingBit_[2*dLayerId + 1] = locationOfRoutingBit_[2*dLayerId];
    }
}

void AcspNeighbors::constructPermCS(){
   
    for (const auto& currPair : permutationCS_usageVector_){
        const auto& currPoly = currPair.first;
        const auto& currUsageVector = currPair.second;
        const size_t testLocation = testLocations_.indexOfConstraint_Permuation(currPoly);
        const size_t varsNum = currUsageVector.size();
        locationOfPermCS_[currPoly].resize(varsNum);

        for(size_t varId=0; varId < varsNum; varId++){
            if(currUsageVector[varId] == true){
                locationOfPermCS_[currPoly][varId] = addNeighbor(constructPermCS(testLocation, varId));
            }
        }
    }
}

LinearPolynomial AcspNeighbors::constructPermCS(const size_t nonPermElemId, const size_t varId)const{
    const FieldElement testLocation = instanceMapping_.mapNonPermutationElement(nonPermElemId);
    const size_t varsInRow = commonDef_.variablesPerm().size() + commonDef_.variablesNonPerm().size();
    
    //if the variable is from current row get it from current row
    if(varId < varsInRow){
        return moveFromPointToVarId(testLocation,varId);
    }
    //else get it from next configuration
    else{
        const auto location = commonDef_.getVarLocation(varId - varsInRow);
        //it must be a permutation variable, because this is how permutation variables defined
        _COMMON_ASSERT(location.inPerm, "A variable not from current row, used in permutation constraint system must be a permutation variable");
        const FieldElement destination = instanceMapping_.mapPermutationElement_spaceElement(0,2*location.index + 1);

        const FieldElement delta = destination - testLocation;

        return LinearPolynomial(delta,one());
    }
}

void AcspNeighbors::constructAssignmentCS(){
    for (const auto& currPair : assignmentCS_usageVector_){
        const auto& currPoly = currPair.first;
        const auto& currUsageVector = currPair.second;
        const size_t testLocation = testLocations_.indexOfConstraint_Assignment(currPoly);
        const size_t varsNum = currUsageVector.size();
        locationOfAssignmentCS_[currPoly].resize(varsNum);

        for(size_t varId=0; varId < varsNum ;varId++){
            if(currUsageVector[varId] == true){
                locationOfAssignmentCS_[currPoly][varId][0] = addNeighbor(constructAssignmentCS(testLocation, varId, false));
                locationOfAssignmentCS_[currPoly][varId][1] = addNeighbor(constructAssignmentCS(testLocation, varId, true));
            }
        }
    }
}

LinearPolynomial AcspNeighbors::constructAssignmentCS(const size_t nonPermElemId, const size_t varId, const bool withCarry)const{
    using NTL::to_GF2E;
    
    const FieldElement testLocation = instanceMapping_.mapNonPermutationElement(nonPermElemId);
    const size_t varsInRow = commonDef_.variablesPerm().size() + commonDef_.variablesNonPerm().size();
    
    //if the variable is from current row get it from current row
    if(varId < varsInRow){
        return moveFromPointToVarId(testLocation,varId);
    }
    
    //else get it from next configuration
    else{
        const FieldElement nextRowOffset = (withCarry) ? FieldElement(to_GF2E(commonDef_.rowsModulus())): zero();

        const LinearPolynomial goToRowID_linearSpace(zero() - testLocation,one());
        const LinearPolynomial goToNextRow(nextRowOffset, getGenerator());
        const LinearPolynomial goToVarLocation = moveFromPointToVarId(zero(),varId - varsInRow);

        return goToVarLocation.compose(goToNextRow.compose(goToRowID_linearSpace));
    }
}

map<polynomialIndicator_t , vector<bool>> AcspNeighbors::getCsUsageVector(const vector<unique_ptr<PolynomialInterface>>& cs,const size_t varsAmount){
    map<polynomialIndicator_t , vector<bool>> usageVector;

    for(size_t poly_indx=0; poly_indx < cs.size(); poly_indx++){
        vector<bool> currPolyUsage;
        for(size_t varId=0; varId< varsAmount; varId++){
            if(cs[poly_indx]->isEffectiveInput(varId)){
                currPolyUsage.push_back(true);
            }
            else{
                currPolyUsage.push_back(false);
            }
        }
   
        //add curr poly usage vector to CS usage vector
        usageVector.emplace(poly_indx,currPolyUsage);
    }

    return usageVector;
}

LinearPolynomial AcspNeighbors::moveFromPointToVarId(const FieldElement src, const size_t varId)const{
    const FieldElement destination = instanceMapping_.mapVariable(varId);

    const FieldElement delta = destination - src;

    return LinearPolynomial(delta,one());
}

} //namespace BREXtoACSP
} //namespace PCP_Project
