#include "witnessMappings.hpp"

namespace PCP_Project{
namespace BREXtoACSP{

using Algebra::FieldElement;
using Algebra::mapIntegerToFieldElement;
using Algebra::mapFieldElementToInteger;
using Algebra::elementsSet_t;
using Infrastructure::Log2;

witnessMappings::witnessMappings(const common& commonInfo):
    commonMappings(commonInfo),
    firstRoutingBitsLayer_(2*commonInfo.variablesPerm().size()),
    rowsModulus_(commonInfo.rowsModulus()),
    overflow_mask_(commonInfo.imageHeight()),
    imageSpaceDim_(calculateImageSpaceDim(commonInfo)){};

FieldElement witnessMappings::mapNonPermutationElement(const size_t vecId, const size_t varIndex)const{
    return map_x_power_modulu_poly(vecId, rowsModulus_) + commonMappings::mapNonPermutationElement(varIndex);
}

witnessMappings::spaceIndex_t witnessMappings::getNextRow_spaceIndex(const spaceIndex_t& row_spaceIndex)const{
    spaceIndex_t nextVal = row_spaceIndex<<1;
    if(nextVal & overflow_mask_){
        nextVal ^= rowsModulus_bin_;
    }
    return nextVal;
}

witnessMappings::spaceIndex_t witnessMappings::mapIndexOfNonPermutationVariable_spaceIndex(const spaceIndex_t& row_spaceIndex, const size_t& varIndex)const{
    return row_spaceIndex ^ mapNonPermutationVariable_spaceIndex(varIndex);
}

witnessMappings::spaceIndex_t witnessMappings::mapNetworkElement_spaceIndex(const size_t rowId, const size_t column, const size_t layerId)const{
    return rowId ^ mapPermutationElement_spaceIndex(column,layerId);
}

witnessMappings::spaceIndex_t witnessMappings::mapNetworkRoutingBit_spaceIndex(const size_t rowId, const size_t column, const size_t layerId)const{
    const size_t routingBitslLayer = (layerId%2) + firstRoutingBitsLayer_;
    return mapNetworkElement_spaceIndex(rowId,column,routingBitslLayer);
}

vector<FieldElement> witnessMappings::getImageSpaceOrderedBasis()const{
    vector<FieldElement> imageSpaceBasis;
    for(size_t i=0; i< imageSpaceDim_; i++){
        //the i'th element of the standard basis
        FieldElement basisElement = mapIntegerToFieldElement(i,1,1);

        //add to basis
        imageSpaceBasis.push_back(basisElement);
    }
    return imageSpaceBasis;
}

size_t witnessMappings::calculateImageSpaceDim(const commonDeffinitions& commonDef){
    const size_t numOfVariablesLayers = ceil(double(commonDef.variablesNonPerm().size())/double(commonDef.imageWidth()));
    const size_t numOfImageLayers = commonDef.amountOfPermutationLayers() + numOfVariablesLayers;
    const size_t dimOfImageDepth = ceil(Log2(numOfImageLayers));
    const size_t dimOfImage = dimOfImageDepth + commonDef.widthSpaceDimension() + commonDef.heightSpaceDimension();
    return dimOfImage;
}
    
} //namespace BREXtoACSP
} //namespace PCP_Project
