#include "commonMappings.hpp"

using Algebra::FieldElement;
using Algebra::mapIntegerToFieldElement;
using Infrastructure::Log2;
using std::max;

namespace PCP_Project{
namespace BREXtoACSP{

commonMappings::commonMappings(const common& commonInfo):
    bitsForRowId_(commonInfo.heightSpaceDimension()),
    bitsForColumnId_(commonInfo.widthSpaceDimension()),
    bitsForLayerId_(commonInfo.contextField().degree() - (commonInfo.heightSpaceDimension()+commonInfo.widthSpaceDimension())),
    bitsForNonPermutationElementId_(commonInfo.contextField().degree() - commonInfo.heightSpaceDimension()),
    firstNonPermUsableIndex_(commonInfo.imageWidth()*commonInfo.amountOfPermutationLayers()),
    column_spaceIndex_(column_spaceIndex_init(commonInfo.columnsModulus(),commonInfo.heightSpaceDimension())),
    rowsModulus_(commonInfo.rowsModulus()),
    rowsModulus_bin_(GF2X_to_int(commonInfo.rowsModulus())){};

FieldElement commonMappings::mapNonPermutationElement(const size_t elementId)const{
    return mapIntegerToFieldElement(bitsForRowId_,bitsForNonPermutationElementId_, firstNonPermUsableIndex_ + elementId);
}

commonMappings::spaceIndex_t commonMappings::mapNonPermutationVariable_spaceIndex(const size_t elementId)const{
//#define USE_MICHAELS_NEIGHBORS_TRICK
#ifdef USE_MICHAELS_NEIGHBORS_TRICK
    const size_t elementId_clearedLSB = (elementId>>1)<<1;
    const spaceIndex_t offset = (elementId%2 == 0? 0 : rowsModulus_bin_);
    return ((firstNonPermUsableIndex_ + elementId_clearedLSB)<<bitsForRowId_) ^ offset;

#else
    return ((firstNonPermUsableIndex_ + elementId)<<bitsForRowId_);
#endif
}

FieldElement commonMappings::mapNonPermutationVariable_spaceElement(const size_t elementId)const{
    const auto spaceIndex = mapNonPermutationVariable_spaceIndex(elementId);
    return map_spaceIndex_to_fieldElement(spaceIndex);
}

commonMappings::spaceIndex_t commonMappings::mapPermutationColumnId_spaceIndex(const size_t columnId)const{
    return column_spaceIndex_[columnId];
}

FieldElement commonMappings::mapPermutationColumnId_spaceElement(const size_t columnId)const{
    const auto spaceIndex= mapPermutationColumnId_spaceIndex(columnId);
    return map_spaceIndex_to_fieldElement(spaceIndex);
}

commonMappings::spaceIndex_t commonMappings::mapPermutationLayerId_spaceIndex(const size_t layerId)const{
    const size_t bitsOffset = bitsForRowId_ + bitsForColumnId_;
    return layerId<<bitsOffset;
}

FieldElement commonMappings::mapPermutationLayerId_spaceElement(const size_t layerId)const{
    const auto spaceIndex= mapPermutationLayerId_spaceIndex(layerId);
    return map_spaceIndex_to_fieldElement(spaceIndex);
}

commonMappings::spaceIndex_t commonMappings::mapPermutationElement_spaceIndex(const size_t columnId, const size_t layerId)const{
    return mapPermutationColumnId_spaceIndex(columnId) ^ mapPermutationLayerId_spaceIndex(layerId);
}

FieldElement commonMappings::mapPermutationElement_spaceElement(const size_t columnId, const size_t layerId)const{
    const auto spaceIndex= mapPermutationElement_spaceIndex(columnId,layerId);
    return map_spaceIndex_to_fieldElement(spaceIndex);
}
    
FieldElement commonMappings::map_x_power_modulu_poly(const size_t x_power, const NTL::GF2X modulus){
    
    using NTL::GF2X;
    using NTL::SetCoeff;
    using NTL::to_GF2E;

    GF2X x_power_poly;
    SetCoeff(x_power_poly, x_power);

    return FieldElement(to_GF2E(x_power_poly % modulus));

}

FieldElement commonMappings::map_spaceIndex_to_fieldElement(const spaceIndex_t& i)const{
    const size_t numBits = bitsForRowId_ + bitsForColumnId_ + bitsForLayerId_;
    return mapIntegerToFieldElement(0 , numBits, i);
}

size_t commonMappings::GF2X_to_int(const NTL::GF2X& src){
   const auto asGF2E = to_GF2E(src);
   const auto asInt = mapFieldElementToInteger(0,deg(src)+1, FieldElement(asGF2E));
   return asInt;
}

vector<commonMappings::spaceIndex_t> commonMappings::column_spaceIndex_init(const NTL::GF2X& columnsModulus, const size_t shiftSize){
    const spaceIndex_t modulus_bin(GF2X_to_int(columnsModulus));
    const spaceIndex_t overflow_mask = Infrastructure::POW2(deg(columnsModulus));
    
    vector<spaceIndex_t> res;
    res.push_back(1L<<shiftSize);
    
    //if the modulus is of degree at most 0,
    //there should be only one element (the field is GF(2))
    if(deg(columnsModulus) < 1){
        return res;
    }
    
    while(true){
        //next element is multiplication by x (done by shift left)
        //modulu the modulus (done by xor if overflow)
        spaceIndex_t nextElem = (res[res.size()-1]>>shiftSize)<<1;
        if (nextElem & overflow_mask){
            nextElem ^= modulus_bin;
        }

        //if we got back to 1, we finished going over all the elements
        if(nextElem == 1){
            return res;
        }

        //add curr result
        res.push_back(nextElem<<shiftSize);
    }
}
    
} //namespace BREXtoACSP
} // namespace PCP_Project
