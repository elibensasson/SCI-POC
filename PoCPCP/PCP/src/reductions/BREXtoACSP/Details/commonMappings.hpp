#ifndef _COMMON_BREXTOACSP_COMMONMAPPINGS_HPP__
#define _COMMON_BREXTOACSP_COMMONMAPPINGS_HPP__

#include "common.hpp"

namespace PCP_Project{
namespace BREXtoACSP{
    
class commonMappings{

public:
     typedef size_t spaceIndex_t;

     commonMappings(const common& commonInfo);
     Algebra::FieldElement mapNonPermutationElement(const size_t elementId)const;
     
     spaceIndex_t mapPermutationColumnId_spaceIndex(const size_t columnId)const;
     Algebra::FieldElement mapPermutationColumnId_spaceElement(const size_t columnId)const;
     
     spaceIndex_t mapPermutationLayerId_spaceIndex(const size_t layerId)const;
     Algebra::FieldElement mapPermutationLayerId_spaceElement(const size_t layerId)const;
     
     spaceIndex_t mapPermutationElement_spaceIndex(const size_t columnId, const size_t layerId)const;
     Algebra::FieldElement mapPermutationElement_spaceElement(const size_t columnId, const size_t layerId)const;
     ///
     /// Let \f$ \eta \f$ be the rows modulus, of degree \f$ d \f$
     /// This function maps the elementId \f$ \sum_{i=0}^{n} b_i 2^i \f$
     /// to the field element \f$ b_0 \cdot \eta + x^d \cdot \sum_{i=1}^{n} b_i x^i \f$
     ///
     /// For more details ask Michael Riabzev (this method reduces maximal number of neighbors)
     ///
     spaceIndex_t mapNonPermutationVariable_spaceIndex(const size_t elementId)const;
     Algebra::FieldElement mapNonPermutationVariable_spaceElement(const size_t elementId)const;
    
     Algebra::FieldElement map_spaceIndex_to_fieldElement(const spaceIndex_t& i)const;
    
protected:
    const size_t bitsForRowId_;
    const size_t bitsForColumnId_;
    const size_t bitsForLayerId_;
    const size_t bitsForNonPermutationElementId_;
    const size_t firstNonPermUsableIndex_;
    const NTL::GF2X rowsModulus_;
    const spaceIndex_t rowsModulus_bin_;
    const std::vector<spaceIndex_t> column_spaceIndex_;

    static Algebra::FieldElement map_x_power_modulu_poly(const size_t x_power, const NTL::GF2X modulus);
    static size_t GF2X_to_int(const NTL::GF2X& src);
    static std::vector<spaceIndex_t> column_spaceIndex_init(const NTL::GF2X& columnsModulus, const size_t shiftSize);
};

} //namespace BREXtoACSP
} //namespace PCP_Project


#endif // _COMMON_BREXTOACSP_COMMONMAPPINGS_HPP__
