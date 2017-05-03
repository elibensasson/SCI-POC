#ifndef _COMMON_BREXTOACSP_INSTANCEMAPPINGS_HPP__
#define _COMMON_BREXTOACSP_INSTANCEMAPPINGS_HPP__

#include "commonMappings.hpp"

namespace PCP_Project{
namespace BREXtoACSP{
    
class instanceMappings : public commonMappings{
public:
    instanceMappings(const common& commonInfo);
    Algebra::FieldElement mapVariable(const size_t varId)const;
    Algebra::FieldElement mapNonPermutation_zeroRow(const size_t elementId)const;
    Algebra::FieldElement mapNonPermutation_lastRow(const size_t elementId)const;
protected:
    const common& commonInfo_;

private:
    Algebra::FieldElement getLastRowIndicator()const;
};

} //namespace BREXtoACSP
} //namespace PCP_Project


#endif // #ifndef _COMMON_BREXTOACSP_INSTANCEMAPPINGS_HPP__
