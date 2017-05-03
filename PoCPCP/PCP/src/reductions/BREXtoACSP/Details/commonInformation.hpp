/**
 * @file BEXtoACSP_commonInforamtion.hpp
 *
 * Declares the BREXtoACSP reduction common information,
 * as descried at ACSP Spec document.
 */

#ifndef _COMMON_BREXTOACSP_COMMONINFORMATION_HPP_
#define _COMMON_BREXTOACSP_COMMONINFORMATION_HPP_

#include "commonDeffinitions.hpp"
#include "common/Algebra/details/Polynomials.hpp"

namespace PCP_Project{
namespace BREXtoACSP{

/**
 * @class commonInformation
 * @brief represents common values that are restricted, but may be chosen with some freedom
 *
 * We use the BREXPartialInstance : \f$(\mathbb{F},d,\mathcal{V},\mathcal{C}_\mathcal{A}, \mathcal{C}_\pi, \mathscr{I})\f$
 */
class commonInformation{
public:
    commonInformation(const commonDeffinitions& commonDef):
        rowsPrimitivePoly_(Algebra::details::findPrimitive(commonDef.heightSpaceDimension())),
        columnsPrimitivePoly_(Algebra::details::findPrimitive(commonDef.widthSpaceDimension()))
        { };

    /// The modulus for row ids field simulation
    const NTL::GF2X& rowsModulus()const {return rowsPrimitivePoly_;}

    /// The modulus for column ids field simulation
    const NTL::GF2X& columnsModulus()const {return columnsPrimitivePoly_;}

private:
    const NTL::GF2X rowsPrimitivePoly_;
    const NTL::GF2X columnsPrimitivePoly_;
    
};
    
} //namespace BREXtoACSP
} //namespace PCP_Project

#endif //_COMMON_BREXTOACSP_COMMONINFORMATION_HPP_
