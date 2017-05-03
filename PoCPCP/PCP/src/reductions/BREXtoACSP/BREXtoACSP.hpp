/************************************** TinyRAMtoBREX.hpp ****************************************/

/**
 * @file BEXtoACSP.hpp
 *
 * Declares the BREXtoACSP reduction class
 */

  /***********************************************************************************************/

#pragma once
#ifndef __BREX_TO_ACSP_HPP
#define __BREX_TO_ACSP_HPP

#include "languages/BREX/BrexInstance.hpp"
#include "languages/BREX/BrexWitness.hpp"
#include "languages/ACSP/ACSPInstance.hpp"
#include "languages/ACSP/ACSPWitness.hpp"

#include <memory>

namespace PCP_Project{


class CBREXtoACSP{
public:
    /**
    * Reduces an BREXInstance to an ACSPInstance.
    * @param[in]  partialInstance
    * @return a pointer to an ACSPInstance.
    */
    static ::std::unique_ptr<ACSPInstance> reduceInstance(const BREXInstance&  instance);

    /**
    * Reduces a BREXExecutedWitness to an ACSPWitness.
    * @param[in]  fullInstance
    * @param[in]  witness
    * @return a pointer to an ACSPWitness.
    */
    static ::std::unique_ptr<ACSPWitness> reduceWitness(const BREXInstance& instance, const BREXWitness& witness);

}; // class BREXtoACSP

}//namespace PCP_Project

#endif // __BREX_TO_ACSP_HPP
