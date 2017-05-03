/**
 * @file BEXtoACSP_common.hpp
 *
 * Declares the BREXtoACSP reduction common,
 * as descried at ACSP Spec document.
 */

#ifndef _COMMON_BREXTOACSP_COMMON_HPP_
#define _COMMON_BREXTOACSP_COMMON_HPP_

#include "commonDeffinitions.hpp"
#include "commonInformation.hpp"
//#define ONE_CONSTRAINT //ARIEL-For debugging purposes to get to pcp part faster 
namespace PCP_Project{
namespace BREXtoACSP{

/*
 * @class common
 * @brief common values, both well defined and those that have some freedom, united in one class
 */
class common : public commonDeffinitions, public commonInformation{
    public:
        common(const BREXInstance& instance):
            commonDeffinitions(instance),
            commonInformation(*((commonDeffinitions*)this))
            {};
        
        common(const common&) = delete;
};
    
} //namespace BREXtoACSP
} //namespace PCP_Project

#endif //_COMMON_BREXTOACSP_COMMONINFORMATION_HPP_
