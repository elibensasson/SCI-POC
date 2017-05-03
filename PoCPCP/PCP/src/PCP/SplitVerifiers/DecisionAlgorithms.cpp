/***************************************** DecisionAlgorithms.cpp *********************************************/
/**
 * @file.
 *
 * The file DecisionAlgorithms.cpp contains the implementation of PCP decision algorithms, starting from the RS= 
 * verifier all the way up to ACSP. The decision algorithm is a part of the PCP verifier, which is in charge of 
 * getting the results to the queries made by the query algorithm and verify their validity and correctness.
 * 
 * For more information - Read the documentation of the header file DecisionAlgorithms.hpp.
 */
  /************************************************************************************************************/
#include "PCP/VerifierOnly/RSPCPP_queriesGenerator.hpp"
#include "DecisionAlgorithms.hpp"
#include "common/Infrastructure/Infrastructure.hpp"
#include "../PCP_common.hpp"
#include "common/Algebra/AlgebraCommon.hpp"
#include "common/Utils/TaskReporting.hpp"
#include <algebraLib/UnivariatePolynomialGeneral.hpp>

using namespace std;
using Infrastructure::POW2;
using Infrastructure::Log2;
using Algebra::zero;
using Algebra::UnivariatePolynomialGeneral;
using Algebra::PolynomialDegree;

namespace PCP_Project {
namespace DecisionAlgorithm{

bool Decision_Algorithm_ACSP(const results_t& results) {
    
    
    bool res = true;//proof valid until proven non-valid
    
    {
        TASK("verifying low degree of boundary proof");
        bool localRes = true;
        for(const auto& res : results.RS_boundary){
            localRes &= res.verify();
        }
        
        if(localRes == true){
            std::cout<<"OK"<<std::endl;
        }
        else {
            std::cout<<"Failed!!!"<<std::endl;
            res =  false; 
        }
    }
    
    {
        TASK("verifying low degree of composition proof");
        bool localRes = true;
        for(const auto& res : results.RS_composition){
            localRes &= res.verify();
        }
        
        if(localRes == true){
            std::cout<<"OK"<<std::endl;
        }
        else {
            std::cout<<"Failed!!!"<<std::endl;
            res =  false; 
        }
    }
    
    {
        const auto& consistencyRes = results.ACSP_consistency;
        TASK("verifying consistency of composition with ACSP witness on " + std::to_string(consistencyRes.size()) + " queries");

        bool localRes = true;
        for(const auto& res : consistencyRes){
            try{
                localRes &= res.verify();
            }
            catch(...){
                localRes = false;
            }
        }
        
        if(localRes == true){
            std::cout<<"OK"<<std::endl;
        }
        else {
            std::cout<<"Failed!!!"<<std::endl;
            res =  false; 
        }
    }

	return res;
}


} //namespace DecisionAlgorithm
}	//Of namespace PCP_Project
