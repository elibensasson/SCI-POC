/******************************************* DecisionAlgorithms.hpp *******************************************/
/**
 * @file.
 *
 * The file DecisionAlgorithms.hpp contains the implementation of PCP decision algorithms, starting from the RS= 
 * verifier all the way up to ACSP. The decision algorithm is a part of the PCP verifier, which is in charge of 
 * getting the results to the queries made by the query algorithm and verify their validity and correctness.
 * 
 * Main classes defined in the file:
 * DecisionAlgorithm - Implements the main decision algorithm functionality, specifically RS=, RS<, RS>, RS, VRS and ACSP.
 * 
 * Main date types defined in the file:
 * (Empty)
 *
 * Main functions defined in the file:
 * Query_Algorithm_XXX - Uses the given randomness in order to collect queries to the proximity proof and evaluation, 
 * and stores the queries using the query builder.
 * 
 * @Note - Each decision algorithm is given a result container as an argument. The result container holds the 
 * results for the queries made by the query algorithm. Notice that most decision algorithms receive an object
 * of type RSXXX_SingleResultsContainer which holds the results for a single verifier, while the amplified decision
 * containers receive XXX_Results_Container which holds the results for all verifiers.
 * 
 */
  /************************************************************************************************************/
#ifndef DECISION_ALGORITHMS_HPP_
#define DECISION_ALGORITHMS_HPP_

#include "languages/ACSP/ACSPInstance.hpp"
#include "PCP/VerifierOnly/RSPCPP_queriesGenerator.hpp"

namespace PCP_Project {

/***********************************************************/
/********** Decision Algorithm Class Definition ************/
/***********************************************************/
/*
 * DecisionAlgorithm class - contains functionality of the PCP decision algorithms.
 */
namespace DecisionAlgorithm {

    typedef struct results_t{
        std::vector<Verifier::RS_PCPP_result> RS_boundary;
        std::vector<Verifier::RS_PCPP_result> RS_composition;
        std::vector<Verifier::ACSP_CONSISTENCY_result> ACSP_consistency;
    };

	/** 
	 * The ACSP decision algorithm.
	 * @param Instance is the instance to the problem, which must be consistent with the polynomial.
	 * @param RandVector is the randomness vector used to test random elements by all the verifiers.
	 * @param results is the result container that holds the queries results.
	 * @param ACSPParams is the set of ACSP parameters.
	 */
	bool Decision_Algorithm_ACSP(const results_t& results);

} //namespace DecisionAlgorithm




}	//Of namespace

#endif /* DECISION_ALGORITHMS_HPP_ */
