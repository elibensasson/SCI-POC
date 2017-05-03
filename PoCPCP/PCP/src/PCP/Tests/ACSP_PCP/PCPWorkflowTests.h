
#ifndef PCP_WORKFLOW_TESTS_H_
#define PCP_WORKFLOW_TESTS_H_

#include "languages/ACSP/ACSPInstance.hpp"
#include "languages/ACSP/ACSPWitness.hpp"
#include "PCP/Prover/Prover.hpp"
#include "PCP/Prover/QueriesResFetcher.hpp"
#include <string>

namespace PCP_Project {

bool PCP_Prove_and_Verify_ACSP(const std::pair<ACSPInstance,ACSPWitness>& ACSPPair, const Prover::ProverAlgorithm& prover, const QueriesResFetcher::resultsFillingAlgorithm_interface& resultsFiller);

}	//Of naemspace

#endif /* PCP_WORKFLOW_TESTS_H_ */

