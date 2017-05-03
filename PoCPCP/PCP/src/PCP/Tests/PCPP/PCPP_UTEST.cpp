//#include "PCPWorkflowTests.h"
//
//#include "common/Utils/TaskReporting.hpp"
//#include "common/common_Utils.hpp"
//#include "PCP/PCP_common.hpp"
//#include "PCP/SplitVerifiers/QueryAlgorithms.hpp"
//#include "PCP/SplitVerifiers/DecisionAlgorithms.hpp"
//#include "PCP/Verifier/QueriesResFetcher.hpp"
//#include <gtest/gtest.h>
//#include <memory>
//namespace PCP_Project {
//
//using std::string;
//using std::max;
//
//size_t size_of_VRS_results(VRS_Results_Container& results, const unsigned long m){
//	
//    size_t numElements = 0;
//
//    VRS_SingleResultsContainer singleVerifierResults(&results);
//	
//    for (long i=0; i<m; i++) {
//		
//		singleVerifierResults.buildResultsOfSingleVerifier(&results, i);
//        numElements += (singleVerifierResults.RS_singleContainer->singleResultsUnion.singleEqualsResults)->size();
//	}
//    return sizeof(FieldElement)*numElements;
//}
//
//size_t size_of_results(const ACSP_Results_Container& results, const unsigned long m){
//
//    size_t queried_bytes = size_of_VRS_results(*results.VRS_results_Inst,m) + size_of_VRS_results(*results.VRS_results,m);
//    const size_t elementsForConsystency = 1 + results.consistencyResults_p0p1.p0_queries.size();
//
//    queried_bytes += sizeof(FieldElement)*elementsForConsystency;
//    return queried_bytes;
//}
//
//void printSpecs(const size_t total_proof_generated_bytes, const ACSP_Results_Container& results, const QueriesResFetcher::paths_t& paths, const unsigned long m){
//    TASK("Print PCP specs");
//    std::cout<<"The total proof size the prover generated: "<<(total_proof_generated_bytes>>20)<<" MB"<<std::endl;
//
//    //count how many bytes the verifier requested
//    const size_t queried_bytes = size_of_results(results,m);
//    std::cout<<"Total queried data is "<<(queried_bytes>>10)<<" KB"<<std::endl;
//
//    //count total communication complexity
//    size_t communication_complexity = 0;
//    std::vector<CryptoCommitment::hashDigest_t> seen_roots;
//    for(const auto& path : paths){
//        communication_complexity += sizeof(CryptoCommitment::hashDigest_t);
//        for(const auto& h : path.path){
//            if(std::find(seen_roots.begin(),seen_roots.end(),h) == seen_roots.end()){
//                communication_complexity += sizeof(CryptoCommitment::hashDigest_t);
//                seen_roots.push_back(h);
//            }
//        }
//        if(std::find(seen_roots.begin(),seen_roots.end(),path.root) == seen_roots.end()){
//            //root unseen yet, add the root to the communication complexity
//            communication_complexity += sizeof(CryptoCommitment::hashDigest_t);
//            seen_roots.push_back(path.root);
//        }
//    }
//    
//    std::cout<<"Total communication complexity is "<<(communication_complexity>>10)<<" KB"<<std::endl;
//}
//
////WARNING:If using singleQuery = true need to use SINGLE_QUERY flag in SoundnessParams.cpp
//bool PCP_Prove_and_Verify_ACSP(const std::pair<ACSPInstance, ACSPWitness>& ACSPPair, const string randomness_str, bool singleQuery){
//
//	//a counter to count the total proof size
//	//the prover generates in bytes
//	size_t total_proof_generated_bytes = 0;
//	auto update_total_proof_size = [&total_proof_generated_bytes](const BivariateExpansionProof& proof) -> void {
//#pragma omp critical
//		total_proof_generated_bytes += 2 * sizeof(FieldElement)*proof.size();
//	};
//
//	Prover::proof_t proof;
//	{
//		TASK("Executing PCP prover");
//		proof = Prover::PCP_Lang_Producer(ACSPPair);
//
//		//the factor 4 is in order to take in count both proofs with their merkle trees
//		update_total_proof_size(*(proof.compositionRS_proof));
//		update_total_proof_size(*(proof.boundaryRS_proof));
//	}
//
//
//	size_t fieldDim = PCP_common::basisForPCPP(ACSPPair.first).size();
//	//Generating the randomness vector:
//	//There are two ways to do it, 
//	//1- a new random vector, with the proper length
//	//0- a singeltone, containing some constant randomness string for debugging purposes
//	//we use an #if #else switch to chose
//
//	vector<Algebra::details::Randomness> RandVector;
//	unsigned long m = 1;
//	srand(time(0));
//	if (randomness_str == ""){
//		if (!singleQuery) {
//			int neighborNum = ACSPPair.first.neighborPolys().size();
//			m = SoundnessParameters::computeRepetitions(fieldDim, DELTA / neighborNum);//required proximity parameter, and therefore, repretition number of PCPP verifer, can depend on number of neighbors
//		}
//		for (unsigned long i = 0; i<m; i++){
//			RandVector.push_back(Algebra::details::Randomness(RANDOMNESS_RS(fieldDim)));
//		}
//	}
//	else{
//		RandVector.push_back(Algebra::details::Randomness(randomness_str));
//	}
//
//	//#define PRINT_RANDOMNESS
//#ifdef PRINT_RANDOMNESS
//	//print randomness vector for debug purposes
//	std::cout << "Randomness initializes, it is:" << std::endl;
//	for (unsigned long i = 0; i<m; i++){
//		for (unsigned long j = 0; j<RandVector[i].getSize(); j++){
//			std::cout << RandVector[i].getBitByIndex(j);
//		}
//		std::cout << ",";
//	}
//	std::cout << std::endl;
//#endif
//
//	//Creating the queries object:
//	Query_Builder_ACSP builder;
//	FullProofQueryContainer_ACSP* queries = builder.getQueries();
//
//	{
//		TASK("Running " + std::to_string(static_cast<uint64_t>(m)) + " query algorithms"); // the cast is because VS10 doesn't implement all std required overloads
//		QueryAlgorithm::Query_Algorithm_ACSP(&RandVector, &builder, ACSPPair.first);
//	}
//
//	ACSP_Results_Container results(queries);
//	//Ariel:some documentation could help. What exactly is done in the line below? 
//	DecisionAlgorithm::instanceOnlyData verifierPrecomputed(ACSPPair.first, &RandVector);
//
//	QueriesResFetcher::paths_t paths;
//	{
//
//		{
//			TASK("Filling queries results");
//			paths = QueriesResFetcher::fillResults(proof, ACSPPair.first, *queries, results, update_total_proof_size);
//			//delete queries;
//			//REMOVE COMMENT ABOVE!
//
//			printSpecs(total_proof_generated_bytes, results, paths, m);
//		}
//
//	}
//
//	bool res = true;
//	{
//		TASK("Verifying merkle paths");
//		for (const auto& p : paths){
//			const bool currRes = CryptoCommitment::verifyPathToBlock(&p.data[0], p.root, p.path, p.blockIndex);
//			if (!currRes){
//				std::cout << "Merkle path verification failed!" << std::endl;
//				res = false;
//			}
//		}
//	}
//
//	 {
//		 TASK("Calling the ACSP decision algorithms");
//		 res &= DecisionAlgorithm::Decision_Algorithm_ACSP(results, ACSPPair.first, verifierPrecomputed);
//	 }
//
//	 return res;
//}
//
//}	//Of naemspace
