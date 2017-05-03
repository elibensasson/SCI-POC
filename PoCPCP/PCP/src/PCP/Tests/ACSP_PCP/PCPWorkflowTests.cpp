#include "PCPWorkflowTests.h"

#include "common/Utils/TaskReporting.hpp"
#include "PCP/PCP_common.hpp"
#include "PCP/SplitVerifiers/DecisionAlgorithms.hpp"
#include "PCP/Prover/QueriesResFetcher.hpp"

#include <memory>
#include <unordered_set>
#ifndef _MSC_VER
#define NOEXCEPT noexcept
#else
#define NOEXCEPT
#endif
namespace PCP_Project {

using std::string;
using std::max;

size_t size_of_results(const DecisionAlgorithm::results_t& results){

    size_t queried_elements= 0;
    for(const auto& r : results.RS_boundary){
        queried_elements += r.results.size();
    }
    for(const auto& r : results.RS_composition){
        queried_elements += r.results.size();
    }
    for(const auto& r : results.ACSP_consistency){
        queried_elements += 1 + r.boundaryPoly_res.size();
    }

    return queried_elements*sizeof(FieldElement);
}

void printSpecs(const size_t total_proof_generated_bytes, const DecisionAlgorithm::results_t& results, const QueriesResFetcher::paths_t& paths){
    std::cout<<"The total proof size the prover generated: "<<(total_proof_generated_bytes>>20)<<" MB"<<std::endl;

    //count how many bytes the verifier requested
    const size_t queried_bytes = size_of_results(results);
    std::cout<<"Total queried data is "<<(queried_bytes>>10)<<" KB"<<std::endl;

    return;

    //count total communication complexity
    size_t communication_complexity = 0;

    struct hashFunc{
        size_t operator()(const CryptoCommitment::hashDigest_t& x)const NOEXCEPT{
            return *((size_t *) &x);
        }
    };
    std::unordered_set<CryptoCommitment::hashDigest_t,hashFunc> seen_roots;
    
    for(const auto& path_pair : paths){
        const auto& path = path_pair.second;
        communication_complexity += sizeof(CryptoCommitment::hashDigest_t);
        for(const auto& h : path.path){
            if(seen_roots.count(h) == 0){
                communication_complexity += sizeof(CryptoCommitment::hashDigest_t);
                seen_roots.insert(h);
            }
        }
        if(seen_roots.count(path.root) == 0){
            //root unseen yet, add the root to the communication complexity
            communication_complexity += sizeof(CryptoCommitment::hashDigest_t);
            seen_roots.insert(path.root);
        }
    }
    
    std::cout<<"Total communication complexity is "<<(communication_complexity>>10)<<" KB"<<std::endl;
}

bool PCP_Prove_and_Verify_ACSP(const std::pair<ACSPInstance, ACSPWitness>& ACSPPair, const Prover::ProverAlgorithm& prover, const QueriesResFetcher::resultsFillingAlgorithm_interface& resultsFiller){

    //
    // Check if the field is big enough
    //
    {
        const short PCPP_spaceDim = PCP_common::basisForPCPP(ACSPPair.first).size();
        const short ContextField_Dim = ACSPPair.first.contextField().degree();
        _COMMON_ASSERT((PCPP_spaceDim <= ContextField_Dim),
            string("Can't construct a PCP proof for ACSP instance, because the context field is too small\n") + 
            string("The context field is of dimension ") + std::to_string(ContextField_Dim) + string("\n") +
            string("The sub-space for PCPP proof must be at least of dimension ") + std::to_string(PCPP_spaceDim)
        );
    }

	//a counter to count the total proof size
	//the prover generates in bytes
	size_t total_proof_generated_bytes = 0;
	auto update_total_proof_size = [&total_proof_generated_bytes](const BivariateExpansionProof& proof) -> void {
#pragma omp critical
		total_proof_generated_bytes += 2 * sizeof(FieldElement)*proof.size();
	};

    //
    //Constructing proof for first layer
    //
	Prover::proof_t proof;
	{
		TASK("Executing PCP prover");
		proof = prover.PCP_Lang_Producer(ACSPPair);

		//the factor 4 is in order to take in count both proofs with their merkle trees
		update_total_proof_size(*(proof.compositionRS_proof));
		update_total_proof_size(*(proof.boundaryRS_proof));
	}

    //
    // NEW QUERY ALGORITHM
    //
    QueriesResFetcher::queries_t queries;
    DecisionAlgorithm::results_t results;
    {
        TASK("VERIFIER QUERIES GENERATION");
        
        const auto PCPP_Basis(PCP_common::basisForPCPP(ACSPPair.first));
        const size_t neighborNum = ACSPPair.first.neighborPolys().size();
        const size_t numPCPP_queries = SoundnessParameters::computeRepetitions(PCPP_Basis.size(), DELTA / neighborNum);
        {
            TASK("Generating " + std::to_string(static_cast<uint64_t>(ACSP_CONSISTENCY_SOUNDNESS)) + " ACSP consistency queries");
            
            results.ACSP_consistency.resize(ACSP_CONSISTENCY_SOUNDNESS);
            const auto basisConsistency = PCP_common::basisForConsistency(ACSPPair.first);
            for(size_t i=0; i< (ACSP_CONSISTENCY_SOUNDNESS) ; i++){
                const FieldElement consistencyPoint = Algebra::getSpaceElementByIndex(basisConsistency, Algebra::zero(), rand());

                results.ACSP_consistency[i].init(ACSPPair.first,consistencyPoint);
                queries.ACSP_consistency.push_back(Verifier::ACSP_CONSISTENCY_query(consistencyPoint,results.ACSP_consistency[i]));
            }
        }
        {
            TASK("Generating " + std::to_string(static_cast<uint64_t>(numPCPP_queries)) + " RS_PCPP queries for boundary polynomial:");
            results.RS_boundary.resize(numPCPP_queries);
            const auto deg_boundary = Algebra::PolynomialDegree::integral_t(PCP_common::witness_div_Z_Boundery_degreeBound(ACSPPair.first));
            const short deg_log_boundary = ceil(Infrastructure::Log2(deg_boundary));
            std::cout<<"Polynomial evaluation space dimension is "<<PCPP_Basis.size()<<std::endl;
            std::cout<<"Polynomial log of degree bound is "<<deg_log_boundary<<std::endl;
            std::cout<<"PCPP recursion depth is "<<RECURSION_DEPTH<<std::endl;
            for(size_t i=0; i< numPCPP_queries; i++){
                Verifier::addRandomQuery(deg_log_boundary, PCPP_Basis, Algebra::zero(), queries.RS_boundary, results.RS_boundary[i],RECURSION_DEPTH,true);
            }
        }
        {
            TASK("Generating " + std::to_string(static_cast<uint64_t>(numPCPP_queries)) + " RS_PCPP queries for composition polynomial:");
            results.RS_composition.resize(numPCPP_queries);
            const auto deg_composition = Algebra::PolynomialDegree::integral_t(PCP_common::composition_div_ZH_degreeBound(ACSPPair.first));
            const short deg_log_composition = ceil(Infrastructure::Log2(deg_composition));
            std::cout<<"Polynomial evaluation space dimension is "<<PCPP_Basis.size()<<std::endl;
            std::cout<<"Polynomial log of degree bound is "<<deg_log_composition<<std::endl;
            std::cout<<"PCPP recursion depth is "<<RECURSION_DEPTH<<std::endl;
            for(size_t i=0; i< numPCPP_queries; i++){
                Verifier::addRandomQuery(deg_log_composition, PCPP_Basis, Algebra::zero(), queries.RS_composition, results.RS_composition[i],RECURSION_DEPTH,true);
            }
        }
    }
	
	QueriesResFetcher::paths_t paths;
	{

		{
			TASK("Filling queries results");
			paths = resultsFiller.fillResults(proof, prover, ACSPPair.first, queries,update_total_proof_size);
			printSpecs(total_proof_generated_bytes, results, paths);
		}

	}

	bool res = true;
	{
		TASK("Verifying merkle paths");
		for (const auto& p_pair : paths){
            const auto& p = p_pair.second;
			const bool currRes = CryptoCommitment::verifyPathToBlock(&p.data[0], p.root, p.path, p.blockIndex);
			if (!currRes){
				res = false;
			}
		}
	}

	 {		 TASK("Calling the ACSP decision algorithms");	
		 res &= DecisionAlgorithm::Decision_Algorithm_ACSP(results);
	 }

	 return res;
}

}	//Of naemspace
