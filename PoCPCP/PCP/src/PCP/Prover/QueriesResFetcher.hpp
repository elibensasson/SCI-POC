#ifndef QUERIESRESFETCHER_HPP__
#define QUERIESRESFETCHER_HPP__

#include "PCP/Prover/Proof/BivariateExpansionProof.hpp"
#include "PCP/Prover/Prover.hpp"
#include "languages/ACSP/ACSPInstance.hpp"
#include "PCP/VerifierOnly/RSPCPP_queriesGenerator.hpp"

namespace PCP_Project{
namespace QueriesResFetcher{

typedef struct queryCommitmentPath{
    size_t blockIndex;
    CryptoCommitment::hashDigest_t root;
    std::vector<Algebra::FieldElement> data;
    CryptoCommitment::path_t path;
};

typedef std::pair<CryptoCommitment::hashDigest_t, size_t> pathIndicator_t;
typedef std::map<pathIndicator_t,queryCommitmentPath> paths_t;

typedef struct queries_t{
    Verifier::RS_PCPP_queriesTree RS_boundary; 
    Verifier::RS_PCPP_queriesTree RS_composition; 
    std::vector<Verifier::ACSP_CONSISTENCY_query> ACSP_consistency;
};


class resultsFillingAlgorithm_interface{
    public:
    virtual paths_t fillResults(const Prover::proof_t& proof, const Prover::ProverAlgorithm& prover,
    const ACSPInstance& instance, const queries_t& queries,
    std::function<void(const BivariateExpansionProof& proof)> collect_proof_spec)const = 0;
};

class resultsFillingAlgorithm : public resultsFillingAlgorithm_interface
{
public:
paths_t fillResults(const Prover::proof_t& proof, const Prover::ProverAlgorithm& prover,
    const ACSPInstance& instance, const queries_t& queries,
    std::function<void(const BivariateExpansionProof& proof)> collect_proof_spec)const;
};

class resultsFillingAlgorithm_dummy : public resultsFillingAlgorithm_interface
{
public:
paths_t fillResults(const Prover::proof_t& proof, const Prover::ProverAlgorithm& prover,
    const ACSPInstance& instance, const queries_t& queries,
    std::function<void(const BivariateExpansionProof& proof)> collect_proof_spec)const;
};

} //namespace QueriesResFetcher
} //namespace PCP_Project

#endif //#ifndef QUERIESRESFETCHER_HPP__
