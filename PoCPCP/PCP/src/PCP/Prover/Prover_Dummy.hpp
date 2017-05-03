#ifndef PROVER_DUMMY_HPP__
#define PROVER_DUMMY_HPP__

#include "Prover.hpp"

namespace PCP_Project {
namespace Prover {

class Prover_Dummy : public ProverAlgorithm{
	public:
    /**
	 * The main prover of the language.
	 * @param ACSPPair is the set tuple of ACSP (instance,witness).
	 */
    proof_t PCP_Lang_Producer(const std::pair<ACSPInstance,ACSPWitness>& ACSPPair)const;

    /**
     * Generation of low degree (aka RS) proof for some evaluation of a function
     */
    std::unique_ptr<BivariateExpansionProof> expandUnivariateToBivariate(Algebra::FieldElement const * const evaluation, const std::vector<Algebra::FieldElement>& evaluationBasis, const Algebra::FieldElement& affineShift)const;
    
    /**
     * Generation of low degree (aka RS) proof for a row of a bivariate expansion
     */
    std::unique_ptr<BivariateExpansionProof> proofRow_lowDegree(const BivariateExpansionProof& src, const size_t row_index ,const std::vector<Algebra::FieldElement>& rowBasis, const Algebra::FieldElement& affineShift)const;
    
    /**
     * Generation of low degree (aka RS) proof for a column of a bivariate expansion
     */
    std::unique_ptr<BivariateExpansionProof> proofColumn_lowDegree(const BivariateExpansionProof& src, const size_t column_index ,const std::vector<Algebra::FieldElement>& columnBasis, const Algebra::FieldElement& affineShift)const;
};


} //namespace Prover
} //namespace PCP_Project
#endif //PROVER_DUMMY_HPP__
