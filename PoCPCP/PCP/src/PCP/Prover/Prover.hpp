/********************************************** Prover.hpp ************************************************/
/**
 * @file.
 *
 * The file Prover.hpp contains the implementation of PCP prover, starting from the RS= prover all the way up to proverLang.
 * 
 * Main classes defined in the file:
 * Prover - Implements the main prover functionality, specifically RS=, RS<, RS>, RS, VRS, ACSP and proverLang.
 * 
 * Main date types defined in the file:
 * (Empty)
 *
 * Main functions defined in the file:
 * Each prover (not including ACSP and proverLang) is comprised of two functions - XXX_Producer and Produce_XXX_proximityProof.
 * Notice that the proof is actually a pair (p,pi), where p is the evaluation of a polynomial and pi is the proximity proof.
 * XXX_Producer is in charge of computing the evaluation p and calling the Produce_XXX_proximityProof function.
 * Produce_XXX_proximityProof is in charge of computing the proximity proof pi. 
 * 
 * @Note that each prover is given a proof consumer. We use the producer-consumer design pattern in the following manner:
 * The prover produces a proof consumed by the verifier. There are several different consumers available (and the prover is 
 * oblivious as to which consumer he is given). For example - The InMemory consumer receives the proof (produced by the prover)
 * and stores it to memory. The Query consumer stores only the proof relevant for the specific queries etc.
 * @Note : There is a lot of `bueracracy' involved consisting of many methods , but where most of the actions happens is 
 * Produce_RSEquals_ProximityProof_Eval where the actual PCP recursion takes place
 * 
 */
  /************************************************************************************************************/
#ifndef PROVER_HPP_
#define PROVER_HPP_

//#define SABOTAGE_PROOF	//If defined, the RS= prover will sabotage a few points in the proof. (Used for soundness tests).

#include "languages/ACSP/ACSPInstance.hpp"
#include "languages/ACSP/ACSPWitness.hpp"
#include "common/Algebra/LightUniPolyEval.hpp"
#include "Proof/BivariateExpansionProof.hpp"
#include "common/Utils/TaskReporting.hpp"
#include <memory>

namespace PCP_Project {

namespace Prover {
    
    typedef struct proof_t{
        std::unique_ptr<BivariateExpansionProof> compositionRS_proof;
        std::unique_ptr<BivariateExpansionProof> boundaryRS_proof;

		  proof_t(){};
		  proof_t(proof_t&& src) :compositionRS_proof(std::move(src.compositionRS_proof)), boundaryRS_proof(std::move(src.boundaryRS_proof)){};
		  void operator=(proof_t&& src){
			  compositionRS_proof = std::move(src.compositionRS_proof);
			  boundaryRS_proof = std::move(src.boundaryRS_proof);
		  }
    };

class ProverAlgorithm{
    public:	
	/**
	 * The main prover of the language.
	 * @param ACSPPair is the set tuple of ACSP (instance,witness).
	 */
    virtual proof_t PCP_Lang_Producer(const std::pair<ACSPInstance,ACSPWitness>& ACSPPair)const;

    /**
     * Generation of low degree (aka RS) proof for some evaluation of a function
     */
    virtual std::unique_ptr<BivariateExpansionProof> expandUnivariateToBivariate(Algebra::FieldElement const * const evaluation, const std::vector<Algebra::FieldElement>& evaluationBasis, const Algebra::FieldElement& affineShift)const;
    
    /**
     * Generation of low degree (aka RS) proof for a row of a bivariate expansion
     */
    virtual std::unique_ptr<BivariateExpansionProof> proofRow_lowDegree(const BivariateExpansionProof& src, const size_t row_index ,const std::vector<Algebra::FieldElement>& rowBasis, const Algebra::FieldElement& affineShift)const;
    
    /**
     * Generation of low degree (aka RS) proof for a column of a bivariate expansion
     */
    virtual std::unique_ptr<BivariateExpansionProof> proofColumn_lowDegree(const BivariateExpansionProof& src, const size_t column_index ,const std::vector<Algebra::FieldElement>& columnBasis, const Algebra::FieldElement& affineShift)const;
};

}
}

#endif /* PROVER_HPP_ */
