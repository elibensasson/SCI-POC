/********************************************** SoundenessParams.hpp ************************************************/
/**
 * @file.
 *
 * The file SoundenessParams.hpp contains the implementation of the parameters relevant to the soundness analysis, i.e. 
 * the parameters that are used to compute the amplified verifier's number of repetitions.
 * 
 * The parameters are computed according to the paper ``Improved concrete and feasible soundness of the
 *  Ben-Sasson-Sudan PCPP for practical input lengths'', referred to in this file and SoundnessParams.cpp as ``the paper''
 * Some commented out code contains calculations of parameters according to the paper ``On the concrete efficiency of PCPs''
 * referred to as ``concrete paper''
 * Main classes defined in the file:
 * SoundnessParameters - Holds the relevant soundness parameters and is in charge of computing number of repetitions.
 * 
 * Main date types defined in the file:
 * ProximityParam and SoundnessParam - data types for the corresponding parameters.
 *
 * Main functions defined in the file:
 * SoundnessParameters Functions - given a set of soundness parameters, computes the number of repetitions.
 * Randomness Macros - Compute the randomness complexity (i.e. number of random bits required) for the PCP verifiers.
 */
  /************************************************************************************************************/
#ifndef VERIFIERPARAMS_HPP_
#define VERIFIERPARAMS_HPP_

#include "common/Algebra/details/FiniteFields.hpp"

namespace PCP_Project {

/** The following are macros of soundness parameters: */
#define Gamma 				SoundnessParameters::gamma
#define Mu 					SoundnessParameters::mu
#define Eta 				SoundnessParameters::eta
#define FIXED_DEPTH SoundnessParameters::useFixedDepth
#define RECURSION_DEPTH  SoundnessParameters::recursionDepth
#define USE_FEASIBLE SoundnessParameters::useFeasible
#define DELTA   SoundnessParameters::deltaVRS
#define ACSP_CONSISTENCY_SOUNDNESS SoundnessParameters::consistencyRepNum
#define BIASED_RECURSION ;//when this flag is defined, verifier chooses a column or row in a biased way according to chooseRowOrCol function
/** Represents the verifier type for the amplified versions. */
typedef enum {witness,composition} AMP_Type;  
typedef enum{ ROW, COLUMN } queryType_t;


/**
 * The class holds and computes the soundness parameters relevant for the RS & VRS verifiers.
 */
typedef double ProximityParam;		///Used for representing delta in the range (0,1]

class SoundnessParameters {

public:

	/**
	 * We use the following parameters in order to manipulate the behavior of the algorithm with the purpose of
	 * finding the best trade-off between query complexity and number of repetitions needed to achieve a certain
	 * constant soundness. The meaning of each parameter is described in practical-PCPs-notes section 6.1.
	 */

	static int gamma;					///L0 contains [k/2] basis elements minus gamma.
	static int mu;						///L0' has mu more dimensions than L0.
	static int eta;						///The degree we test in RS= is |L|/(2^eta)-1, e.g. |L|/8-1 for eta=3.
	static bool useFeasible;            //should we use the feasible(conjectured) soundness bound for computing # of verifer reps
	static ProximityParam deltaVRS;		///The proximity parameter for AMP-VRS
	static int consistencyRepNum; //number of reps for ACSP consistency check to achieve desired consistency error
	static int securityParam; // the log, base 1/2, of desired soundness error
	static bool useFixedDepth;   //a flag for using a fixed PCPP recursion depth rather than 
	static int recursionDepth; //The desired recursion depth (when desiring to reach a certain depth, rather than a certain size base case)
	static void setupSoundnessParams();
	static void setRecursionDepth(const short& basisSize); //setting recursion depth according to dimension of PCPP space.
	//if dim is k, then recursion depth should be largest d s.t. 2d2^d<k. This will ensure every recursive row has intersection with original poly (see Remark 4.11 in paper)

	/**
	* The function computes the number of repetitions done by AMP-VRS verifiers.
	* The formulas are according to the paper
	* @param k is the size of the basis.
	* @param delta is the parameter we want to test proximity to (need to get this as parameter as can depend on number of ACSP neighbors)
	*/
	static unsigned long computeRepetitions(const int k, const ProximityParam& delta);
	static queryType_t chooseRowOrCol(int depth); //choose whether to recurse on row or column to maximize rejection probability (currently againt row\col compliant attackes on random functions)
};



/** 
* Assertions on the parameters, as required according to practical-PCPs-notes section
* 6.1.3 (7).
*/
extern bool areSoundnessParametersValid();


}  //Of namespace PCP_Project

#endif /* VERIFIERPARAMS_HPP_ */
