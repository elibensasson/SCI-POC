/************************************************* FFT.hpp **************************************************/
/**
 * @file.
 *
 * The file FFT.hpp contains the implementation of interpolation and evaluation functions for univariate 
 * polynomials, implemented using FFT algorithms.
 * 
 * Main classes defined in the file:
 * (Empty)
 * 
 * Main date types defined in the file:
 * FFT algorithms can sometimes be less efficient than the naive algorithms, specifically for small instances. 
 * We define the constants INTERPOLATION_FFT_LOWER_BOUND and EVALUATION_FFT_LOWER_BOUND as the minimal dimensions 
 * from which we use the FFT algorithms. Note that the cosntants were measured on a specific architecture, and might 
 * be different in different computers. So it is advised to verify the constants once in a while.
 *
 * Main functions defined in the file:
 * GathenGerhardInterpolation - Implements an interpolation algorithm that uses the special characteristics of GF2 extension spaces. The 
 * algorithms were taken from the Mateer's thesis which can be found under the folder \references\computation in finite fields.
 * InterpolatePoly - The interpolation algorithm implemented by Arnab. The function is for testing purposes only and is generally not in use.
 *
 * GathenGerhardEvaluation - Implements an evaluation algorithm that uses the special characteristics of GF2 extension spaces. The 
 * algorithms were taken from the Mateer's thesis which can be found under the folder \references\computation in finite fields.
 * fft_additive - The evaluation algorithm implemented by Arnab. The function is for testing purposes only and is generally not in use.
 *
 * @Note - The pseudo code for the Gathen-Gerhard algorithms can be found in SVN under \notes\Implementation Notes\Algorithms\GathenGerhardFFT
 */
  /************************************************************************************************************/
#ifndef FFT_HPP_
#define FFT_HPP_

#include "Polynomials.hpp"
#include <NTL/GF2EX.h>

namespace Algebra {
namespace details {

///Forward Declaration:
class VanishingPolynomial;
class UniPolyEval_forProof;

/******************************************************************************/
/*********************** Interpolation Algorithms *****************************/
/******************************************************************************/

	/**
	* If the basis size >= INTERPOLATION_FFT_LOWER_BOUND: Use Mateer's FFT interpolation algorithm.
	* Otherwise use NTL's interpolation algorithm.
	*/
	#define INTERPOLATION_FFT_LOWER_BOUND 6

	/**
	* If the basis size >= INTERPOLATION_MULTITHREADED_FFT_LOWER_BOUND: Use Gathen-Gerhard FFT interpolation algorithm with multi-threading.
	* Otherwise use single thread.
	*/
	#define INTERPOLATION_MULTITHREADED_FFT_LOWER_BOUND 12 

	/** 
	* The function GathenGerhardInterpolation is a wrapper for Mateer's interpolation algorithm, which does the necessary pre-computations,
	* and then calls the recursive algorithm that does the actual FFT.
	* @param result is the univariate polynomial (coefficients) that's the result of the interpolation function. 
	* @param evaluation is the evaluation of the polynomial, i.e. the list of (point value) pairs over the entire subspace.
	* @param Irr is the irreducible polynomial that defines the field extension.
	* @param basis is the basis that spans the space over which we interpolate.
	* @param spanned is the spanned space (shifted in case the shift is non-zero).
	* @param shift is the affine shift in case the subspace is shifted.
	* @param numThreads is the number of threads OpenMP will execute the FFT on. Default is single thread.
    * @param curThreadNum is a unique thread id used for NTL memory management in NTL/GF2X1.cpp. Thread ID management is done using this parameter.
	* @Note - The pseudo code for the Gathen-Gerhard algorithms can be found in SVN under \notes\Implementation Notes\Algorithms\GathenGerhardFFT
	*/

	void GathenGerhardInterpolation(UnivariatePolynomial& result, const UniPolyEval_forProof* evaluation, const NTL::GF2X& Irr, const Basis& basis, const ExplicitLinearSubset& spanned, const FieldElement& shift = Algebra::zero());
	void GathenGerhardInterpolation_Aux(UnivariatePolynomial& result, const UniPolyEval_forProof* evaluation, const Basis& basis, const BasisIndex i, const FieldIndex j, VanishingPolynomial* vanishingPolys, 
											const ExplicitLinearSubset& spanned, const FieldElement& shift, const unsigned int numThreads, const unsigned int curThreadNum);

	/**
	* Arnab's Interpolation FFT:
	* The following functions are used for efficient interpolation of a univariate polynomial over GF(2^k).
	* @Time complexity is O(n*(log(n)^2)).
	*/
	//void interpolate_poly(NTL::GF2EX& ans, const UniPolynomialEvaluation& values, const Basis& basis, const ExplicitLinearSubset spanned);
	void combine_poly(NTL::GF2EX& p, const NTL::GF2EX& p0, const NTL::GF2EX& p1, const NTL::GF2EX& q, long n);
	void expt2power(FieldElement& y, const FieldElement& x, int degree);

/******************************************************************************/
/************************** Evaluation Algorithms *****************************/
/******************************************************************************/
	
	/**
	* If the basis size >= EVALUATION_FFT_LOWER_BOUND: Use the FFT evaluation algorithm.
	* Otherwise use point-by-point evaluation.
	*/
	#define EVALUATION_FFT_LOWER_BOUND 12

	/**
	* If the basis size >= EVALUATION_MULTITHREADED_FFT_LOWER_BOUND: Use Gathen-Gerhard FFT evaluation algorithm with multi-threading
	* Otherwise use single thread.
	*/
	#define EVALUATION_MULTITHREADED_FFT_LOWER_BOUND 8

	/**
	* Mateer's Evaluation FFT:
	* The following functions are used for efficient evaluation of a univariate polynomial over GF(2^k).
	*/

	static const int GATHEN_GERHARD_OPTIMAL_BASE = 1;		///The optimal base case size for the Gathen-Gerhard FFT algorithm

	/** 
	* The function GathenGerhardEvaluation is a wrapper for Mateer's evaluation algorithm, which does the necessary pre-computations,
	* and then calls the recursive algorithm that does the actual FFT.
	* @param result_Array is the array of field elements that are the result of the evaluation (NULL if we are using maps for evaluation)
	* @param result_Map is the mapping of field elements that are the result of the evaluation (NULL if we are using arrays for evaluation)
	* @param Irr is the irreducible polynomial that defies the field extension.
	* @param poly is the polynomial to be evaluated.
	* @param basis is the basis that spans the field we're evaluation on.
	* @param shift is the affine shift in case the subspace is shifted.
	* @param numThreads is the number of threads OpenMP will execute the FFT on. Default is single thread.
    * @param curThreadNum is a unique thread id used for NTL memory management in NTL/GF2X1.cpp. Thread ID management is done using this parameter.
	* @Note - The pseudo code for the Gathen-Gerhard algorithms can be found in SVN under \notes\Implementation Notes\Algorithms\GathenGerhardFFT
	*/
	void GathenGerhardEvaluation(FieldElement* result_Array, UniPolynomialEvaluation* result_Map, const NTL::GF2X& Irr, const UnivariatePolynomial& poly, const Basis& basis, const FieldElement& shift = Algebra::zero());
	void GathenGerhardEvaluation_Aux(FieldElement* result_Array, UniPolynomialEvaluation* result_Map, const UnivariatePolynomial& f, const Basis& basis, const BasisIndex i, const FieldIndex j, 
												VanishingPolynomial* vanishingPolys, ExplicitLinearSubset& spanned, const FieldElement& shift, const unsigned int numThreads, const unsigned int curThreadNum);

	/**
	* Arnab's Evaluation FFT:
	* The following functions are used for efficient evaluation of a univariate polynomial over GF(2^k).
	*/
	void fft_additive(UniPolynomialEvaluation& ans, const UnivariatePolynomial& poly, const Basis& basis);
	void find_couple(NTL::GF2EX& p0, NTL::GF2EX& p1, const NTL::GF2EX& poly, const NTL::GF2EX& q, long n);



} // of namespace details
} // of namespace Algebra

#endif	//FFT_HPP_
