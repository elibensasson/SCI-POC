/************************************************* FFT.cpp ***************************************************/
/**
 * @file.
 *
 * The file FFT.cpp contains the implementation of interpolation and evaluation functions for univariate 
 * polynomials, implemented using FFT algorithms.
 * 
 * For more information - Read the documentation of the header file FFT.hpp.
 */
  /************************************************************************************************************/
#include <omp.h>
#include "FFT.hpp"
using namespace std;
using namespace NTL;



/**************************************************************************************************
 **************************************************************************************************
 ***                                                                                            ***
 ***   The following are global variables used for integration of multi-threaded FFT with NTL   ***
 ***                                                                                            ***
 **************************************************************************************************
 *************************************************************************************************/

#ifdef _MSC_VER // Visual C++
/* Workaround (shaul): This #ifdef is due to (yet another) violation of the OpenMP standard in the VC++ compiler.
   The threadprivate clause in VC++ is implemented directly as the application of the thread __declspec attribute */
extern "C" __declspec(thread) unsigned int ntlGlobalThreadNum;  /**< Used to pass a unique thread number to NTL for memory management purposes
                                                                    Used by GF2XRegister in NTL/GF2X1.cpp                                                               
                                                                    Defined as extern "C" to avoid name mangling when used with C++.  */
#else
extern "C" unsigned int ntlGlobalThreadNum;                     /**< Used to pass a unique thread number to NTL for memory management purposes
                                                                    Used by GF2XRegister in NTL/GF2X1.cpp                                                               
                                                                    Defined as extern "C" to avoid name mangling when used with C++.  */
#pragma omp threadprivate(ntlGlobalThreadNum)
#endif

extern unsigned int globalFFTNumthreads = 0;     /**< should be set to 0 (automatically sets to number of processors),
                                                      change to a positive integer to force a number of threads in
                                                      GathenGerhardInterpolation() and GathenGerhardEvaluation(). For 
                                                      debug and analysis purposes. */

extern "C" const unsigned int ntlMaxNumThreads;  /**< defines the maximum number of threads supported by multi-threaded version of NTL (change this in GF2X1.cpp) */

/*******************************  end of global multi-threading variables  ***********************/


namespace Algebra {
namespace details {

/******************************************************************************/
/*********************** Interpolation Algorithms *****************************/
/******************************************************************************/

/*****************************************************/
/****** Mateer's FFT interpolation algorithms: *******/
/*****************************************************/

// TODO delete Irr if not used (shaul)
void GathenGerhardInterpolation(UnivariatePolynomial& result, const UniPolyEval_forProof* evaluation, const GF2X& Irr, const Basis& basis, const ExplicitLinearSubset& spanned, const FieldElement& shift) {
	BasisIndex k = basis.getSizeOfBasis();

	///Step 1 - Computing vanishing polynomials for each partial basis:
	VanishingPolynomial* vanishingPolys = new VanishingPolynomial[k+1];
	Basis zeroBasis;
	zeroBasis.addElement(zero());
	vanishingPolys[0] = VanishingPolynomial(zeroBasis, shift);

	Basis partialBasis;
	for (int i=0; i<k; i++) {
		partialBasis.addElement(basis.getBasisElementByIndex(i));	///Expanding the partial basis by another element
		vanishingPolys[i+1] = VanishingPolynomial(partialBasis, shift);	///Computing the vanishing polynomial
	}

	///Step 2 - Calling the recursive FFT function:
	unsigned int numThreads = (globalFFTNumthreads) ? globalFFTNumthreads : omp_get_max_threads();
    numThreads = numThreads < ntlMaxNumThreads ? numThreads : ntlMaxNumThreads;
	int ompNested;
	if (numThreads > 1) {
		omp_set_num_threads(numThreads);
		ompNested = omp_get_nested();
		omp_set_nested(true);
	}

	GathenGerhardInterpolation_Aux(result, evaluation, basis, k-1, 0, vanishingPolys, spanned, shift, numThreads, 0);

	if (numThreads > 1)
		omp_set_nested(ompNested);	//set nested back to original value.

	///Clearing memory:
	delete[] vanishingPolys;
}

void GathenGerhardInterpolation_Aux(UnivariatePolynomial& result, const UniPolyEval_forProof* evaluation, const Basis& basis, const BasisIndex i, const FieldIndex j, 
										VanishingPolynomial* vanishingPolys, const ExplicitLinearSubset& spanned, const FieldElement& shift,
										const unsigned int numThreads, const unsigned int curThreadNum) {
                                        
    ntlGlobalThreadNum = curThreadNum; // needed for NTL memory management in NTL/GF2X1.cpp

	///Base case:
	if (i == -1) {
		result = UnivariatePolynomial(evaluation->queryAtPoint(spanned.getFieldElementByIndex(j)));
		return;
	}

	///Recursive Calls:
	FieldIndex m = ::Infrastructure::POW2(i);
	UnivariatePolynomial res1, res2;
	if(numThreads > 1) { // run each part of the recursion in a different thread. numThreads keeps count of the thread pool while curThreadNum gives a thread ID for NTL memory management in GF2X1.cpp
		#pragma omp parallel sections
		{
			#pragma omp section
			GathenGerhardInterpolation_Aux(res1, evaluation, basis, i-1, j, vanishingPolys, spanned, shift, numThreads/2, curThreadNum);		///res1 = r
			#pragma omp section
			GathenGerhardInterpolation_Aux(res2, evaluation, basis, i-1, j + m, vanishingPolys, spanned, shift, numThreads-(numThreads/2), curThreadNum+(numThreads/2));	///res2 = q*s_i(b_i+1) + r
		} // implied OMP barrier, wait for both threads.
        ntlGlobalThreadNum = curThreadNum; //reset thread number after (implied) barrier.
	} else {
		GathenGerhardInterpolation_Aux(res1, evaluation, basis, i-1, j, vanishingPolys, spanned, shift, 1, curThreadNum);		///res1 = r
		GathenGerhardInterpolation_Aux(res2, evaluation, basis, i-1, j + m, vanishingPolys, spanned, shift, 1, curThreadNum);	///res2 = q*s_i(b_i+1) + r
	}

	UnivariatePolynomial q = (res2 - res1) / UnivariatePolynomial(vanishingPolys[i].evalLinearizedPoly(basis.getBasisElementByIndex(i) + shift));

	result = q.multiplyWithSparse(vanishingPolys[i]) - q*UnivariatePolynomial(vanishingPolys[i].evalLinearizedPoly(spanned.getFieldElementByIndex(j))) + res1;
}

/*****************************************************/
/****** Arnab's FFT interpolation algorithms: ********/
/*****************************************************/
//commented out by Ariel as not used and need some conversions added now that FieldElement typedefd to FieldElement rather than GF2E
//void interpolate_poly(GF2EX& ans, const UniPolynomialEvaluation& values, const Basis& basis, const ExplicitLinearSubset spanned){
//
//	long k = basis.getSizeOfBasis();
//	GF2EX poly, q, p0, p1;
//	vector<FieldElement> newBasis;
//	UniPolynomialEvaluation p0eval, p1eval;
//	FieldElement qu, temp;
//
//	if(k == 1) {
//		NTL::GF2E val_at_0 = ::NTL::GF2E(values.queryAtPoint(zero()));
//		ans =
//				GF2EX(0,val_at_0) +
//				GF2EX(1, (::NTL::GF2E(values.queryAtPoint(basis.getBasisElementByIndex(0)))
//				+ val_at_0) / ::NTL::GF2E(basis.getBasisElementByIndex(0)));
//		return;
//	}
//
//	ans = GF2EX(0,0);
//	NTL::GF2E a_k = ::NTL::GF2E(basis.getBasisElementByIndex(k - 1));
//	q = GF2EX(2,1) + GF2EX(1,a_k);
//	newBasis.SetLength(k-1);
//	for(long i=0; i<(k-1); i++) {
//		NTL::GF2E a_i = ::NTL::GF2E(basis.getBasisElementByIndex(i));
//		newBasis[i] = FieldElement(a_i * (a_i + a_k));
//	}
//
//	//get_span(u, basis);
//	for(FieldIndex i=0; i < spanned.getSizeOfField(); i++){
//		FieldElement elem_i = spanned.getFieldElementByIndex(i);
//		FieldElement elem_i_val = values.queryAtPoint(elem_i);
//
//		qu = elem_i * (elem_i + a_k);
//		temp = (values.queryAtPoint(elem_i + a_k) -	elem_i_val) / a_k;
//		p1eval.addPoint(qu, temp);
//		p0eval.addPoint(qu, elem_i_val - elem_i * temp);
//	}
//
//	ExplicitLinearSubset newSpanned(newBasis);
//	interpolate_poly(p0, p0eval, newBasis, newSpanned);
//	interpolate_poly(p1, p1eval, newBasis, newSpanned);
//
//	return combine_poly(ans, p0, p1, q, power_long(2,k));
//}
//
//void combine_poly(GF2EX& p, const GF2EX& p0, const GF2EX& p1, const GF2EX& q, long n){
//
//	if(n <= 2){
//		p = GF2EX(0,coeff(p0,0)) + GF2EX(1,coeff(p1,0));
//		return;
//	}
//
//	GF2EX a0, b0, a1, b1, a, b;
//	FieldElement c, d;
//
//	for(long i=0; i<n/4; i++){
//		b0 += GF2EX(i, coeff(p0,i));
//		b1 += GF2EX(i, coeff(p1,i));
//	}
//	for(long i=n/4; i<n/2; i++){
//		a0 += GF2EX(i-n/4, coeff(p0,i));
//		a1 += GF2EX(i-n/4, coeff(p1,i));
//	}
//
//	combine_poly(b, b0, b1, q, n/2);
//	combine_poly(a, a0, a1, q, n/2);
//
//	expt2power(d, coeff(q,0),n/4);
//	expt2power(c, coeff(q,1),n/4);
//
//	p = GF2EX(0,0);
//
//	for(long i=0; i<n/4; i++){
//		//SetCoeff(p, i, coeff(a,i)*d + coeff(b,i));
//		p += GF2EX(i, coeff(a,i)*d + coeff(b,i));
//	}
//	for(long i=n/4; i<n/2; i++){
//		//SetCoeff(p, i, coeff(a,i)*d + coeff(a,i-n/4)*c + coeff(b,i));
//		p += GF2EX(i, coeff(a,i)*d + coeff(a,i-n/4)*c + coeff(b,i));
//	}
//	for(long i=n/2; i<3*n/4; i++){
//		//SetCoeff(p, i, coeff(a,i-n/2) + coeff(a,i-n/4)*c);
//		p += GF2EX(i, coeff(a,i-n/2) + coeff(a,i-n/4)*c);
//	}
//	for(long i=3*n/4; i<n; i++){
//		//SetCoeff(p, i, coeff(a,i-n/2));
//		p += GF2EX(i, coeff(a,i-n/2));
//	}
//}
//
//// log(degree) field operations
//// Sets y to x^degree
void expt2power(FieldElement& y, const FieldElement& x, int degree) {
	long i = 1;
	y = x;

	while(i != degree){
		y = FieldElement::sqr(y);
		i *= 2;
	}
}

/******************************************************************************/
/************************* Evaluation Algorithms ******************************/
/******************************************************************************/

/*****************************************************/
/*******  Mateer's FFT Evaluation algorithms:  *******/
/*****************************************************/

/** Mateer's Evaluation FFT: */
void GathenGerhardEvaluation(FieldElement* result_Array, UniPolynomialEvaluation* result_Map, const GF2X& Irr, const UnivariatePolynomial& poly, const Basis& basis, const FieldElement& shift) {
	BasisIndex k = basis.getSizeOfBasis();
	
	///Step 1 - Computing vanishing polynomials for each partial basis:
	VanishingPolynomial* vanishingPolys = new VanishingPolynomial[k+1];
	Basis zeroBasis;
	zeroBasis.addElement(zero());
	vanishingPolys[0] = VanishingPolynomial(zeroBasis, shift);	///the polynomial that vanishes on W0 is P(x)=x.
	
	Basis partialBasis;
	for (int i=0; i<k; i++) {
		partialBasis.addElement(basis.getBasisElementByIndex(i));	///Expanding the partial basis by another element
		if (i < GATHEN_GERHARD_OPTIMAL_BASE) 
			continue;	///No need to compute vanishing polynomial because of the base case

		vanishingPolys[i+1] = VanishingPolynomial(partialBasis, shift);
	}
	
	///Step 2 - Spanning the field:
	ExplicitLinearSubset spanned(basis);

	///Step 3 - Calling the recursive FFT function:
	UnivariatePolynomial quo, rem;
	poly.divideBySparse(vanishingPolys[k], quo, rem);
 
	int numThreads = 1;
	int ompNested;
	if (basis.getSizeOfBasis() >= EVALUATION_MULTITHREADED_FFT_LOWER_BOUND) {
		numThreads = (globalFFTNumthreads) ? globalFFTNumthreads : omp_get_max_threads();
        numThreads = numThreads < static_cast<int>(ntlMaxNumThreads) ? numThreads : static_cast<int>(ntlMaxNumThreads);
		ompNested = omp_get_nested();
		omp_set_nested(true);
	}

	GathenGerhardEvaluation_Aux(result_Array, result_Map, rem, basis, k-1, 0, vanishingPolys, spanned, shift, numThreads, 0);
	
	if (basis.getSizeOfBasis() >= EVALUATION_MULTITHREADED_FFT_LOWER_BOUND) 
		omp_set_nested(ompNested);

	///Clearing memory:
	delete[] vanishingPolys;
}

/** The recursive function that computes the evaluation using Gathen-Gerhard FFT */
void GathenGerhardEvaluation_Aux(FieldElement* result_Array, UniPolynomialEvaluation* result_Map, const UnivariatePolynomial& f, const Basis& basis, const BasisIndex i, const FieldIndex j, 
																		VanishingPolynomial* vanishingPolys, ExplicitLinearSubset& spanned, const FieldElement& shift, const unsigned int numThreads,
                                                                        const unsigned int curThreadNum) {
    ntlGlobalThreadNum = curThreadNum;
	//Base case: evaluating the polynomial and filling in the result:
	if (i == GATHEN_GERHARD_OPTIMAL_BASE) {
		if (result_Array != NULL) {		///The results go an array
			for (int d=0; d < ::Infrastructure::POW2(i+1); d++)
				result_Array[j+d] = f.queryAtPoint(spanned.getFieldElementByIndex(j+d) + shift);
			return;
		}

		else {							///The results go a univariate evaluation represented as a map
			for (int d=0; d < ::Infrastructure::POW2(i+1); d++) {
				FieldElement point = spanned.getFieldElementByIndex(j+d) + shift;
				FieldElement value = f.queryAtPoint(spanned.getFieldElementByIndex(j+d) + shift);
				// All the points across all threads get added to the same map. The map data structure
				// isn't thread-safe, so we use a critical section to protect it.
				//TODO: Add points in big batches to reduce overhead. --Eran
				#pragma omp critical (ggEvalAddPoint)
				{ result_Map->addPoint(point, value); }
			}
			return;
		}
	}

	FieldIndex m = ::Infrastructure::POW2(i);
	assert(f.getDegree() < 2*m);
	UnivariatePolynomial q,r;
	f.divideBySparse(vanishingPolys[i] - vanishingPolys[i].evalLinearizedPoly(spanned.getFieldElementByIndex(j) + shift), q, r);

	UnivariatePolynomial f1 = r;
	UnivariatePolynomial f2 = r + q*UnivariatePolynomial(vanishingPolys[i].evalLinearizedPoly(basis.getBasisElementByIndex(i) + shift));

	if (numThreads > 1) {// run each part of the recursion in a different thread. numThreads keeps count of the thread pool while curThreadNum gives a thread ID for NTL memory management in GF2X1.cpp
		#pragma omp parallel sections
		{
			#pragma omp section
			GathenGerhardEvaluation_Aux(result_Array, result_Map, f1, basis, i-1, j, vanishingPolys, spanned, shift, numThreads/2, curThreadNum);
			#pragma omp section		
            GathenGerhardEvaluation_Aux(result_Array, result_Map, f2, basis, i-1, j + m, vanishingPolys, spanned, shift, numThreads-(numThreads/2), curThreadNum+(numThreads/2));
		}
        ntlGlobalThreadNum = curThreadNum; // reset thread number after implied barrier
	} else {
		GathenGerhardEvaluation_Aux(result_Array, result_Map, f1, basis, i-1, j, vanishingPolys, spanned, shift, 1, curThreadNum);
		GathenGerhardEvaluation_Aux(result_Array, result_Map, f2, basis, i-1, j + m, vanishingPolys, spanned, shift, 1, curThreadNum);
	}
}

/*****************************************************/
/*******  Arnab's FFT Evaluation algorithms:  ********/
/*****************************************************/
//
//void fft_additive(UniPolynomialEvaluation& ans, const UnivariatePolynomial& poly, const Basis& basis){
//	long i,k,n;
//	UnivariatePolynomial q;
//	UniPolynomialEvaluation lans, rans;
//	Basis newBasis;
//	FieldElement elt, lval, rval, currElement, elementK;
//
//	k = basis.getSizeOfBasis();
//	elementK = basis.getBasisElementByIndex(k-1);	//The last element of the basis
//	n = power_long(2,k);
//
//	if(k==1){
//		ans.addPoint(zero(), poly.getCoeff(0));
//		ans.addPoint(basis.getBasisElementByIndex(0), basis.getBasisElementByIndex(0)*poly.getCoeff(1) + poly.getCoeff(0));
//		return;
//	}
//
//	q.setCoeff(2,one());
//	q.setCoeff(1, elementK);
//
//	for (i=0; i<(k-1); i++) {
//		currElement = basis.getBasisElementByIndex(i);
//		newBasis.addElement(currElement * (currElement + elementK));
//	}
//
//	GF2EX p0_rep, p1_rep;
//	find_couple(p0_rep, p1_rep, poly.getPoly(), q.getPoly(), n);
//	UnivariatePolynomial p0(p0_rep);
//	UnivariatePolynomial p1(p1_rep);
//	fft_additive(lans, p0, newBasis);
//	fft_additive(rans, p1, newBasis);
//
//	ExplicitLinearSubset span(basis);
//	
//	for(i=0; i<n; i++){
//		elt = span.getFieldElementByIndex(i);
//		lval = lans.queryAtPoint(elt * (elt + elementK));
//		rval = rans.queryAtPoint(elt * (elt + elementK));
//		ans.addPoint(elt, lval+elt*rval);
//	}
//
//}
//
//// If p(x) is a degree n polynomial and q is a monic quadratic polynomial,
//// find_couple returns <p0,p1> where p(x) = p0(q(x)) + x*p1(q(x)).
//// Assume n is a power of 2.
//void find_couple(GF2EX& p0, GF2EX& p1, const GF2EX& poly, const GF2EX& q, long n) {
//	GF2EX as0, as1, bs0, bs1;
//
//	if(n <= 2){
//		p0 = GF2EX(0, coeff(poly, 0));
//		p1 = GF2EX(0, coeff(poly, 1));
//		return;
//	}
//
//	FieldElement c, d;
//	GF2EX a, b;
//	pair<GF2EX,GF2EX> as, bs;
//	GF2EX alpha, beta, gamma;
//
//	for(long i=0; i<n/2; i++)
//		gamma += GF2EX(i, coeff(poly,i));
//	for(long i=n/2; i<(3*n/4); i++)
//		beta += GF2EX(i-n/2, coeff(poly,i));
//	for(long i=3*n/4; i<n; i++)
//		alpha += GF2EX(i-3*n/4, coeff(poly,i));
//
//	c = power(coeff(q,1),n/4);
//	d = power(coeff(q,0),n/4);
//
//	a = (alpha << (n/4)) + (beta - c*alpha);
//	b = (((c*c-d)*alpha + c*beta) << (n/4)) + c*d*alpha - d*beta + gamma;
//
//	find_couple(as0, as1, a,q,n/2);
//	find_couple(bs0, bs1, b,q,n/2);
//
//	p0 = (as0 << (n/4)) + bs0;
//	p1 = (as1 << (n/4)) + bs1;
//}

} // of namespace details
} // of namespace Algebra

#include <algebraLib/UnivariatePolynomialGeneral.hpp>
#include "common/Infrastructure/Infrastructure.hpp"
#include "common/lightCircLib/lightCircPoly.hpp"
#include "common/Algebra/details/FiniteFields.hpp"
//#include "ACSPWitnessChecker.hpp"
//#include "ACSPConverter.hpp"
#include "../MultiVarPoly.hpp"
#include "algebraLib/SubspacePolynomial.hpp"

using namespace std;
using namespace PCP_Project;
using namespace Algebra;
using namespace Algebra::details;
using namespace NTL;

////testing we can make copies of same ACSPCompositionPolynomial
//TEST(ACSP, copyComposed){
//	initGF2E(64);
//	ntlGlobalThreadNum = 2;
//	/* setting P=X*Y, N=X, A = 1+X.
//	* so Q=X*(A(N(X))=X^2 +X, and Q(x) =x^2 + x
//	*/
//	//creating the polynomial X*Y
//	MultiVarPoly P;
//	P.AddMonomial({ 0, 1 });
//
//	MultiVarPoly P2;
//	P2.AddMonomial({ 0, 1 });
//
//	vector<FieldElement> coeffs = { zero(), one() };
//	unique_ptr<UnivariatePolynomialInterface> N(new UnivariatePolynomialGeneral(coeffs));
//	{
//		ACSPPartialInstance::polynomialsVec neighborPolys;
//		neighborPolys.push_back(move(N));
//
//	
//		const FieldElement x = xElement();
//		FieldElement res = x*x + x;
//#pragma omp parallel sections num_threads(2)
//		{
//#pragma omp section 
//			{
//				int a;
//				for (int i = 0; i < 250000; i++){
//					//#pragma omp critical
//					{
//						P.eval({ x, one() });
//					}
//
//				}
//			}
//#pragma omp section 
//			{
//			int a;
//			for (int i = 0; i < 250000; i++){
//				//#pragma omp critical
//
//				{
//					P2.eval({ one(), one() });
//				}
//
//			}
//		}
//
//		}
//	}
//
//}
