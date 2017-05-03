/********************************************** Polynomials.cpp ************************************************/
/**
 * @file.
 *
 * The file Polynomials.cpp contains the implementation of various polynomial - univariate and multivariate, polynomials
 * as coefficients and polynomial evaluations etc.
 * 
 * For more information - Read the documentation of the header file Polynomials.hpp.
 */
  /************************************************************************************************************/
#include <stdio.h>
#include <assert.h>
#include <boost/serialization/access.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "Polynomials.hpp"
#include "common/Algebra/LightUniPolyEval.hpp"
#include "algebraLib/FieldElement.hpp"

using namespace std;
using namespace NTL;

namespace Algebra {
namespace details {

/*****************************************************************************************/
/************************* UnivariatePolynomial Class ************************************/
/*****************************************************************************************/
/*
 * Implementation of UnivariatePolynomial class constructors:
 */
/** The irreducible polynomial that is used to build the field extension. */
GF2X UnivariatePolynomial::IrreduciblePolynomial;


/** Class C'tor - Builds a polynomial with the given constant value. */
UnivariatePolynomial::UnivariatePolynomial(const FieldElement& constant) {
	SetCoeff(this->poly, 0, ::NTL::GF2E(constant));
}

/** The class's copy C'tor */
UnivariatePolynomial::UnivariatePolynomial(const UnivariatePolynomial& copy_from) {
	this->poly = GF2EX(copy_from.poly);
}

/** Class C'tor - Builds a polynomial out of NTL's polynomial */
UnivariatePolynomial::UnivariatePolynomial(const GF2EX& poly) {
	this->poly = poly;
}

/*****************************************************************************************/
/*
 * Implementation of wrapper function for basic arithmetic operations:
 */

void UnivariatePolynomial::add(const UnivariatePolynomial& other) {
	this->poly += other.poly;
}

void UnivariatePolynomial::sub(const UnivariatePolynomial& other) {
	this->poly -= other.poly;
}

void UnivariatePolynomial::mul(const UnivariatePolynomial& other) {
	this->poly *= other.poly;
}

void UnivariatePolynomial::div(const UnivariatePolynomial& other) {
	this->poly /= other.poly;
}

void UnivariatePolynomial::rem(const UnivariatePolynomial& other) {
	this->poly %= other.poly;
}

///q = a/b, r = a%b (Wrapper for NTL's DivRem)
void UnivariatePolynomial::divRem(UnivariatePolynomial& q, UnivariatePolynomial& r, const UnivariatePolynomial& a, const UnivariatePolynomial& b) {
	
	NTL::DivRem(q.poly, r.poly, a.poly, b.poly);
}

/** NTL's plainDivRem tailored for sparse b. */
void NTLDivSparse(GF2EX& q, GF2EX& r, const GF2EX& a, const VanishingPolynomial& b) {
	long da, db, dq, i, j;
	GF2E *qp;
	GF2X *xp;

	GF2E LCInv, t;
	GF2X s;

	da = deg(a);
	db = b.getDegree();

	if (db < 0) Error("GF2EX: division by zero");
	if (da < db) {
		r = a;
		clear(q);
		return;
	}

	GF2XVec x(da + 1, 2*GF2E::WordLength());

	for (i = 0; i <= da; i++)
		x[i] = rep(a.rep[i]);

	xp = x.elts();

	dq = da - db;
	q.rep.SetLength(dq+1);
	qp = q.rep.elts();

	for (i = dq; i >= 0; i--) {
		conv(t, xp[i+db]);
		qp[i] = t;

		for (j = db/2; j >= 0; j/=2) {
			if (j==0) {
				mul(s, rep(t), rep(::NTL::GF2E(b.getConstantFactor())));	//For affine shifts
				add(xp[i+j], xp[i+j], s);
				break;
			}

			mul(s, rep(t), rep(::NTL::GF2E(b.getCoeff((unsigned long)::Infrastructure::Log2(j)))));
			//mul(s, rep(t), rep(b.getCoeff((unsigned long)::Infrastructure::Log2(j*2))));	
			add(xp[i+j], xp[i+j], s);

		///The original NTL code:
		//for (j = db-1; j >= 0; j--) {
		//	mul(s, rep(t), rep(bp[j]));
		//	add(xp[i+j], xp[i+j], s);
		//}
		}
	}

	r.rep.SetLength(db);
	for (i = 0; i < db; i++)
		conv(r.rep[i], xp[i]);
	r.normalize();
}

/*
 * Implementation of arithmetic operators:
 */
UnivariatePolynomial operator+(const UnivariatePolynomial& a, const UnivariatePolynomial& b) {
	UnivariatePolynomial result(a);
	result.add(b);
	return result;
}

UnivariatePolynomial operator-(const UnivariatePolynomial& a, const UnivariatePolynomial& b) {
	UnivariatePolynomial result(a);
	result.sub(b);
	return result;
}

UnivariatePolynomial operator*(const UnivariatePolynomial& a, const UnivariatePolynomial& b) {
	UnivariatePolynomial result(a);
	result.mul(b);
	return result;
}

UnivariatePolynomial operator/(const UnivariatePolynomial& a, const UnivariatePolynomial& b) {
	//TODO - Assert b!=0
	UnivariatePolynomial result(a);
	result.div(b);
	return result;
}

UnivariatePolynomial operator%(const UnivariatePolynomial& a, const UnivariatePolynomial& b) {
	//TODO - Assert b!=0
	UnivariatePolynomial result(a);
	result.rem(b);
	return result;
}

bool operator == (const UnivariatePolynomial & p,const UnivariatePolynomial & q) {
	return p.equals(q);
}

/******************************************************************/
/*** Implementation of univariate polynomial main functionality: **/
/******************************************************************/

/** Returns the i'th coefficient. */
FieldElement UnivariatePolynomial::getCoeff(int64_t i) const {
	return FieldElement(coeff(poly, ::Infrastructure::safeConvert(i)));
}

vector<FieldElement> UnivariatePolynomial::getCoeffsVec()const {
    vector<FieldElement> res(1+getDegree());
    for (size_t i =0; i< res.size(); i++){
        res[i] = getCoeff(i);
    }
    return res;
}

/** Sets the i'th coefficient to coeff */
void UnivariatePolynomial::setCoeff(const int64_t i, const FieldElement& coeff) {
	SetCoeff(poly, ::Infrastructure::safeConvert(i), ::NTL::GF2E(coeff));
}


void UnivariatePolynomial::setPoly(const NTL::GF2EX& p){
	poly = p;
}

/** Returns this polynomial multiplied by x^n */
UnivariatePolynomial UnivariatePolynomial::shiftLeft(const int64_t& n) const {
	UnivariatePolynomial res;
	LeftShift(res.poly, this->poly, ::Infrastructure::safeConvert(n));
	return res;
}

/** Returns this polynomial divided by x^n */
UnivariatePolynomial UnivariatePolynomial::shiftRight(const int64_t& n) const {
	UnivariatePolynomial res;
	RightShift(res.poly, this->poly, ::Infrastructure::safeConvert(n));
	return res;
}

/** Returns this polynomial modulo x^n */
UnivariatePolynomial UnivariatePolynomial::truncate(const int64_t& n) const {
	UnivariatePolynomial res;
	trunc(res.poly, this->poly, ::Infrastructure::safeConvert(n));
	return res;
}

/** Returns the polynomial degree */
unsigned long UnivariatePolynomial::getDegree() const {
	return (deg(poly) == -1) ? 0 : deg(poly);
	//return deg(poly);
}

/**Evaluates the polynomial at a given point and returns the result. */
FieldElement UnivariatePolynomial::queryAtPoint(const FieldElement& x) const {
	return FieldElement(eval(this->poly, ::NTL::GF2E(x)));
}

/**Evaluates the SPARSE polynomial at a given point and returns the result. */
FieldElement UnivariatePolynomial::evalLinearizedPoly(const FieldElement& x) const{
  unsigned long deg = getDegree();
  FieldElement z, y;
  
  for(long i=deg; i>=1; i/=2){
    expt2power(z,x,i);
    y += FieldElement(coeff(this->poly, i)) * z;
  }
  return y;
}

/** The function multiplies this polynomial by a sparse polynomial and returns the result */
UnivariatePolynomial UnivariatePolynomial::multiplyWithSparse(const VanishingPolynomial& sparsePoly) const {

	UnivariatePolynomial res;
	GF2EX originalPoly(this->poly);
	int sparseSize = sparsePoly.getSize();

	FieldElement currentCoeff;
	for (int i=0; i<sparseSize; i++) {
		currentCoeff = sparsePoly.getCoeff(i);

		UnivariatePolynomial temp;
		LeftShift(temp.poly, originalPoly, ::Infrastructure::safeConvert(::Infrastructure::POW2(i)));	///temp = P*x^2^i
		res.poly += (temp.poly * NTL::GF2E(currentCoeff));
	}

	///In case of a vanishing polynomial over an affine subspace, multiply the constant factor as well.
	FieldElement constantFactor = sparsePoly.getConstantFactor();
	if (constantFactor != zero())
		res.poly += originalPoly * NTL::GF2E(constantFactor);

	return res;
}


/** The function divides this polynomial by a sparse polynomial and returns the result.	*/ 
void UnivariatePolynomial::divideBySparse(const VanishingPolynomial& sparsePoly, UnivariatePolynomial& quo, UnivariatePolynomial& rem) const {
	
	///Calling NTL's divrem that was tailored to think of the denominator as sparse.
	NTLDivSparse(quo.poly, rem.poly, this->poly, sparsePoly);
}

/** 
 * The function returns a polynomial that vanishes on the field spanned by the given basis. 
 * NOTE - This function is currently used by the old code only (old verifier, tests etc.). Vanishing
 * polynomials should be represented by the class VanishingPolynomial.
 */
void UnivariatePolynomial::findVanishingPoly(const GF2Poly& Irr, const Basis& b, const FieldElement& affineShift) {

	///Initializing GF2E elements for the irreducible polynomial:
	//GF2E::init(Irr);

	int i;
	long basisSize = b.getSizeOfBasis();
	FieldElement* coeffArray = new FieldElement[basisSize + 1];
	for (i=0; i<basisSize+1; i++) {
		coeffArray[i] = zero();
	}

	///Calling the auxiliary function that fills out the array.
	vanishingPolyAux(b, basisSize, coeffArray);

	///Creating the result polynomial such that coeffArr[i] = coefficient of x^(2^i)
	for (i=0; i<basisSize+1; i++) 
		setCoeff(::Infrastructure::POW2( i), coeffArray[i]);
	
	///For the affine case: Z_H_aff(x) = Z_H(x) + Z_H(aff)
	if (affineShift != zero())	
		setCoeff(0, evalLinearizedPoly(affineShift));

	delete[] coeffArray;
}

/** 
 * The function is an auxiliary function that recursively computes vanishing polynomial over a non-linear space X.
 * Base case: X={a_1}  =>  P(x) = (x+a_1)
 * Recursion: X={a_1,...,a_n}   =>   P(X) = P({a_1...a_n/2}_ * P({a_n/2,...a_n).
 */

//Warning: seems incorrect when n not power of 2. That's why changed findNonSparseVanishingPoly not to use it
GF2EX recursiveFindNonSparseVanishingPoly(const vector<FieldElement>& X, const int64_t first, const int64_t last) {

	//Base case: P = (x+a)
	if (first == last) {
		GF2EX res;
		if (first >= (int64_t)X.size()) {	//Out of boundary - poly = 1F.
			SetCoeff(res, 0, ::NTL::GF2E(one()));
		}
		else {						//Poly = (x+a_i)
			SetCoeff(res, 0, ::NTL::GF2E(X[::Infrastructure::safeConvert(first)]));
			SetCoeff(res, 1, ::NTL::GF2E(one()));
		}
		return res;
	}

	//Recursion: X={a_1,...,a_n}   =>   P(X) = P({a_1...a_n/2}_ * P({a_n/2,...a_n).
	GF2EX res1 = recursiveFindNonSparseVanishingPoly(X, first, (last+first)/2);
	GF2EX res2 = recursiveFindNonSparseVanishingPoly(X, (int64_t)ceil((last+first)/2.0), last);
	return res1 * res2;
}

/**
* The function computes a polynomial Z_X that vanishes on a group of field points X.
* Unlike findVanishingPoly, the group of points is not a linear subspace, and therefore the algorithm that
* computes Z_X is a simple implementation of Z_X = (x+a_1)*...*(x+a_n).
* @param X is the group of field elements.
*/
void UnivariatePolynomial::findNonSparseVanishingPoly(const vector<FieldElement>& X) {
	assert(this->poly.rep.length() == 0);			//Making sure the polynomial was not initialized yet.
	//assert(IsPower2(::Infrastructure::safeConvert(X.size())));

	//if (X.size()==0) {	//No input was read. P(x)=0F
	//	GF2EX res;
	//	poly = res;
	//	return;
	//}

	//poly = recursiveFindNonSparseVanishingPoly(X, 0, X.size()-1);


	//Old (NON EFFICIENT) implementation: brought back to life by Ariel as one above  works incorrectly for n
	// that's not power of 2. Also - disagree in both case n-1 multplications, and here no recursion
	
	//first setting poly to 1
	clear(poly);
	SetCoeff(poly, 0); //sets 0'th coeff to 1
	//init currTerm as X
	GF2EX currTerm;
	SetX(currTerm);
	for (int i=0; i<X.size(); i++) {
		SetCoeff(currTerm, 0, ::NTL::GF2E(X[i]));	//Setting currTerm(x) = (x+a_i).
	poly *= currTerm;
	}
}



/**
 * The main interpolation function - Linear subspaces.
 */
void UnivariatePolynomial::interpolation(const UniPolyEval_forProof& values, const Basis& basis, const ExplicitLinearSubset& spanned, const FieldElement& affineShift) {

	///Choosing between NTL's interpolation and Mateer's FFT algorithm:
	if (basis.getSizeOfBasis() >= INTERPOLATION_FFT_LOWER_BOUND) 	///Using FFT (Gathen-Gerhard)
		GathenGerhardInterpolation(*this, &values, UnivariatePolynomial::IrreduciblePolynomial, basis, spanned, affineShift);
	
	else 			
		interpolateNTL(values, spanned);							///Using NTL's algorithm	
}

/**
 * The main interpolation function - Non-Linear subspaces.
 */
void UnivariatePolynomial::interpolation(const UniPolyEval_forProof& values, const vector<FieldElement> points) {

	FieldIndex size = points.size();
	vec_GF2E x_values, y_values;
	x_values.SetLength(::Infrastructure::safeConvert(size));
	y_values.SetLength(::Infrastructure::safeConvert(size));
	FieldElement currentPoint;

	for(FieldIndex i=0; i<size; i++) {
		currentPoint = points[i];
		x_values[::Infrastructure::safeConvert(i)] = ::NTL::GF2E(currentPoint);
		y_values[::Infrastructure::safeConvert(i)] = ::NTL::GF2E(values.queryAtPoint(currentPoint));
	}

	interpolate(poly, x_values, y_values);
}

/** 
* Returns a primitive polynomial 
* The primitive polynomial is found in a lookup table, according to the requested degree.
*/
GF2Poly findPrimitive(const int degree) {
	GF2Poly res;
	int64_t polyRepresentation;

	//Step 1 - Finding the representation of the polynomial in the lookup table:
	switch(degree) {
	case 0: polyRepresentation = 01; break;
	case 1: polyRepresentation = 02; break;
	case 2: polyRepresentation = 07; break;
	case 3: polyRepresentation = 013; break;
	case 4: polyRepresentation = 023; break;
	case 5: polyRepresentation = 045; break;
	case 6: polyRepresentation = 0103; break;
	case 7: polyRepresentation = 0203; break;
	case 8: polyRepresentation = 0435; break;
	case 9: polyRepresentation = 01021; break;
	case 10: polyRepresentation = 02011; break;
	case 11: polyRepresentation = 04005; break;
	case 12: polyRepresentation = 010123; break;
	case 13: polyRepresentation = 020033; break;
	case 14: polyRepresentation = 042103; break;
	case 15: polyRepresentation = 0100003; break;
	case 16: polyRepresentation = 0210013; break;
	case 17: polyRepresentation = 0400011; break;
	case 18: polyRepresentation = 01000201; break;
	case 19: polyRepresentation = 02000047; break;
	case 20: polyRepresentation = 04000011; break;
	case 21: polyRepresentation = 010000005; break;
	case 22: polyRepresentation = 020000003; break;
	case 23: polyRepresentation = 040000041; break;
	case 24: polyRepresentation = 0100000207; break;
	case 25: polyRepresentation = 0200000011; break;
	case 26: polyRepresentation = 0400000107; break;
	case 27: polyRepresentation = 01000000047; break;
	case 28: polyRepresentation = 02000000011; break;
	case 29: polyRepresentation = 04000000005; break;
	case 30: polyRepresentation = 010040000007; break;
	case 31: polyRepresentation = 020000000011; break;
	case 32: polyRepresentation = 040020000007; break;
	case 33: polyRepresentation = 0100000020001; break;
	case 34: polyRepresentation = 0201000000007; break;
	case 35: polyRepresentation = 0400000000005; break;
	case 36: polyRepresentation = 01000000004001; break;
	case 37: polyRepresentation = 02000000012005; break;
	case 38: polyRepresentation = 04000000000143; break;
	case 39: polyRepresentation = 010000000000021; break;
	case 40: polyRepresentation = 020000012000005; break;
    case 41: polyRepresentation = 0x20000000009; break;
    case 42: polyRepresentation = 0x404c0000001; break;
    case 43: polyRepresentation = 0x80008400021; break;
    case 44: polyRepresentation = 0x108800040001; break;
    case 45: polyRepresentation = 0x208010000011; break;
    case 46: polyRepresentation = 0x410080040001; break;
    case 47: polyRepresentation = 0x800000000021; break;
    case 48: polyRepresentation = 0x1000000080203; break;
    case 49: polyRepresentation = 0x2000000000201; break;
    case 50: polyRepresentation = 0x4000480020001; break;
    case 51: polyRepresentation = 0x8400001008001; break;
    case 52: polyRepresentation = 0x10000000000009; break;
    case 53: polyRepresentation = 0x24020000100001; break;
    case 54: polyRepresentation = 0x62000020000001; break;
    case 55: polyRepresentation = 0x80000001000001; break;
    case 56: polyRepresentation = 0x100028020000001; break;
    case 57: polyRepresentation = 0x200000000000081; break;
    case 58: polyRepresentation = 0x400000000080001; break;
    case 59: polyRepresentation = 0x840400004000001; break;
    case 60: polyRepresentation = 0x1000000000000003; break;

	default: _COMMON_FATAL("No primitive polynomial found for the requested degree");
	}

	//Step 2 - Building the primitive polynomial according to the binary representation.
	for (int digitIndex = 0; polyRepresentation != 0; digitIndex++, polyRepresentation >>= 1) {
		int digit = polyRepresentation & 1;
		SetCoeff(res, digitIndex,digit);
	}

	return res;
}

/******************************************************************************/
/** Implementation of univariate polynomial helper and auxiliary functions: ***/
/******************************************************************************/

/** Finds t such that (2^t <= deg < 2^(t+1)) and t>=k */
int UnivariatePolynomial::findBoundingT(const int64_t deg, const int k) {

	int t = k;
	int64_t tExponent = ::Infrastructure::POW2(t);
	while (tExponent < deg) {
		t++;
		tExponent *= 2;
	}
	t -= 1;
	if (deg == ::Infrastructure::POW2(t+1))	t += 1;

	return t;
}


/** The recursive function that divides this polynomial by a sparse polynomial and returns the result; */
void UnivariatePolynomial::divideBySparse_Aux(const UnivariatePolynomial& P, const VanishingPolynomial& sparsePoly, const int k, VanishingPolynomial* B_Powers, VanishingPolynomial* Zs_Squares, 
																												UnivariatePolynomial& quo, UnivariatePolynomial& rem) {
	unsigned long degP = P.getDegree();

	///Base case - the degree of the nominator is smaller than the denominator's.
	if (degP < (unsigned long)::Infrastructure::POW2(k)) {
		///quo is automatically set to 0. 
		rem.poly = P.poly;
		return;
	}

	int t = findBoundingT(degP, k);
	
	///Finding P0 and P1 such that P(x) = P1(x)*x^2^t + P0(x).
	UnivariatePolynomial P0, P1;
	P.splitByDegree(::Infrastructure::POW2(t), P0, P1);

	///Computing temoQuo = P1(x) * Zs * Zs^2 * Zs^4 * ... * Zs^2^(t-k-1):
	UnivariatePolynomial tempQuo = P1;
	for (int i=0; i <= t-k-1; i++) {
		tempQuo = tempQuo.multiplyWithSparse(Zs_Squares[i]);
	}

	///Making the recursive call:
	UnivariatePolynomial resQuo, resRem;
	divideBySparse_Aux(P1.multiplyWithSparse(B_Powers[t-k]) + P0, sparsePoly, k, B_Powers, Zs_Squares, resQuo, resRem);
	quo = resQuo + tempQuo;
	rem = resRem;
}

/** The function divides this polynomial by a sparse polynomial and returns the result */
void UnivariatePolynomial::sparseDivisionRecursive(const VanishingPolynomial& sparsePoly, const Basis& b, UnivariatePolynomial& quo, UnivariatePolynomial& rem) const {

	unsigned long degP = this->getDegree();
	int k = b.getSizeOfBasis();
	UnivariatePolynomial oldP(*this);

	///Edge case - the degree of the nominator is smaller than the denominator's.
	if (degP < (unsigned long)::Infrastructure::POW2(k)) {
		///quo is automatically set to 0. 
		rem.poly = this->poly;
		return;
	}

	///Finding the number t such that 2^t <= deg(P) < 2^(t+1)
	int t = findBoundingT(degP, k);

	///Computing B, B^2, B^4, ... , B^(t-k) when Zs = x^2^^k + B
	///Similarly, computing Zs, Zs^2, Zs^4, ... , Zs^(t-k).
	VanishingPolynomial* B_Powers = new VanishingPolynomial[t-k+1];
	VanishingPolynomial* Zs_Squares = new VanishingPolynomial[t-k+1];
	Zs_Squares[0] = sparsePoly;
	B_Powers[0] = sparsePoly.cutLastMonomial();		///Zs = x^2^k + B

	for (int i=1; i <= t-k; i++) {		
		Zs_Squares[i] = Zs_Squares[i-1].square();	
		B_Powers[i] = Zs_Squares[i].cutLastMonomial();//B_Powers[i-1].square();
	}

	///Calling the recursive function that computes quotient and remainder:
	divideBySparse_Aux(*this, sparsePoly, k, B_Powers, Zs_Squares, quo, rem);

	assert(oldP == *this);	///Remove!

	///Cleaning up:
	delete[] B_Powers;
	delete[] Zs_Squares;
}

/** The function receives a degree, and returns the polynomials P0(x), P1(x) such that: P(x) = P1(x)*x^deg + P0(x). */
void UnivariatePolynomial::splitByDegree(const int64_t& degree, UnivariatePolynomial& P0, UnivariatePolynomial& P1) const {

	UnivariatePolynomial X_Deg;		/// = x^degree
	X_Deg.setCoeff(degree, one());

	trunc(P0.poly, this->poly, ::Infrastructure::safeConvert(degree));					///P0 is now P modulo x^deg
	RightShift(P1.poly, this->poly, ::Infrastructure::safeConvert(degree));			///P1 is now P / X^deg
	//P1.poly = (this->poly - P0.poly) / X_Deg.poly;
}

/** Compare two univariate polynomials */
bool UnivariatePolynomial::equals(const UnivariatePolynomial& other) const {
	return (this->poly == other.poly) ? true : false;
}

/** Returns the polynomial represented as NTL's GF2EX */
GF2EX UnivariatePolynomial::getPoly() const {
	return this->poly;
}

/** An auxiliary function for findVanishingPoly that computes the array of coefficients of the subspace polynomial vanishing on. 
 * the subspace with basis b={e1,..ek}. The algorithm inductively computes the (coeffs of the) subspace polynomial P_i of the subspace with basis {e_1,..,e_i}. 
 * for i=1, this is the polynomial x^2 + e1*x. Assume we have computed P_{i-1}. 
 * it turns out that P_i is P_{i-1} composed from the outside with x^2 + P_{i-1}(e_i)*x
 * this is the formula the code implements.
 */

void UnivariatePolynomial::vanishingPolyAux(const Basis& b, const long k, FieldElement* coeffArr) {
	int i;

	///Base case: result = x^2 + e1*x
	if (k == 1) {
		coeffArr[0] = b.getBasisElementByIndex(0);
		coeffArr[1] = one();
		return;
	}

	///The array will be filled by the inner call
	vanishingPolyAux(b, k-1, coeffArr);

	///Calculating C' = Sigma(di*ek^(2^i)) + ek^(2^(k-1))
	FieldElement ek = b.getBasisElementByIndex(k-1);	///Basis indices start with 0
	FieldElement cPrime = Algebra::zero();
    FieldElement ekPower = ek;
	for (i=0; i<=(k-2); i++) {
		cPrime += coeffArr[i]*ekPower;
		ekPower = power(ekPower,2);
	}
	cPrime += ekPower;

	///Filling out the new array according to the old one:
	coeffArr[k] = one();
	for (i=k-1; i>=1; i--) {
		coeffArr[i] = cPrime*coeffArr[i] + power(coeffArr[i-1], 2);
	}
	coeffArr[0] = cPrime*coeffArr[0];
}

/** Prints the polynomial's coefficients - Testing purpooses only */
void UnivariatePolynomial::print() const {
	cout << "Polynomial: " << this->poly << endl;
}

/** Wrapper function for NTL's interpolation. */
void UnivariatePolynomial::interpolateNTL(const UniPolyEval_forProof& evaluation, const ExplicitLinearSubset& S) {
	vec_GF2E x_values, y_values;
	FieldIndex i;		
	FieldIndex size = S.getSizeOfField();
	FieldElement currentPoint;
	x_values.SetLength(::Infrastructure::safeConvert(size));
	y_values.SetLength(::Infrastructure::safeConvert(size));

	for(i=0; i<size; i++) {
		currentPoint = S.getFieldElementByIndex(i);
		x_values[::Infrastructure::safeConvert(i)] = ::NTL::GF2E(currentPoint);
		y_values[::Infrastructure::safeConvert(i)] = ::NTL::GF2E(evaluation.queryAtPoint(currentPoint));
		//cout << "Point " << i << ": " << x_values[::Infrastructure::safeConvert(i)] << endl;
		//cout << "Value " << i << ": " << y_values[::Infrastructure::safeConvert(i)] << endl;
	}

	interpolate(poly, x_values, y_values);
}

/******************************************************************************************************/
/************************************ VanishingPolynomial Class ***************************************/
/******************************************************************************************************/

/** Class Default Constructor */
VanishingPolynomial::VanishingPolynomial(){
	this->size = 0;
	this->coefficients = NULL;
	this->expoShift = 0;
	this->constantFactor = zero();
}

/** Class Constructor - Computes the polynomial that vanishes over a shifted space spanned by basis. */
VanishingPolynomial::VanishingPolynomial(const Basis& basis, const FieldElement& elementShift) {
	this->size = basis.getSizeOfBasis() + 1;	///For coefficients 2^0, 2^1, ... , 2^k
	this->coefficients = new FieldElement[size];
    for(size_t i=0; i<size;i++) coefficients[i] = Algebra::zero();
	this->constantFactor = zero();//Ariel added for safety. In release mode this might not be default

	///Computing the basis b+shift and then calling computeVanishingPoly on the new basis:
	computeVanishingPoly(basis, basis.getSizeOfBasis());
	this->expoShift = 0;
	this->constantFactor = (elementShift == zero()) ? elementShift : evalLinearizedPoly(elementShift);//ARIEL:Code assumes constantFactor is 0 up to this line
#ifdef USE_AFFINEPOLYNOMIALNTL
	p = AffinePolynomialNTL(fieldPointArrayToPointVector(coefficients,size), constantFactor);
#endif
}

/** Allocates space according to the given size. */
void VanishingPolynomial::allocateSize(const unsigned long& size) {
	this->coefficients = new FieldElement[size];
    for(size_t i=0; i<size;i++) coefficients[i] = Algebra::zero();
	this->size = size;
}

/** Class Constructor - Assignes the given coefficint array to the newly created polynomial. */
VanishingPolynomial::VanishingPolynomial(FieldElement* coefficients, const unsigned long& size, const FieldElement& constantFactor) {
	this->coefficients = coefficients;
	this->size = size;
	this->expoShift = 0;
	this->constantFactor = constantFactor;
}

/** Copy C'tor */
VanishingPolynomial::VanishingPolynomial(const VanishingPolynomial& copy_from) {

	this->size = copy_from.size;
	this->coefficients = new FieldElement[size];
	for (unsigned long i=0; i<size; i++) {
		coefficients[i] = copy_from.coefficients[i];
	}
	this->expoShift = copy_from.expoShift;
	this->constantFactor = copy_from.constantFactor;
#ifdef USE_AFFINEPOLYNOMIALNTL
	this->p = copy_from.p;
#endif
}

/** Assignment operator */
VanishingPolynomial& VanishingPolynomial::operator= (const VanishingPolynomial& other)  {

	if (this != &other) // protect against invalid self-assignment
	{
		// 1: allocate new memory and copy the elements
		FieldElement* newCoefficient = new FieldElement[other.size];
		for (unsigned long i=0; i<other.size; i++)
			newCoefficient[i] = other.coefficients[i];

		// 2: deallocate old memory
		delete[] coefficients;

		// 3: assign the new memory to the object
		coefficients = newCoefficient;
		size = other.size;
		expoShift = other.expoShift;
		constantFactor = other.constantFactor;
#ifdef USE_AFFINEPOLYNOMIALNTL
		p = other.p;
#endif

	}
	// by convention, always return *this
	return *this;
}

/** The function that computes the coefficients of the vanishing polynomial. */
void VanishingPolynomial::computeVanishingPoly(const Basis& b, const int k) {
	int i;

	///Base case: result = x^2 + e1*x
	if (k == 1) {
		if(b.getBasisElementByIndex(0) == zero()) {
			coefficients[0] = one();
			return;
		}
		coefficients[0] = b.getBasisElementByIndex(0);
		coefficients[1] = one();
		return;
	}

	///The array will be filled by the inner call
	computeVanishingPoly(b, k-1);

	///Calculating C' = Sigma(di*ek^(2^i)) + ek^(2^(k-1))
	FieldElement ek = b.getBasisElementByIndex(k-1);	///Basis indices start with 0
	FieldElement cPrime = Algebra::zero();
    FieldElement ekPower = ek;
	for (i=0; i<=(k-2); i++) {
		cPrime += coefficients[i]*ekPower;
		ekPower = power(ekPower,2);
	}
	cPrime += ekPower;

	///Filling out the new array according to the old one:
	coefficients[k] = one();
	for (i=k-1; i>=1; i--) {
		coefficients[i] = cPrime*coefficients[i] + power(coefficients[i-1], 2);
	}
	coefficients[0] = cPrime*coefficients[0];
}

/** The function sets the (2^i)'th coefficient of the polynomial to a specific value. */
void VanishingPolynomial::setCoeff(const long i, const FieldElement coeff) {
#ifdef USE_AFFINEPOLYNOMIALNTL
	_COMMON_FATAL("tell Ariel He needs to update VanishingPolynomial and AffinePolynomialNTL or comment out USE_AFFINEPOLYNOMIALNTL");
#endif
	coefficients[i] = coeff;
}


/** The function evaluates the sparse polynomial at a given point and returns the result. */
FieldElement VanishingPolynomial::evalLinearizedPoly(const FieldElement& x) const {

	unsigned long expo = 1;
	FieldElement res = zero();
	//If it's an (affine) linearized polynomial use poly for fast evaluation
#ifdef USE_AFFINEPOLYNOMIALNTL
	if (expoShift == 0)
		return p.eval(x);
#endif
	//The general case - use NTL's power function, not very efficient
	if (expoShift!=0 && expoShift!=-1) {
		for (unsigned long i=0; i<size; i++) {
			res += coefficients[i] * power(x, expo + expoShift);
			expo *= 2;
		}
		return res + constantFactor;
	}

	//Special case - x=0F and expoShift=-1  =>  Return the constant coefficient:
	if (expoShift==-1 && x==zero())
		return coefficients[0] + constantFactor;

	//Optimization for expoShift=0 or expoShift=-1 which does not use NTL's power:
	FieldElement currPower, power_x = x;
	for (unsigned long i=0; i<size; i++) {
		currPower = coefficients[i] * power_x;
		if (expoShift==-1)
			currPower /= x;

		res += currPower;
		power_x = FieldElement::sqr(power_x);
	}

	return res + constantFactor;
}

/** The function composes the sparse polynomial with a given polynomial and returns the result. */
NTL::GF2EX VanishingPolynomial::evalLinearizedPoly(const NTL::GF2EX& p) const {

	unsigned long expo = 1;
	NTL::GF2EX res;

	//The general case - use NTL's power function, not very efficient
	for (unsigned long i=0; i<size; i++) {
		res += NTL::GF2E(coefficients[i]) * power(p, expo + expoShift);
		expo *= 2;
	}
	return res + NTL::GF2E(constantFactor);
}

/** Sets the value of the shift	*/
void VanishingPolynomial::setExpoShift(const int expoShift) {
	this->expoShift = expoShift;
}

/** Divides all the coefficients by a given constant value */
void VanishingPolynomial::divideByConst(const FieldElement& value) {
#ifdef USE_AFFINEPOLYNOMIALNTL
	_COMMON_FATAL("tell Ariel He needs to update VanishingPolynomial and AffinePolynomialNTL or comment out USE_AFFINEPOLYNOMIALNTL");
#endif

	for (unsigned long i=0; i<size; i++) 
		coefficients[i] /= value;

	constantFactor /= value;
}

/** Subtracts a constant from the polynomial's constant term */
VanishingPolynomial operator-(const VanishingPolynomial& vanishPoly, const FieldElement& term) {
	VanishingPolynomial result(vanishPoly);
	result.subtractConstantTerm(term);
	return result;
}

/** Subtracts a constant from the polynomial's constant term */
void VanishingPolynomial::subtractConstantTerm(const FieldElement& sub) { constantFactor -= sub; 
#ifdef USE_AFFINEPOLYNOMIALNTL
p.subtractConstantTerm(sub);
#endif 
}
/** 
* The function returns this polynomial without its last cell in the array 
* i.e. for Zs(x) = c0*x^1 + c1*x^2 + ... + c_(k-1)*x^2^(k-1) + c_k*x^2^k the function 
* returns Zs'(x) = c0*x^1 + c1*x^2 + ... + c_(k-1)*x^2^(k-1).
*/
VanishingPolynomial VanishingPolynomial::cutLastMonomial() const {
#ifdef USE_AFFINEPOLYNOMIALNTL
	_COMMON_FATAL("tell Ariel He needs to update VanishingPolynomial and AffinePolynomialNTL or comment out USE_AFFINEPOLYNOMIALNTL");
#endif

	int newSize = this->size - 1;
	FieldElement* newCoefficients = new FieldElement[newSize];
	for (int i=0; i<newSize; i++) {
		newCoefficients[i] = FieldElement(this->coefficients[i]);
	}

	return VanishingPolynomial(newCoefficients, newSize, constantFactor);
}

/** Returns this polynomial in power 2 */
VanishingPolynomial VanishingPolynomial::square() const {
#ifdef USE_AFFINEPOLYNOMIALNTL
	_COMMON_FATAL("tell Ariel He needs to update VanishingPolynomial and AffinePolynomialNTL or comment out USE_AFFINEPOLYNOMIALNTL");
#endif

	int newSize = size + 1;						
	VanishingPolynomial res;					
	res.allocateSize(newSize);					
	res.setCoeff(0, zero());				
	for (int i=1; i<newSize; i++) {
		res.setCoeff(i, FieldElement::sqr(this->coefficients[i-1]));
	}

	res.constantFactor = power(this->constantFactor, 2);	

	return res;
}

/** Testing purposes */
void VanishingPolynomial::print() const {
	cout << "Vanishing Polynomial: " << endl;

	if (coefficients == NULL) {
		cout << "Empty object!" << endl;
		return;
	}

	for (unsigned long i=0; i<size; i++) {
		cout << i << ": " << coefficients[i] << endl;
	}

	if (constantFactor != zero())
		cout << "Constant Factor: " << constantFactor << endl;
}

/** Class Destructor */
VanishingPolynomial::~VanishingPolynomial() {
	delete[] coefficients;
}

/***********************************************************************************************************/
/********************************** UniPolynomialEvaluation Class ******************************************/
/***********************************************************************************************************/
/*
 * Class constructors:
 */
UniPolynomialEvaluation::UniPolynomialEvaluation() {
	evaluationTable = new map<FieldElement, FieldElement, Algebra::classCompElements>();
}

UniPolynomialEvaluation::UniPolynomialEvaluation(const UniPolynomialEvaluation& copy_from) {
	evaluationTable = new map<FieldElement, FieldElement, Algebra::classCompElements>(*copy_from.evaluationTable);
}

UniPolynomialEvaluation::UniPolynomialEvaluation(const std::vector<FieldElement>& eval, const std::vector<FieldElement>& orderedBasis, const FieldElement& shift){
    evaluationTable = new map<FieldElement, FieldElement, Algebra::classCompElements>();
    for(size_t i=0; i< eval.size() ; i++){
        const FieldElement x = getSpaceElementByIndex(orderedBasis,shift,i);
        addPoint(x ,eval[i]);
    }
}

UniPolynomialEvaluation::UniPolynomialEvaluation(const UnivariatePolynomial poly, const Basis& BasisS, const FieldElement& affineShift) {
	evaluationTable = new map<FieldElement, FieldElement, Algebra::classCompElements>();
	this->constructEvaluation(poly, BasisS, affineShift, etAuto);
}

void UniPolynomialEvaluation::constructEvaluation(const UnivariatePolynomial& poly, const Basis& BasisS, const FieldElement& affineShift, EvalType evalType) {
		LightUniPolyEval lightEval(poly, BasisS, affineShift);
		lightEval.fillOldEvaluation(*this);
}

/** Class Destructor */
UniPolynomialEvaluation::~UniPolynomialEvaluation() {
	delete evaluationTable;
}

/*****************************************************************************************/
/*
 * Implementation of univariate polynomial evaluation functionality:
 */
const FieldElement& UniPolynomialEvaluation::queryAtPoint(const FieldElement& x) const {
	const auto iter = evaluationTable->find(x);
    if(iter == evaluationTable->end()){
		cout << "query point: " << x << "is missing" << endl;
        throw("unfound");
    }
    return iter->second;
}

/** Adds a (point, value) pair to the evaluation */
void UniPolynomialEvaluation::addPoint(const FieldElement& point, const FieldElement& value) {
	(*evaluationTable)[point] = value;
}


UniPolyEval_forProof* UniPolynomialEvaluation::clone()const{
    return new UniPolynomialEvaluation(*this);
}


FieldIndex UniPolynomialEvaluation::getSize() const {
	return (FieldIndex)evaluationTable->size();
}

/** Testing purposes only - Prints the entire evaluation */
void UniPolynomialEvaluation::print() const {
	map<FieldElement, FieldElement, Algebra::classCompElements>::iterator p;

	cout << endl << "Printing Univariate Evaluation of size: " << getSize() << endl;
	for(p = evaluationTable->begin(); p != evaluationTable->end(); p++) {
		cout << "(" << p->first << "," << p->second << ")" << endl;
	}
}


} // of namespace details
} // of namespace Algebra
