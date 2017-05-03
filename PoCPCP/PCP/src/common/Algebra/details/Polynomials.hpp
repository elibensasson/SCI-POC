/********************************************** Polynomials.hpp ************************************************/
/** @file
 *
 * The file Polynomials.hpp contains the implementation of various polynomial - univariate and multivariate, polynomials
 * as coefficients and polynomial evaluations etc.
 * 
 * Main classes defined in the file:
 * UnivariatePolynomial  - Implements a univariate polynomial represented as coefficients.
 * VanishingPolynomial - Implements a univariate polynomial that vanishes on a certain subspace (and thus is very sparse).
 * UniPolynomialEvaluation - Implements an evaluation of a univariate polynomial (represented as a map from element to element).
 * 
 * Main date types defined in the file:
 * (Empty)
 *
 * Main functions defined in the file:
 * UnivariatePolynomial Functions - Arithmetic operations, mult or div by a sparse polynomial, FFT interpolation etc.
 * VanishingPolynomial Functions - Compute a vanishing polynomial, Arithmetic operations etc.
 * UniPolynomialEvaluation Functions - Generate the evaluation using a given polynomial and space, query etc.
 */
/************************************************************************************************************/

/**************************************************************************************************************/
/** @origin
 * Arnon says: Used Arnab�s implementation of FFT algorithms for Univariate evaluation and interpolation.
 * Tech report URL: http://publications.csail.mit.edu/tmp/MIT-CSAIL-TR-2005-051.pdf
 * License: ?
 * Version: There is only one version
 * Author: Arnab Bhattacharyya
 **************************************************************************************************************/

#ifndef POLYNOMIALS_HPP_
#define POLYNOMIALS_HPP_

#include <map>
#include "FiniteFields.hpp"
#include "FFT.hpp" 
#include "../../Infrastructure/Infrastructure.hpp" // TODO - cyclical includes, sort this out
#include "../../CXX11_macros.hpp"
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/map.hpp>
#include <NTL/GF2EX.h>
#include "../AlgebraCommon.hpp"
#include "../AffinePolynomialNTL.hpp"
//#include "PCP/PCP_common.hpp"

using std::vector;
namespace Algebra {
namespace details {
///Forward Declaration:
class VanishingPolynomial;
class UniPolyEval_forProof;

/******************************************************************************/
/********************* UnivariatePolynomial Class *****************************/
/******************************************************************************/

/**
 * The class UnivariatePolynomial implements a univariate polynomial represented as coefficients.
 * @Implementation - The polynomial's representation is NTL's GF2EX type. The code can be found in
 * the NTL library under GF2EX.h.
 */
class UnivariatePolynomial {

private:

	NTL::GF2EX poly;			///NTL's type for polynomials over GF(2) extensions.

	/**
	 * An auxiliary function for findVanishingPolynomial. Uses recursion to fill in the coefficients
	 * array, such that coeffArr[i] is the polynomial's coefficient for x^(2^i).
	 * @param b - The basis that the polynomial vanishes on its span.
	 * @param k - The size of the basis currently handled.
	 * @param coeffArr - The result array.
	 */
	void vanishingPolyAux(const Basis& b, const long k, FieldElement* coeffArr);

	/** A helper function for the function that divides this polynomial by a sparse polynomial and returns the result. */
	static void divideBySparse_Aux(const UnivariatePolynomial& P, const VanishingPolynomial& sparsePoly, const int k, VanishingPolynomial* B_Powers, VanishingPolynomial* Zs_Squares, 
																												UnivariatePolynomial& quo, UnivariatePolynomial& rem);
	
	/** Finds t such that (2^t <= deg < 2^(t+1)) and t>=k */
	static int findBoundingT(const int64_t deg, const int k);	

	/** The function receives a degree, and returns the polynomials P0(x), P1(x) such that: P(x) = P1(x)*x^deg + P0(x).	*/
	void splitByDegree(const int64_t& degree, UnivariatePolynomial& P0, UnivariatePolynomial& P1) const;

public:

	/** The irreducible polynomial that is used to build the field extension. */
	static NTL::GF2X IrreduciblePolynomial;		

	/** The class default C'tor */
	UnivariatePolynomial() {};

	/** Class C'tor - Builds a polynomial with the given constant value. */
	UnivariatePolynomial(const FieldElement& constant);

	/** Class C'tor - Builds a polynomial out of NTL's polynomial */
	UnivariatePolynomial(const NTL::GF2EX& poly);

	/** Conversion from polynomial to NTL::GF2EX */
	_COMMON_CXX11_EXPLICIT operator const NTL::GF2EX&()const {return poly;}

	/** The class's copy C'tor */
	UnivariatePolynomial(const UnivariatePolynomial& copy_from);
	
	/** Basic arithmetic operations: */
	void add(const UnivariatePolynomial& other);	///Addition
	void sub(const UnivariatePolynomial& other);	///Subtraction
	void mul(const UnivariatePolynomial& other);	///Multiplication
	void div(const UnivariatePolynomial& other);	///Division
	void rem(const UnivariatePolynomial& other);	///Remainder
	
	/** 
	 * Computes the quotient and remainder of the polynomial a and b: a = q*b + r 
	 * @param a and b are the two polynomials.
	 * @param q is the quotient, r the remainder.
	 */
	static void divRem(UnivariatePolynomial& q, UnivariatePolynomial& r, const UnivariatePolynomial& a, const UnivariatePolynomial& b);
		
	/**
	 * The function computes a vanishing polynomial over b. 
	 * @param Irr is the irreducible polynomial over GF2 that defines the field extension.
	 * @param b is the Basis for the linear space on which the result polynomial vanishes.
	 * @param affineShift is the shift in case the linear space is shifted.
	 * NOTE - This function is currently used by the old code only (old verifier, tests etc.). Vanishing
	 * polynomials should be represented by the class VanishingPolynomial.
	 */
	void findVanishingPoly(const GF2Poly& Irr, const Basis& b, const FieldElement& affineShift = Algebra::zero());

	/**
	 * The function computes a polynomial Z_X that vanishes on a group of field points X.
	 * Unlike findVanishingPoly, the group of points is not a linear subspace, and therefore the algorithm that
	 * computes Z_X is a simple implementation of Z_X = (x+a_1)*...*(x+a_n).
	 * @param X is the group of field elements.
	 */
	void findNonSparseVanishingPoly(const std::vector<FieldElement>& X);

	/**
	 * The function evaluates the polynomial at a given point and returns the result.
	 * @param x is the field point that the polynomial is evaluated at.
	 */
	FieldElement queryAtPoint(const FieldElement& x) const;

	/**
	 * The function evaluates the SPARSE polynomial at a given point and returns the result.
	 * The function can be used only when the polynomial is sparse, i.e. non-zero in coefficients x,x^2,x^4,...
	 * @Implementation - The function only looks at coefficients of the form x^(2^k).
	 * @param x is the field point that the polynomial is evaluated at.
	 */
	FieldElement evalLinearizedPoly(const FieldElement& x) const;

	/**
	* The function multiplies this polynomial by a sparse polynomial and returns the result.
	* @param sparsePoly is the sparse polynomial we multiply by.
	*/
	UnivariatePolynomial multiplyWithSparse(const VanishingPolynomial& sparsePoly) const;

	/**
	* The function divides this polynomial by a sparse polynomial and returns the result.
	* @param sparsePoly is the sparse polynomial we divide by.
	* @param b is the basis that spans the space on which sparsePoly vanishes.
	* @param quo is the quotient, rem is the remainder.
	*/ 
	void divideBySparse(const VanishingPolynomial& sparsePoly, UnivariatePolynomial& quo, UnivariatePolynomial& rem) const;

	/** 
	* The function that divides this polynomial by a sparse polynomial and returns the result. 
	* @IMPORTANT - This function is slower than the function which uses tailored NTL (DivideBySparse), and therefore is not used.
	* @param sparsePoly is the sparse polynomial we divide by.
	* @param b is the basis that spans the space on which sparsePoly vanishes.
	* @param quo is the quotient, rem is the remainder.
	*/
	void sparseDivisionRecursive(const VanishingPolynomial& sparsePoly, const Basis& b, UnivariatePolynomial& quo, UnivariatePolynomial& rem) const;

	/**
	 * The function returns the i'th coefficient of the polynomial.
	 * @param i - The index.
	 */
	FieldElement getCoeff(const int64_t i) const;
    std::vector<FieldElement> getCoeffsVec()const;

	/**
	 * The function sets the i'th coefficient of the polynomial to a specific value.
	 * @param i - The index.
	 * @param coeff - The required coefficient.
	 */
	void setCoeff(const int64_t i, const FieldElement& coeff);

	/**
	* Sets poly to p.
	*/
	void setPoly(const NTL::GF2EX& p);

	/** Returns this polynomial multiplied by x^n */
	UnivariatePolynomial shiftLeft(const int64_t& n) const;

	/** Returns this polynomial divided by x^n */
	UnivariatePolynomial shiftRight(const int64_t& n) const;

	/** Returns this polynomial modulo x^n */
	UnivariatePolynomial truncate(const int64_t& n) const;
	
	/**
	 * Returns the degree of the polynomial.
	 */
	unsigned long getDegree() const;

	/**
	 * Returns true iff the polynomials are equal.
	 */
	bool equals(const UnivariatePolynomial& other) const;

	/**
	 * Returns the polynomial representation as NTL's GF2EX.
	 */
	NTL::GF2EX getPoly() const;

	

	/** 
	 * Wrapper function for NTL's interpolation 
	 * @param evaluation is the polynomial's evaluation over the space.
	 * @param S is the linear space over which the interpolation is made.
	 * TODO - Move to private.
	 */
	void interpolateNTL(const UniPolyEval_forProof& evaluation, const ExplicitLinearSubset& S);

	/**
	 * The function computes an interpolation using the evaluation table "values", made on the linear space whose basis
	 * is "basis". We have 2 interpolation algorithms:
	 * if k >= INTERPOLATION_FFT_LOWER_BOUND then we use Gathen-Gerhard interpolation algorithm.
	 * Otherwise we use NTL's interpolation algorithm.
	 * @param values is the evaluation of the points.
	 * @param basis is the field's basis.
	 * @param spanned is the linear space spanned by the basis.
	 * @param affineShift is the shift in case the linear space is shifted.
	 */
	void interpolation(const UniPolyEval_forProof& values, const Basis& basis, const ExplicitLinearSubset& spanned, const FieldElement& affineShift = Algebra::zero());
	void interpolation(const UniPolyEval_forProof& values, const std::vector<FieldElement> points);

	/** Prints the polynomial's coefficients (For testing purposes) */
	void print() const;
};

/** Univariate polynomial operators: */
UnivariatePolynomial operator+(const UnivariatePolynomial& a, const UnivariatePolynomial& b);

UnivariatePolynomial operator-(const UnivariatePolynomial& a, const UnivariatePolynomial& b);

UnivariatePolynomial operator*(const UnivariatePolynomial& a, const UnivariatePolynomial& b);

UnivariatePolynomial operator/(const UnivariatePolynomial& a, const UnivariatePolynomial& b);

UnivariatePolynomial operator%(const UnivariatePolynomial& a, const UnivariatePolynomial& b);

/** 
* NTL's plainDivRem tailored for sparse b. 
* We assume that the polynomial b has the form: (c_k)*x^2^k + (c_k-1)*x^2^(k-1) + ... + (c_2)*x^2 + (c_1)*x + c_0
*/
void NTLDivSparse(NTL::GF2EX& q, NTL::GF2EX& r, const NTL::GF2EX& a, const VanishingPolynomial& b);	

/** Comparison operator */
bool operator == (const UnivariatePolynomial & p,const UnivariatePolynomial & q);

/** Returns a primitive polynomial of a specified degree */
GF2Poly findPrimitive(const int degree);	


/******************************************************************************/
/************************ VanishingPolynomial Class ***************************/
/******************************************************************************/
/**
* The class represent univariate polynomial that vanishes on a certain field.
* These polynomial have a unique characteristic - Only the coefficients of the exponents
* x, x^2, x^4, x^8... are not zero.
* These coefficients are all kept in an array.
*/
//#define USE_AFFINEPOLYNOMIALNTL
class VanishingPolynomial {

private:

	FieldElement* coefficients;	//coefficients[i] = coefficient c of the monomial c*x^(2^i)
	unsigned int size;			// = array size
	int expoShift;				//The shift of the exponents. example: shift = -1 => 1, x, x^3, x^7... instead of x, x^2, x^4, x^8 (but array is the same)
	FieldElement constantFactor;  //The constant factor of the vanishing polynomial, relevant for affine subspaces and not linear.
#ifdef USE_AFFINEPOLYNOMIALNTL
	AffinePolynomialNTL p; //The affine linearized polynomial corresponding to this one in case expoShift  = 0. A patch to use AffinePolynomial's faster evaluation
#endif
	/**
	 * The function that computes the coefficients of the vanishing polynomial.
	 * @param b - The basis that the polynomial vanishes on its span.
	 * @param k - The size of the basis currently handled.
	 */
	void computeVanishingPoly(const Basis& b, const int k);

public:

	/** Class Default Constructor */
	VanishingPolynomial();

	/** 
	 * Class Constructor - Computes the polynomial that vanishes over a shift of a linear space spanned by basis. 
	 * @param Irr is the irreducible polynomial over GF2 that defines the field extension.
	 * @param basis is the basis that spans the relevant space.
	 * @param elementShift is non-zero only when the relevant linear space is shifted.
	 */
	VanishingPolynomial(const Basis& basis, const FieldElement& elementShift = Algebra::zero());

	/** Class Constructor - Assignes the given coefficint array to the newly created polynomial. */
	VanishingPolynomial(FieldElement* coefficients, const unsigned long& size, const FieldElement& constantFactor);

	/** Copy C'tor */
	VanishingPolynomial(const VanishingPolynomial& copy_from);

	/** Assignment operator */
	VanishingPolynomial& operator= (const VanishingPolynomial& other);

	/** Allocates space according to the given size. */
	void allocateSize(const unsigned long& size);

	/** Returns the i'th element, i.e. coefficient of x^2^i */
	FieldElement getCoeff(const unsigned long i) const { return coefficients[i]; }

	/** The function sets the (2^i)'th coefficient of the polynomial to a specific value. */
	void setCoeff(const long i, const FieldElement coeff);

	/** Returns the constant factor of the polynomial (relevant for affine subspaces only) */
	FieldElement getConstantFactor() const { return constantFactor; }

	/** Subtracts a constant from the polynomial's constant term */
	void subtractConstantTerm(const FieldElement& sub); 

	/**
	 * The function evaluates the sparse polynomial at a given point and returns the result.
	 * @param x is the field point that the polynomial is evaluated at.
	 */
	FieldElement evalLinearizedPoly(const FieldElement& x) const;
	
	/**
	 * The function composes the sparse polynomial with a given polynomial and returns the result.
	 * @param p is the polynomial to compose with
	 */
	NTL::GF2EX evalLinearizedPoly(const NTL::GF2EX& p) const;

	/** Sets the value of the shift	*/
	void setExpoShift(const int expoShift);

	/** Divides all the coefficients by a given constant value */
	void divideByConst(const FieldElement& value);

	/** 
	* The function returns this polynomial without its last cell in the array 
	* i.e. for Zs(x) = c0*x^1 + c1*x^2 + ... + c_(k-1)*x^2^(k-1) + c_k*x^2^k the function 
	* returns Zs'(x) = c0*x^1 + c1*x^2 + ... + c_(k-1)*x^2^(k-1).
	*/
	VanishingPolynomial cutLastMonomial() const;

	/** Returns this polynomial in power 2 */
	VanishingPolynomial square() const;

	/** Returns the size of the coefficient array */
	int getSize() const { return size; }

	/** Returns the degree of the polynomial */
	unsigned long getDegree() const { return ::Infrastructure::safeConvert(::Infrastructure::POW2(size-1)); }

	/** Testing purposes */
	void print() const;

	/** Class Destructor */
	~VanishingPolynomial();
};

VanishingPolynomial operator-(const VanishingPolynomial& vanishPoly, const FieldElement& sub);

/******************************************************************************/
/********************** UniPolynomialEvaluation Class *************************/
/******************************************************************************/

class UniPolyEval_forProof{
private:
	friend class boost::serialization::access;
	template <class Archive>
	void save(Archive & ar, const unsigned int version)const{
		_COMMON_FATAL("Not implemented");
	}
	template <class Archive>
	void load(Archive & ar, const unsigned version){
		_COMMON_FATAL("Not implemented");
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()

    public:

	/**
	 * The function evaluates the polynomial at a given point and returns the result.
	 * param x is the field point that the polynomial is evaluated at.
	 */
	 virtual const FieldElement& queryAtPoint(const FieldElement& x) const = 0;
	
    /**
	 * The function adds a point-value pair to the evaluation object.
	 * param point is the field point.
	 * param value is the value of the polynomial in the point.
	 */
	 virtual void addPoint(const FieldElement& point, const FieldElement& value) = 0;

	/**
	 * The function returns the number of field points evaluated.
	 */
	 virtual FieldIndex getSize() const=0;

     //cloning
     virtual UniPolyEval_forProof* clone()const = 0;
};

/**
 * The class UniPolynomialEvaluation implements an evaluation of a univariate polynomial.
 * @Implementation - The evaluation is actually an stl map, with keys of type FieldElement
 * (GF2E) and values of type FieldElement (GF2E). The comparison function between GF2E's is
 * implemented in FiniteFields.hpp and uses comparison between the weights of each GF2E.
 */
class UniPolynomialEvaluation : public UniPolyEval_forProof {

private:
	friend class boost::serialization::access;
	template <class Archive>
	void save(Archive & ar, const unsigned int version)const{
		ar & *evaluationTable;
	}
	template <class Archive>
	void load(Archive & ar, const unsigned version){
		ar & *evaluationTable;
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()


	///The evaluation is an stl map from one field point to another.
	std::map<FieldElement, FieldElement, Algebra::classCompElements>* evaluationTable;

public:

    UniPolynomialEvaluation(const std::vector<FieldElement>& eval, const std::vector<FieldElement>& orderedBasis, const FieldElement& shift);
    
	/** The class's Default C'tor */
	UniPolynomialEvaluation();

	/** Copy C'tor */
	UniPolynomialEvaluation(const UniPolynomialEvaluation& copy_from);
    UniPolyEval_forProof* clone()const;

	/**
	 * Constructs an evaluation of a univariate polynomial over a linear space whose basis is given.
	 * We determine which evaluation algorithm to use depending on the size of the field:
	 * For smaller fields (<2^EVALUATION_FFT_LOWER_BOUND) we shall use the naive implementation. Otherwise we use FFT.
	 * @param poly is the univariate polynomial (coefficients) we evaluate.
	 * @param BasisS is the basis of the linear space we're evaluating on.
	 * @param affineShift is the shift in case the linear space is shifted.
	 * @note Michael: tried to use it, cause glib error, don't know why
	 */
	UniPolynomialEvaluation(const UnivariatePolynomial poly, const Basis& BasisS, const FieldElement& affineShift = Algebra::zero());

	/** This function is called by the evaluation C'tor, with the boolean member FFT selected according to parameters.
	 * @param poly is the univariate polynomial represented as coefficients.
	 * @param BasisS is the basis for the linear space we evaluate on.
	 * @param affineShift is the shift in case the linear space is shifted.
	 * @param evalType tells whether we want to use the naive implementation, or the Fast Fourier Tranform, or automatically-chosen
	 * @Note - This function is public only for testing purposes. Will be private in the future.
	 */
	enum EvalType { etNaive, etFFT, etAuto };
	void constructEvaluation(const UnivariatePolynomial& poly, const Basis& BasisS, const FieldElement& affineShift = Algebra::zero(), EvalType evalType = etAuto);


	 inline bool contains(const FieldElement& x) const{
		return (evaluationTable->find(x) != evaluationTable->end());
	}

	/**
	 * The function evaluates the polynomial at a given point and returns the result.
	 * param x is the field point that the polynomial is evaluated at.
	 */
	 const FieldElement& queryAtPoint(const FieldElement& x) const;

	/**
	 * The function adds a point-value pair to the evaluation object.
	 * param point is the field point.
	 * param value is the value of the polynomial in the point.
	 */
	 void addPoint(const FieldElement& point, const FieldElement& value);

	/**
	 * The function returns the number of field points evaluated.
	 */
	 FieldIndex getSize() const;

	/** The class's destructor */
	~UniPolynomialEvaluation();

	/** Testing purposes only. */
	void print() const;

	/** Serialization functions: */
	template <class Archive,class Point, class Value, class Compare>	
	static void mapPointerLoad (Archive & ar, std::map<Point,Value*,Compare> & toLoad){
		toLoad.clear();
		int mapSize;
		ar & mapSize;
		for (int i =0;i<mapSize;i++)
		{
			Point* currentPoint = new Point;
			Algebra::details::UniPolynomialEvaluation currentValue;

			ar & (*currentPoint);
			ar & currentValue;
			Algebra::details::UniPolynomialEvaluation * currentValuePointer= new Algebra::details::UniPolynomialEvaluation;
			(*currentValuePointer->evaluationTable) = (*currentValue.evaluationTable);
			toLoad[*currentPoint]=currentValuePointer;

			int size = currentValue.getSize();
		}
	}
};


} // of namespace details
} // of namespace Algebra


#endif /* POLYNOMIALS_HPP_ */
