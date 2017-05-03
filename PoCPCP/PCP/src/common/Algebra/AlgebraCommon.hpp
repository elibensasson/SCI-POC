#ifndef _COMMON_REFACTORING_ALGEBRA_ALGEBRACOMMON_HPP_
#define _COMMON_REFACTORING_ALGEBRA_ALGEBRACOMMON_HPP_

#include <gtest/gtest.h>
#include <vector>
#include <set>
#include "details/FiniteFields.hpp"
#include <algebraLib/FieldElement.hpp>
#include <NTL/GF2XFactoring.h>

/* The purpose of this class is to put together various convenient methods related to the algebraic classes.
* In particular, many useful conversions between the older algebraic objects used mainly
* in the PCP code, to newer algebraic objects used a lot in the ACSP and BREX->ACSP code (and perhaps elsewhere)
* Class created by Ariel Gabizon - ariel.gabizon@gmail.com
*/

using Algebra::details::Basis;
using Algebra::FieldElement;
using Algebra::details::FieldIndex;
using NTL::GF2E;
using NTL::vec_GF2;
using std::vector;
using NTL::to_GF2X;
namespace Algebra{

	/*****CONVERSIONS *****/
	//converts a basis of FieldElement to a vector of FieldElement
	vector<FieldElement> oldBasisToElementVector(const Algebra::details::Basis& b);
	//converts an array of FieldElement to a vector of FieldElement
	vector<FieldElement> fieldPointArrayToElementVector(const FieldElement* c, size_t size);

	//converts an array of FieldElement to a vector of FieldElement
	vector<FieldElement> fieldPointArrayToPointVector(const FieldElement* c, size_t size);
	//converts a vector of FElem to a vector of FieldElement
	vector<FieldElement> FElemVectorToFieldElementVector(const vector<FieldElement> v);

	//converts an NTL GF2 Vector to a FieldElement
	inline FieldElement vec_GF2ToFieldElement(const NTL::vec_GF2& src){
		return FieldElement(to_GF2E(to_GF2X(src)));
	}
	

	/** returns the vector/coeffs of polynomial over GF2 corresponding to the field element e
	subtlety : if an element e of GF_{2^m} has only the first l coeffs non-zero NTL will store
	e as a vector of length l. On the other hand, this method always returns a vector of length m. e.g., in the case of
	such e, the vector returned will end with m-l zeros. */
	NTL::vec_GF2 fieldPointToGF2Vector(const FieldElement& e);
		


	/** returns the vector/coeffs of polynomial over GF2 corresponding to the field element e
	subtlety : if an element e of GF_{2^m} has only the first l coeffs non-zero NTL will store
	e as a vector of length l. On the other hand, this method always returns a vector of length m. e.g., in the case of
	such e, the vector returned will end with m-l zeros. */
	NTL::vec_GF2 fieldPointToFullVector(const FieldElement& e);


	/**converts a FieldElement to integer.
	Specifically if e corresponds to GF2 poly b_0 + b_1*x +.. b_s*x^s
	returns integers 1 + 2*b_1 +.. 2^s*b_s.
	Currently assumes e contains at most 32 non-zero bits, i.e., s<32.
	Depends heavily on NTL's representation of elements in buckets of 32 bits.
	*/
	//this method not working right now!
	FieldIndex fieldElementToInteger(const FieldElement& e);

	/***OTHER USEFUL METHODS *****/

	//generates FieldELement x
	FieldElement xElement();

//naive powering-don't use for large powers. Returns a^b.
	FieldElement pow(const FieldElement& a, int b);
	/**(credit:A cut and paste of Michael's method converting FieldElement.)
	* Let "numToWrite" = \f$ 2^0  b_0 + 2^1  b_1 + \dots + 2^k  b_k \f$
	* Lets name "shift" as S , and "numBits" as N
	*
	* This functions generated the field element that
	* as a polynomial over GF(2) might be represented as:
	* \f$ x^{s} b_0 + x^{s+1} b_1 + \dots x^{s+n-1} b_{n-1} \f$
	**/
	//return the field extension degree currently being used by NTL
	inline long getNTLFieldDegree(){
		return GF2E::degree();
	}
	inline void setNTLFieldDegree(size_t n){
		::NTL::GF2X Irr;
		::NTL::BuildIrred(Irr, n);

		NTL::GF2E::init(Irr);

	}

}//namespace Algebra
#endif //_COMMON_REFACTORING_ALGEBRA_ALGEBRACOMMON_HPP_
