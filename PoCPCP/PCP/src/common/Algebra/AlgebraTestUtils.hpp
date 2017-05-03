#ifndef _COMMON_REFACTORING_ALGEBRA_ALGEBRATESTUTILS_HPP_
#define _COMMON_REFACTORING_ALGEBRA_ALGEBRATESTUTILS_HPP_

#include "AlgebraCommon.hpp"
#include "details/Polynomials.hpp"
/***AlgebraTestUtils.hpp -- the purpose of this file is to put together various functions that are useful
for the tests of the algebraic classes,e.g., generating a random set of FieldElement  
created by Ariel Gabizon - ariel.gabizon@gmail.com
**/

namespace Algebra{
	/*******************************************************************************/
	/******************************* Field Elements *********************************/
	/*******************************************************************************/

	static ::Algebra::FieldElement randomFieldElement() {
		NTL::GF2E result;
		int i, repetitions = rand() % 5 + 1;
		for (i = 0; i < repetitions; i++)
			random(result._GF2E__rep, getNTLFieldDegree());
		return FieldElement(result);
	}

	static ::Algebra::FieldElement randomFieldElementNoZero()
	{
		NTL::GF2E result;
		int i, repetitions = rand() % 5 + 1;
		while (true) {
			for (i = 0; i < repetitions; i++)
				random(result._GF2E__rep, getNTLFieldDegree());
			FieldElement res(result);
			if (res != zero())
				return res;
		}
	}

	
	//generates an array of length random FieldElements - assumes c has already allocated space
	static void randomFieldElementArray(const int size, FieldElement* c){
		for (int i = 0; i < size; i++){
			c[i] = randomFieldElement();
		}
	}


	//generates a vector of length size of random FieldElements 
	static vector<FieldElement> randomFieldElementVector(const int size){
		vector<FieldElement> v;
		for (int i = 0; i < size; i++){
			v.push_back(randomFieldElement());
		}
		return v;
	}

	/*******************************************************************************/
	/************************** Univariate Polynomials *****************************/
	/*******************************************************************************/

	///Generates a random univariate polynomial and puts in result.
	static void createRandomUnivariatePoly(int IrrPolyDegree, int64_t polyDegree, ::Algebra::details::UnivariatePolynomial& result) {
		for (int64_t i = 0; i <= polyDegree; i++) {
			result.setCoeff(i, randomFieldElementNoZero());
		}
	}

	/*******************************************************************************/
	/********************************** Basis **************************************/
	/*******************************************************************************/

	static ::Algebra::details::Basis createRandomBasis(int size) {
		int i;
		std::vector<FieldElement> basisElements;
		basisElements.resize(size);

		for (i = 0; i < size; i++) {
			while (true) {
				basisElements[i] = generateRandom();
				if (basisElements[i] != ::Algebra::zero())
					break;
			}
		}
		return ::Algebra::details::Basis(basisElements);
	}







}//namespace Algebra

#endif	//_COMMON_REFACTORING_ALGEBRA_ALGEBRATESTUTILS_HPP_

