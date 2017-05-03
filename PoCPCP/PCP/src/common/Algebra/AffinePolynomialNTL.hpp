///********************************************** AffinePolynomialNTL.hpp ************************************************/
///**
//* @file.
//*
//
//* A class of (affine) linearized polynomials - i.e.,
//univariate polynomials whose non-zero coefficients are only of monomials of the form  x^{2^i},
//and possibly, a non-zero constant coefficient.
//
//written by Ariel Gabizon ariel.gabizon@gmail.com and Michal Riabzev
//* Almost identical to AffinePolynomial class - but uses only NTL objects like the FieldElement typedef to avoid unnecessary conversions 
//warning:: implementation - especially, of the compuateMat() function, currently depends heavily on specifics of NTL
//*
//*
//*/
///************************************************************************************************************/
//#ifndef _COMMON_REFACTORING_ALGEBRA_AFFINEPOLYNOMIALNTL_HPP_
//#define _COMMON_REFACTORING_ALGEBRA_AFFINEPOLYNOMIALNTL_HPP_
//
//#include <map>
////#include "details/FiniteFields.hpp"
////#include "details/AlgebraUtils.hpp"
////#include "NTL/GF2E.h"
////typedef Algebra::elementsSet_t elementsSet_t;
//
////#include "../constraints/gadgetlib2/variable.hpp"
////#include "details\AlgebraUtils.hpp"
//#include <vector>
//#include <string>
//#include "AlgebraCommon.hpp"
//#include <algebraLib/PolynomialDegree.hpp>
//namespace Algebra {
//
//
//	/**
//	* As described above, a class for general Affine polynomials.
//	*/
//	class AffinePolynomialNTL{
//
//	public:
//
//		/** Constructors **/
//		//Default const is discouraged for use, currently added cause Ariel using it in VanishingPolynomial class;
//		AffinePolynomialNTL(){}
//		AffinePolynomialNTL(const std::vector<FieldElement>& coefficients, const FieldElement& constantFactor);
//
//	
//		/** The function evaluates the Affine polynomial at a given point and returns the result. */
//		FieldElement eval(const FieldElement& x)const;
//
//		//return the i'th coefficient of this polynomial
//		FieldElement getCoefficient(const unsigned int i)const;
//		inline void subtractConstantTerm(const FieldElement& sub) { constantFactor_ -= sub; }
//		PolynomialDegree getDegree() const;
//
//		/** Class Destructor */
//		virtual ~AffinePolynomialNTL(){};
//
//	private:
//		std::vector<FieldElement> coefficients_;	//coefficients[i] = coefficient c of the monomial c*x^(2^i)
//		FieldElement constantFactor_;  //The constant factor of the vanishing polynomial, relevant for affine subspaces and not linear.
//		NTL::mat_GF2 polyMat_; //the matrix corresponding to evaluating this poly (not including adding constantFactor)
//		/**updates the polyMat field to contain correct evaluation matrix*/
//		void computeMat();
//
//
//	};
//
//
//}// of namespace Algebra
//
//#endif //_COMMON_REFACTORING_ALGEBRA_AFFINEPOLYNOMIALNTL_HPP_
