#ifndef _ALGEBRA_LIGHTUNIPOLYEVAL_HPP
#define _ALGEBRA_LIGHTUNIPOLYEVAL_HPP

#include "details/Polynomials.hpp"
#include "details/FiniteFields.hpp"

#include <algebraLib/PolynomialInterface.hpp>
/****LightUniPolyEval.hpp***/
/*The purpose of this class is to be a much lighter more efficient version of UniPolynomialEvaluation class,
* that corresponds more closely to what is needed in the RS-PCP recursion, especially after the integration of Matan's FFT.
* UniPolynomialEvaluation uses a map, usually implemented by a tree, supporting dynamic deletions and insertions and queries.
* This object is much weaker: You have to say in advance what subspace you are storing evaluations for and an array of appropriate size is alloacted.
* Furthermore, at least for now, no search queries - you have to provide the index of the evaluation you want.
* Written by Ariel Gabizon - ariel.gabizon@gmail.com
*/

namespace Algebra{
	
    class LightUniPolyEval : public details::UniPolyEval_forProof{
	
    public:
		//inits empty evaluation. Not really useful, not having a default constructor was pain in ACSPCompostionPolynomial
		LightUniPolyEval(){}

		/** allocates needed size for array according to space size. This version getting basis and shift as parameter rather than  */
		LightUniPolyEval(const Algebra::details::Basis& basis, const FieldElement& shift) ;

		/** runs Matan's FFT and stores the result */
		LightUniPolyEval(const details::UnivariatePolynomial& poly, const details::Basis& basis, const FieldElement& affineShift) ;
		
		/** runs Matan's FFT and stores the result */
		LightUniPolyEval(const UnivariatePolynomialInterface& poly, const details::Basis& basis, const FieldElement& affineShift) ;
		

		/** Fills evaluation using Matan's FFT. Used incase it's `too late' to use analogous constructor*/
		void fillEvaluation(const details::UnivariatePolynomial& poly, const details::Basis& basis, const FieldElement& affineShift);

		
		/**
		* The function evaluates the polynomial at a given point and returns the result.
		* param x is the field point that the polynomial is evaluated at.
		*/
		const FieldElement& queryAtPoint(const size_t& i)const ;

		const FieldElement& queryAtPoint(const FieldElement& i)const ;


		void  reset(const details::Basis& basis, const FieldElement& affineShift);

		/**
		* The function adds a point-value pair to the evaluation object.
		* param point is the field point.
		* param value is the value of the polynomial in the point.
		*/
		void  addPoint(const size_t& pointIndex, const FieldElement& value);
		
		inline details::FieldIndex getSize() const{ return evaluationTable.size(); }
		FieldElement getPoint(size_t i) const;//returns the point a, such that P(a) is stores at index i
		//generates an UniPolynomialEvaluation object holding same eval, for consistency with old code parts
		void fillOldEvaluation(details::UniPolynomialEvaluation& eval);
        
        /**
         * The function adds a point-value pair to the evaluation object.
         * param point is the field point.
         * param value is the value of the polynomial in the point.
         */
        void addPoint(const FieldElement& point, const FieldElement& value);

        /**
         * The function returns the number of field points evaluated.
         */

        //cloning
        UniPolyEval_forProof* clone()const;
		
		~LightUniPolyEval();
        const std::vector<FieldElement>& getTable()const{return evaluationTable;}
	private:
		std::vector<FieldElement> evaluationTable;
		
        std::vector<FieldElement> orderedBasis_;
        FieldElement spaceShift_;
	};
}//namespace Algebra
#endif //_ALGEBRA_LIGHTUNIPOLYEVAL_HPP
