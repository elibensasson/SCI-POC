#include "LightUniPolyEval.hpp"
#include "common/Infrastructure/Infrastructure.hpp"

#include "FFT/FFTCommon.hpp"

#include <algebraLib/FFT.hpp>

namespace Algebra{

using namespace details;
using Infrastructure::POW2;

	/** allocates needed size for array according to space size. This version getting basis and shift as parameter rather than  */
	LightUniPolyEval::LightUniPolyEval(const Algebra::details::Basis& basis, const FieldElement& shift) : orderedBasis_(basis.asVector()),spaceShift_(shift){
        evaluationTable.resize(POW2(orderedBasis_.size()));
    }


	/** runs Matan's FFT and stores the result */
	LightUniPolyEval::LightUniPolyEval(const UnivariatePolynomial& poly, const details::Basis& basis, const FieldElement& affineShift){
	    fillEvaluation(poly,basis,affineShift);
    }

	/** runs Matan's FFT and stores the result */
	LightUniPolyEval::LightUniPolyEval(const UnivariatePolynomialInterface& poly, const details::Basis& basis, const FieldElement& affineShift):
        orderedBasis_(basis.asVector()),
        spaceShift_(affineShift),
        evaluationTable(poly.eval(basis.asVector(),affineShift)){};

	/** Fills evaluation using Matan's FFT. Used incase it's `too late' to use analogous constructor*/
	void LightUniPolyEval::fillEvaluation(const UnivariatePolynomial& poly, const details::Basis& basis, const FieldElement& affineShift) {
        
        spaceShift_ = affineShift;
        orderedBasis_ = basis.asVector();
		
		//TODO: how to delete evaluationTable when already allocated but not cause error if it hasn't
		int k = orderedBasis_.size();
		if (k >= MATAN_FFT_LOWER_BOUND){

            std::vector<FieldElement> poly_vec(std::max(0,1+int(poly.getDegree())));
            for(size_t i=0; i<poly_vec.size(); i++){
                poly_vec[i] = poly.getCoeff(i);
            }

            evaluationTable = Algebra::FFT(poly_vec,orderedBasis_,spaceShift_);
		}
		else{// if k is small better use naive evaluation
            evaluationTable.resize(POW2(orderedBasis_.size()));
			for (size_t i = 0; i < evaluationTable.size(); i++){
				FieldElement a = getSpaceElementByIndex(orderedBasis_,spaceShift_,i);
				evaluationTable[i] = poly.queryAtPoint(a);
			}

		}
	}
	/**
	* The function evaluates the polynomial at a given point and returns the result.
	* param x is the field point that the polynomial is evaluated at.
	*/
	const FieldElement& LightUniPolyEval::queryAtPoint(const size_t& i) const {
		return evaluationTable[i];
	}
    
    const FieldElement& LightUniPolyEval::queryAtPoint(const FieldElement& e)const{
        return queryAtPoint(getSpaceIndexOfElement(orderedBasis_,spaceShift_,e));
    }



	//generates an UniPolynomialEvaluation object holding same eval, for consistency with old code parts
	void LightUniPolyEval::fillOldEvaluation(UniPolynomialEvaluation& eval){
		for (size_t i = 0; i < evaluationTable.size(); i++){
            FieldElement a = getSpaceElementByIndex(orderedBasis_,spaceShift_,i);
			eval.addPoint(a, evaluationTable[i]);
		}
	}


	void  LightUniPolyEval::reset(const details::Basis& basis, const FieldElement& affineShift){
		_COMMON_ASSERT((basis.getSizeOfBasis() == orderedBasis_.size()), "can reset only to space of same size- point is to save realocations");
        spaceShift_ = affineShift;
        for(int i=0; i< orderedBasis_.size(); i++){
            orderedBasis_[i] = basis.getBasisElementByIndex(i);
        }
	}

	/**
	* The function adds a point-value pair to the evaluation object.
	* param point is the field point.
	* param value is the value of the polynomial in the point.
	*/
	void  LightUniPolyEval::addPoint(const size_t& pointIndex, const FieldElement& value) {
		_COMMON_ASSERT((pointIndex < evaluationTable.size()), "Out of array bounds");
		evaluationTable[pointIndex] = value;
	}


	
    FieldElement LightUniPolyEval::getPoint(size_t i) const{
        return getSpaceElementByIndex(orderedBasis_,spaceShift_,i);
    }
    
    /**
     * The function adds a point-value pair to the evaluation object.
     * param point is the field point.
     * param value is the value of the polynomial in the point.
     */
    void LightUniPolyEval::addPoint(const FieldElement& point, const FieldElement& value){
        addPoint(getSpaceIndexOfElement(orderedBasis_,spaceShift_,point),value);
    }

    //cloning
    UniPolyEval_forProof* LightUniPolyEval::clone()const{
        return new LightUniPolyEval(*this);
    }



	LightUniPolyEval::~LightUniPolyEval()
	{
	}
}//Namespace Algebra
