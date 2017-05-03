/********************************************** AffinePolynomial.cpp ************************************************/
/**
* @file.
*
* A class of (affine) A polynomials - i.e.,
univariate polynomials whose non-zero coefficients_ are only of monomials of the form  x^{2^i},
and possibly, a non-zero constant coefficient.
*
*
*/
/************************************************************************************************************/

#include "algebraLib/AffinePolynomial.hpp"

#include <stdio.h>
#include <assert.h>
#include <math.h>

using namespace std;
using namespace NTL;
namespace Algebra {
    
	/******************************************************************************************************/
	/***************************************** Helper functions *******************************************/
	/******************************************************************************************************/
    
	//TODO:Perhaps move these to AlgebraCommon.hpp (some might exist under different name)
	
	NTL::vec_GF2 to_vec(const FieldElement& src){
            using NTL::GF2E;
            using NTL::to_vec_GF2;
            
            return to_vec_GF2(GF2E(src).LoopHole());
    }

    FieldElement to_FieldElement(const NTL::vec_GF2& src){
            using NTL::to_GF2X;
            using NTL::to_GF2E;
            
            return FieldElement(to_GF2E(to_GF2X(src)));
    }


	long getFieldDegree(){
		return GF2E(zero()).degree();
	}


	/** returns the vector/coeffs of polynomial over GF2 corresponding to the field element e
	subtlety : if an element e of GF_{2^m} has only the first l coeffs non-zero NTL will store
	e as a vector of length l. On the other hand, this method always returns a vector of length m. e.g., in the case of
	such e, the vector returned will end with m-l zeros. */
	NTL::vec_GF2 getElementRep(const FieldElement& e){
		NTL::GF2E f = NTL::GF2E(e);
		NTL::vec_GF2 a = to_vec_GF2(f.LoopHole());
		NTL::vec_GF2 b(NTL::INIT_SIZE, f.degree());
		for (long i = 0; i < a.rep.length(); i++)
			b.rep[i]= a.rep[i];

		return b;
	}





	/* The function evaluates the linear part of the Affine polynomial at a given point and returns the result.
	It is used in the computeMat method to evaluate the linear part of the 
	poly on a basis in order to construct the matrix corresponding
	to the poly's operation.
	After the affine poly is intialized, this matrix will be used to evaluate it.
	Tests showed using the matrix, rather than regular univariate evaluation using coeffs,  can make the evaluation
	16 times faster!
	*/
	FieldElement evalLinearPart(const FieldElement& x, const vector<FieldElement>& coefficients)  {

		unsigned long expo = 1;
		FieldElement res = zero();

		/*we take advantage of the fact squaring is very efficient in GF2E,
		and that x^{2^i} = (x^{2^{i-1})^2. Each iteration we square the current power
		of x- power_x, multiply it by the relevant coefficient, and add to the sum,
		resulting in sum_{i=0}^{size-1} coeffcieints[i]*x^{2^i} */

		FieldElement currPower, power_x = x;
		for (unsigned long i = 0; i<coefficients.size(); i++) {
			currPower = coefficients[i] * power_x;

			res += currPower;
			power_x = FieldElement(sqr(NTL::GF2E(power_x)));
		}

		return res;
	}


    

	/******************************************************************************************************/
	/************************************ AffinePolynomial Class ******************************************/
	/******************************************************************************************************/
	
	/** Class Constructor - Assigns the given coefficient array to the newly created polynomial. */
	AffinePolynomial::AffinePolynomial(const vector<FieldElement>& coefficients, const FieldElement& constantFactor) {
		if (coefficients.size())
		{

			//
			//Initialize coefficients 
			//
			const FieldElement ZERO = zero();

			//find highest non zero index:
			int lastIndex = coefficients.size() - 1;
			while ((lastIndex >= 0) && (coefficients[lastIndex] == ZERO)){
				lastIndex--;
			}
			coefficients_ = vector<FieldElement>(coefficients.begin(), coefficients.begin() + lastIndex + 1);
		}
	
        constantFactor_ = constantFactor;
		computeMat();
	}
    
    unique_ptr<PolynomialInterface> AffinePolynomial::clone()const{
        return unique_ptr<PolynomialInterface>(new AffinePolynomial(coefficients_, constantFactor_));
    }
		
    /** The function evaluates the Affine polynomial at a given point and returns the result. */
    FieldElement AffinePolynomial::eval(const Algebra::FieldElement& x)const{
	     return to_FieldElement(getElementRep(x)*polyMat_)+constantFactor_;
    }

	//return the i'th coefficient of this polynomial
	FieldElement AffinePolynomial::getCoefficient(const unsigned int i)const {
        if ( i == 0) return constantFactor_;
		
        //if "i" is not a power of 2, return zero
        if ((1<<int(log2(i))) != i) return zero();

        const size_t indexLog = log2(i);
        if (indexLog < coefficients_.size()) return coefficients_[indexLog];

		return zero();
    }
    

    PolynomialDegree AffinePolynomial::getDegree() const{
        const FieldElement ZERO = zero();
        
        if (coefficients_.size() > 0){
            return PolynomialDegree(1LL<<(coefficients_.size()-1));
        }
        
        if (constantFactor_ != ZERO){
            return PolynomialDegree(0);
        }

        //otherwise this is the zero polynomial, and its degree is (-infinity)
        return PolynomialDegree(-1);
    }


	void AffinePolynomial::computeMat() {
		long dim = getFieldDegree();

		if (coefficients_.size() == 0){
			polyMat_= NTL::mat_GF2(NTL::INIT_SIZE, dim, dim);
			return;
		}
		

		FieldElement d = one();
		FieldElement x = mapIntegerToFieldElement(0, dim, 2);
		NTL::vec_GF2 c;
		NTL::mat_GF2 m(NTL::INIT_SIZE, dim, dim);
		//m.SetDims(dim, dim);
		for (long i = 0; i < dim; i++){
			c = getElementRep(evalLinearPart(d,coefficients_));
			m[i] = c;
			d *= x;

		}
		polyMat_ = m;
		}
	

	

} // of namespace Algebra
