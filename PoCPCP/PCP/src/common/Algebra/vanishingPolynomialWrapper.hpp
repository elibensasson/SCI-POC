#ifndef ALGEBRA_VANISHING_POLY_WRAPPER_HPP__
#define ALGEBRA_VANISHING_POLY_WRAPPER_HPP__

#include <algebraLib/UnivariatePolynomialGeneral.hpp>
#include "details/Polynomials.hpp"

namespace Algebra{
/**
 * This should be a temporary wrapper that implements univariate polynomial interface
 * Should be used only untill a good solution is written, that
 * includes 'affine polynomials' (linearized + constant) in the new
 * algebra library
 */
class vanishingPolynomialWrapper : public DivisorPolynomial{
public:
	details::VanishingPolynomial model;
	
	vanishingPolynomialWrapper(const details::Basis& basis, const FieldElement& elementShift) : model(basis,elementShift){};
	vanishingPolynomialWrapper(const details::Basis& basis) : model(basis){};
	
	FieldElement eval(const FieldElement& x)const {
		return (FieldElement)model.evalLinearizedPoly((FieldElement)x);
	}
	FieldElement getCoefficient(const unsigned int index)const {
		if (index == 0) return FieldElement(model.getConstantFactor());
		if (!Infrastructure::IsPower2(index)) return zero();
		return FieldElement(model.getCoeff(Infrastructure::Log2(index)));
	}
	PolynomialDegree getDegree()const {return PolynomialDegree(model.getDegree());}
	
	std::unique_ptr<UnivariatePolynomialInterface> divideByMe(const UnivariatePolynomialInterface& dividend)const {
		//generate an instance of the old polynomials
		//not very efficient, but saves programming time,
		//later division can be done directly
		details::UnivariatePolynomial tmp;
		//ARIEL:the next three lines added because of the technicality that the 0's polynomial degree is negative
		size_t degBound;
		if (dividend.getDegree().isInteger())
			degBound = PolynomialDegree::integral_t(dividend.getDegree());
		else
			degBound = 0;
		for(size_t i=0; i<= degBound ; i++){
			tmp.setCoeff(i,(FieldElement)(dividend.getCoefficient(i)));
		}
		if (model.getDegree() > PolynomialDegree::integral_t(dividend.getDegree()))
			_COMMON_FATAL("trying to divide a poly p by a polynmial of larger degree than p in divideByMe method");
		//calculate the division
		details::UnivariatePolynomial q, r;
		tmp.divideBySparse(model, q, r);

		//return result
		return std::unique_ptr<UnivariatePolynomialGeneral>(new UnivariatePolynomialGeneral((NTL::GF2EX)q));	
	}
    
    /**
     * @brief return a clone of the current polynomial
     * @return a unique_ptr of PolynomialInterface,
     * representing a polynomial equivalent to current
     */
    std::unique_ptr<PolynomialInterface> clone()const{
        using std::unique_ptr;
        return unique_ptr<PolynomialInterface>(new vanishingPolynomialWrapper(model));
    }

private:
	vanishingPolynomialWrapper(const details::VanishingPolynomial& otherModel) : model(otherModel){};
};
} // namespace Algebra

#endif// ALGEBRA_VANISHING_POLY_WRAPPER_HPP__
