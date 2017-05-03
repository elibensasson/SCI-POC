#ifndef  __UnivariatePolynomialGeneral_HPP
#define  __UnivariatePolynomialGeneral_HPP

#include "algebraLib/PolynomialInterface.hpp"

#include <set>
#include <map>
#include <NTL/GF2EX.h>

namespace Algebra{

/**
 * @class UnivariatePolynomialGeneral
 * @brief An implementation of a generic univariate polynomial using the NTL library
 * Describes a polynomial over some
 * field, and not over a specific field.
 * This is due to implementation details,
 * specifically the usage of NTL library that
 * defines a global single context field.
 *
 * @todo Update the interface such that each
 * polynomial can tell what is its context field.
 * This way the usage of NTL and its global field
 * won't be reflected on the interface (2013.07.09)
 */
class UnivariatePolynomialGeneral : public DivisorPolynomial {
public:
	/**
	 * @brief   Default constructor, generates the zero polynomial
	 */
	UnivariatePolynomialGeneral(){};
	
	/**
	 * @brief   semi copy constructor
	 */
	UnivariatePolynomialGeneral(const UnivariatePolynomialInterface& src);
    
	/*sets poly according to coeffs*/
	UnivariatePolynomialGeneral(const std::vector<FieldElement>& coeffs);
    
	/*sets poly according to coeffs, just takes them to itself*/
	UnivariatePolynomialGeneral(std::vector<FieldElement>&& coeffs);
    /**
     * @brief return a clone of the current polynomial
     * @return a unique_ptr of PolynomialInterface,
     * representing a polynomial equivalent to current
     */
    std::unique_ptr<PolynomialInterface> clone()const;
	
    /**
     * @brief   construct using interpolation over a vector subspace 
     * (of the Field, as a vector space over GF2)
     *
     * @param evaluationTabel mapping of the form \f$ x \mapsto y \f$ which the interpulant should agree with
     * @param spaceBasis basis of the vector space \f$ V \f$ over GF2 which is a subspace of the current context field, which is the interpolation domain
     */
    typedef std::map<FieldElement,FieldElement,bool(*)(const FieldElement&,const FieldElement&)> evaluation_t;
    UnivariatePolynomialGeneral(const evaluation_t& evaluationTabel, const elementsSet_t& spaceBasis);
    
    //interpolation not a vector space
    UnivariatePolynomialGeneral(const evaluation_t& evaluationTabel);

    //interpolation over a vector space, represented by a vector of elements
    UnivariatePolynomialGeneral(const std::vector<FieldElement>& evaluation, const std::vector<FieldElement>& orderedBasis, const FieldElement& spaceShift);
    
	/**
	 * @brief   Construct an equivalent polynomial to a given NTL polynomial
	 * @param  p NTL polynomial
	 */
	explicit UnivariatePolynomialGeneral(const NTL::GF2EX& p);

	explicit operator const NTL::GF2EX()const;
	/**
	 * @brief   Constructs a constant polynomial
	 * (a mapping that is independent in the variable)
	 * @param   constant the constant mapping value
	 */
	UnivariatePolynomialGeneral(const FieldElement& constant);

	/**
	 * @brief   Constructs the monic lowest degree 
	 * polynomial with the given roots
	 * @param   roots set of roots
	 */
	UnivariatePolynomialGeneral(const elementsSet_t& roots);

	/**
	 * @brief   Evaluates the polynomial instance using specific element
	 * @param   x the assignment
	 * @return  P(x)
	 */
	FieldElement eval(const FieldElement& x)const;
    
    /**
     * @brief evaluates a polynomial over a givven space (represented by a ordered basis and affine shift)
     * @return the evaluation as a vector, element number \f$ \sum_{i=0}^n b_i \cdot 2^i \f$ represents
     *          the value at the point \f$ shift + \sum_{i=0}^n (basis element)_i \cdot 2^i \f$
     */
    std::vector<FieldElement> eval(const std::vector<FieldElement>& orderedBasis, const FieldElement& shift)const;
	
	/**
	 * @brief   Evaluates the polynomial, as a polynomial over \f$ \mathbb{F}[x] \f$
	 * @param   p input univariate polynomial
	 * @return  The evaluation (composition) \f$this \circ p\f$ such that
	 * \f$(this \circ p)(x) = this(p(x))\f$
	 */
	UnivariatePolynomialInterface* eval(const UnivariatePolynomialInterface& p)const;

	/**
	 * @brief   Changes the value of the coefficient of \f$x^\text{index}\f$
	 * @param   index the power of \f$x\f$ the coefficient multiplies 
	 * @param   cefficient the value of the coefficient
	 */
	void setCoefficient(const unsigned int index,const FieldElement& coefficient);
	
	/**
	 * @brief   returns the value of the coefficient of \f$x^\text{index}\f$
	 * @param   index the power of \f$x\f$ the coefficient multiplies 
	 * @return   cefficient the value of the coefficient
	 */
	 FieldElement getCoefficient(const unsigned int index)const;

     //get all coefficients in vector
     const std::vector<FieldElement> getCoefficients()const;
	
	 /**
	 * @brief   returns the degree of the polynomial
	 * @return   polynomial degree
	 */
	 PolynomialDegree getDegree()const;

	/**
	 * @brief Let this polynomial be \f$a \in \mathbb{F}[x]\f$, and the parameter be \f$b \in \mathbb{F}[x]\f$
	 * than by the devision theorem there exists \f$ q,r \in \mathbb{F}[x] \f$ such that \f$ \deg r < \deg a \f$
	 * and \f$ b = aq + r \f$.
	 * This method return \f$q\f$
	 */
	std::unique_ptr<UnivariatePolynomialInterface> divideByMe(const UnivariatePolynomialInterface& dividend)const; 

	 /**
	  * @brief   Multiply this polynomial by other polynomial
	  * @param   other the other polynomial
	  */
	 void multiply(const UnivariatePolynomialGeneral& other);

	 /**
	  * @brief   Adds this polynomial to other polynomial
	  * @param   other the other polynomial
	  */
	 void add(const UnivariatePolynomialGeneral& other);

	 /** Bunch of operators */
	 friend UnivariatePolynomialGeneral operator+(const UnivariatePolynomialGeneral&, const UnivariatePolynomialGeneral&);
	 friend UnivariatePolynomialGeneral& operator+=(UnivariatePolynomialGeneral&, const UnivariatePolynomialGeneral&);
	 friend UnivariatePolynomialGeneral operator-(const UnivariatePolynomialGeneral&, const UnivariatePolynomialGeneral&);
	 friend UnivariatePolynomialGeneral operator*(const UnivariatePolynomialGeneral&, const UnivariatePolynomialGeneral&);
	 friend UnivariatePolynomialGeneral& operator*=(UnivariatePolynomialGeneral&, const UnivariatePolynomialGeneral&);
	 friend bool operator==(const UnivariatePolynomialGeneral&, const UnivariatePolynomialGeneral&);
	 friend bool operator!=(const UnivariatePolynomialGeneral&, const UnivariatePolynomialGeneral&);
	 friend UnivariatePolynomialGeneral power(const UnivariatePolynomialGeneral&, long exponent);
	 friend std::ostream& operator<<(std::ostream&, const UnivariatePolynomialGeneral&);
private:
	 std::vector<FieldElement> polynomial_; /**< @brief instance state */
	
};

} // namespace Algebra

#endif   // __UnivariatePolynomialGeneral_HPP
