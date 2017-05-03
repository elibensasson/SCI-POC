#include "FieldElement.hpp"

#include <NTL/GF2X.h>
#include <set>

#ifndef __FINITE_FIELD_HPP
#define __FINITE_FIELD_HPP

namespace Algebra {

class FiniteField {
public:
	/**
	 * @brief   Default constructor
	 * constructs the prime field GF(2)
	 */
	FiniteField();
	
	/**
	 * @brief   copy constructor
	 */
	FiniteField(const FiniteField&);

    /**
     * @brief generate using NTL irreducible
     */
    FiniteField(const NTL::GF2XModulus& irr):modulus_(irr){};
	
	/**
	 * @brief   Constructor, constructs a field of some characteristic
	 * and extension degree
	 * @param   characteristic the characteristic of the constructed field
	 * @param   extensionDegree the extension degree of the field over
	 * the prime field
	 * @param useMatan a flag that currently enforces field degree 64 to use Matan's FFT (default set to true)
	 */
	FiniteField(const size_t extensionDegree);

	/**
	 * @brief   Returns the extension degree of the field (over the prime field)
	 * @return  extension degree
	 */
	size_t degree()const;

	/**
	 * @brief   returns a basis for this field as a vector space over GF2 
	 * @return  the basis as a vector of FieldElements
	 * @note    I'm not familiar with our FFT which this basis is used for
	 * but maybe the order does not matter and it could be a set of
	 * FieldElements instead. Mailed Arnon about it. 
	 */
	elementsSet_t getBasis()const;

	/**
	 * return the field extension modulus
	 */
	const NTL::GF2XModulus& getModulus()const{return modulus_;}	
	
	/**
	 * Sets this field as the context field for all operations
	 */
	void setContext()const;

    /**
     * Must be called in order to release the context field
     * Best do it before program termination,
     * although leaks are not harmful on that time,
     * but it is nicer to see valgrind is happy
     */
    static void releaseContextField();

private:
	NTL::GF2XModulus modulus_;
};

} //namespace Algebra

#endif // __FINITE_FIELD_HPP
