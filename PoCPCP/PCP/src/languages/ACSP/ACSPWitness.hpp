/**
 *       @file  ACSPWitness.hpp
 *      @brief  Header file for ACSP witness
 *
 *     @author  Michael Riabzev, RiabzevMichael@gmail.com
 * =====================================================================================
 */
#pragma once // This is not needed, here just to get rid of an annoying VS2010 warning.
#ifndef  __ACSP_WITNESS_HPP
#define  __ACSP_WITNESS_HPP

#include <algebraLib/PolynomialInterface.hpp>
#include "common/CXX11_macros.hpp"

#include <memory>

namespace PCP_Project {

/**
 * @class ACSPWitness
 * @brief class for ACSP witness
 *
 * This class describes a witness for ACSPFullInstance.
 * A witness is a polynomial \f$ A \in \mathbb{F}[x] \f$,
 * and named Assignment Polynomial.
 * Such a polynomial shows that\n
 * \f$(\mathbb{F},H,\vec{N},P,witnessDegreeBound,B)\f$ is a satisfiable ACSPInstance
 * if and only if
 * \f{eqnarray*}{
 * \forall z \in H: (P \circ (x \vert A \circ \vec{N}))(z) = 0 \\
 * \deg A <= witnessDegreeBound \\
 * \lambda \cdot \deg(P \circ (x \vert A \circ \vec{N})) \le |\mathbb{F}| \\
 * \forall (x,y) \in B : A(x)=y \\
 * \f}
 *
 * In the code we name the assignment polynomial \f$A\f$
 * simply 'assignmentPoly'
 *
 * Methods:
 * Witness class contains only getters, 
 * an empty constructor and a destructor.
 * The only possible way to change its data members
 * is using the reduction class from EGCP or inside a UTEST.
 */
class ACSPWitness {
public:
	typedef ::Algebra::UnivariatePolynomialInterface polynomial;

	ACSPWitness(std::unique_ptr<const polynomial>&& assignmentPoly):assignmentPoly_(std::move(assignmentPoly)){};
	
	/**
	 * @brief   Move constructor
	 * should be default, but not supported
	 * by MS VS2010, so just implements 'default'
	 */
	ACSPWitness(ACSPWitness&& src):assignmentPoly_(std::move(src.assignmentPoly_)){};// = default;

	inline const polynomial& assignmentPoly() const {
		return *(assignmentPoly_.get());
	}

private:

	std::unique_ptr<const polynomial> assignmentPoly_;

	/**
	 * @brief  Copy constructor
	 * should be deleted, but not supported
	 * by MS VS2010, so it is private and throws exception 
	 */
	ACSPWitness(const ACSPWitness& src) _COMMON_CXX11_DELETED;
};

} //namespace PCP_Project

#endif //__ACSP_WITNESS_HPP
