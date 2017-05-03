#ifndef  __ACSPWitnessChecker_HPP
#define  __ACSPWitnessChecker_HPP

#include "ACSPInstance.hpp"
#include "ACSPWitness.hpp"

namespace PCP_Project {

/**
 * @class ACSPWitnessChecker
 * @brief class for ACSP (instance, witness) pairs checking
 */
class ACSPWitnessChecker {
public:
	/**
	 * @brief   verifies all conditions
	 * @param   instance a tuple \f$(\mathbb{F},H,\vec{N},P,WitnessDegreeBound,B)\f$
	 * @param   witness a polynomial \f$ A \in \mathbb{F}[x]\f$
	 * @return  return true if and only if
	 * \f{eqnarray*}{
	 * \forall z \in H: (P \circ (x \vert A \circ \vec{N}))(z) = 0 \\
	 * \deg A <= witnessDegreeBound \\
	 * \lambda \cdot \deg(P \circ (x \vert A \circ \vec{N})) \le |\mathbb{F}| \\
	 * \forall i \in [n] : x_i = A(I(i)) \\
	 * \f}
	 */
	static bool verify(const ACSPInstance& instance, const ACSPWitness& witness);
	
	/**
	 * @brief   verifies the vanishing subset condition
	 * @param   instance a tuple \f$(\mathbb{F},H,\vec{N},P,WitnessDegreeBound,B)\f$
	 * @param   witness a polynomial \f$ A \in \mathbb{F}[x]\f$
	 * @return  return true if and only if
	 * \f$
	 * \forall z \in H: (P \circ (x \vert A \circ \vec{N}))(z) = 0
	 * \f$
	 */
	static bool verify_vanishing(const ACSPInstance& instance, const ACSPWitness& witness);
	
	/**
	 * @brief   verifies the witness degree condition
	 * @param   instance a tuple \f$(\mathbb{F},H,\vec{N},P,WitnessDegreeBound,B)\f$
	 * @param   witness a polynomial \f$ A \in \mathbb{F}[x]\f$
	 * @return  return true if and only if
	 * \f$
	 * \deg A <= witnessDegreeBound
	 * \f$
	 */
	static bool verify_witness_degree(const ACSPInstance& instance, const ACSPWitness& witness);
	
	/**
	 * @brief   verifies the boundary constraints condition
	 * @param   instance a tuple \f$(\mathbb{F},H,\vec{N},P,WitnessDegreeBound,B)\f$
	 * @param   witness a polynomial \f$ A \in \mathbb{F}[x]\f$
	 * @return  return true if and only if
	 * \f$
	 * \forall (x,y) \in B : A(x)=y$
	 * \f$
	 */
	static bool verify_boundary(const ACSPInstance& instance, const ACSPWitness& witness);
private:
};

} // namespace PCP_Project
#endif   // __ACSPWitnessChecker_HPP
