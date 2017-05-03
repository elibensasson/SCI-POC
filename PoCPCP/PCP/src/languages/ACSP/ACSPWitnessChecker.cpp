#include "ACSPWitnessChecker.hpp"
#include "common/Utils/TaskReporting.hpp"
#include "common/Algebra/details/Polynomials.hpp"
#include <memory>
#include <vector>

using std::unique_ptr;
using Algebra::FieldElement;

namespace PCP_Project {


/**
 * @class notRootPredicate
 * @brief A predicate class
 * that checks if a given field element is
 * not a root of some internally kept polynomial
 */
class notRootPredicate : public ::Algebra::FieldElementPredicate {
public:
	/**
	 * @brief   The constructor
	 * @param   instnace that contains a constraint polynomial
	 * \f$P:\mathbb{F}^n \to \mathbb{F}\f$ and a neighbors polynomials
	 * vector \f$\vec{N}:\mathbb{F}^n \to  \mathbb{F}^n\f$
	 * @param   witness that contains an assignment polynomial \f$A:\mathbb{F} \to \mathbb{F}\f$
	 */
	notRootPredicate(const ACSPInstance& instance, const ACSPWitness& witness): instance_(instance), witness_(witness){};

	/**
	 * @brief   checks that the field element is not a root of some polynomial
	 * @param   x the field element to check
	 * @return  True iff \f$ 0 \neq P(x,A(N_1(x),A(N_2(x),\dots,A(N_n(x)) \f$
	 */
	bool test(const ::Algebra::FieldElement& x)const {
        vector<FieldElement> assignment;
        assignment.push_back(x);
        for(const auto& n : instance_.neighborPolys()){
            assignment.push_back(witness_.assignmentPoly().eval(n->eval(x)));
        }
        return !instance_.constraintPoly().isRoot(assignment);
	}

private:
    const ACSPInstance& instance_;
    const ACSPWitness& witness_;
};



bool ACSPWitnessChecker::verify(const ACSPInstance& instance, const ACSPWitness& witness){
	
	TASK("Executes deterministic  ACSP witness checker");
	
	if (!verify_vanishing(instance,witness)) return false;

	if (!verify_witness_degree(instance,witness)) return false;
		
	if (!verify_boundary(instance,witness)) return false;

    return true;
}

bool ACSPWitnessChecker::verify_vanishing(const ACSPInstance& instance, const ACSPWitness& witness){
	TASK("Tests vanishing constraint");

	const ACSPInstance::set& vanisinigSet = instance.vanishingSet();
	notRootPredicate* predicate = new notRootPredicate(instance,witness);
	
	std::unique_ptr<const ::Algebra::FieldElementPredicate> pred_ptr(predicate);
	return !(vanisinigSet.exist(pred_ptr));

}

bool ACSPWitnessChecker::verify_witness_degree(const ACSPInstance& instance, const ACSPWitness& witness){
	TASK("Tests witness degree constraint");
    return !(instance.witnessDegreeBound() < witness.assignmentPoly().getDegree());
}

bool ACSPWitnessChecker::verify_boundary(const ACSPInstance& instance, const ACSPWitness& witness){
	TASK("Tests boundary constraints");
	const ACSPWitness::polynomial& assignment = witness.assignmentPoly();

	for(const auto& p : instance.boundaryConstraints()){
		if(assignment.eval(p.first) != p.second) {
			return false;
		}
	}
	return true;
}

} // namespace PCP_Project
