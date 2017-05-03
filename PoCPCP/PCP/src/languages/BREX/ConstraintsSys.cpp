#include "ConstraintsSys.hpp"

namespace PCP_Project{

bool ConstraintSys::varUsed(const size_t varId) const{
    for(const auto& poly : constraints()){
        if(poly->isEffectiveInput(varId)){
            return true;
        }
    }
    return false;
}


bool ConstraintSys::verify(const std::vector<Algebra::FieldElement>& assignment) const{
	for (const polyPtr_t& c : constraints()){
		if (c->eval(assignment) != Algebra::zero()){
			return false;
		}
	}
	return true;
}

Algebra::PolynomialDegree ConstraintSys::getMaximalDegree()const{
    using Algebra::PolynomialDegree;
    using std::vector;
    using std::max;

    const vector<PolynomialDegree> inputDegrees(varsAmount(),PolynomialDegree(1));
    PolynomialDegree res = PolynomialDegree::getZeroPolyDegree();
    
    for(const auto& poly : constraints()){
        res = max(res,poly->getDegreeBound(inputDegrees));
    }

    return res;
}

} // namespace PCP_Project
