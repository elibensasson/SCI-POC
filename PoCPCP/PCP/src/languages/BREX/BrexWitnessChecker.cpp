#include "BrexWitnessChecker.hpp"
#include "common/Utils/TaskReporting.hpp"

using std::vector;
using Algebra::FieldElement;

namespace PCP_Project {

static bool isPermutation(size_t const numElements, const BREXWitness::permutation_t& mapping){
	//We check the mapping is a bijection (permutation)
	//by checking that it is injective
	
	//inImage[i] will be set to 'true' iff
	//there exists j that is mapped to it
	vector<bool> inImage(numElements);
	
	//initialization
	for(size_t i=0; i< numElements; i++) inImage[i] = false;

	//Check if is an injection
	for(size_t i=0; i< numElements; i++){
		size_t img = mapping.getElementByIndex(i);

		//check the image is inside {0 ... numElements-1}
		if ((img < 0) || (img >= numElements)) return false;

		//check no element was mapped to img before (validate injectivity)
		if (inImage[img]) return false;

		//mark in image
		inImage[img] = true;
	}
	
	return true;
}

static bool isSatisfied(const ConstraintSys& sys, const BREXWitness::color_t& c1, const BREXWitness::color_t& c2){
	vector<FieldElement> assignment(c1.begin(),c1.begin() + (sys.varsAmount()/2));
    assignment.insert(assignment.end(),c2.begin(),c2.begin() + (sys.varsAmount()/2));
	return sys.verify(assignment);
}

bool BREXWitnessChecker::verify(const BREXInstance& instance, const BREXWitness& witness){
	TASK("Executing BREX determenistic witness checker");
    if(!verify_permutation(instance,witness)) return false;
	if(!verify_constraintsAssignment(instance,witness))return false;
    if(!verify_constraintsPermutation(instance,witness))return false;
	if(!verify_boundary(instance,witness))return false;
    return true;
}

bool BREXWitnessChecker::verify_permutation(const BREXInstance& instance, const BREXWitness& witness){
	TASK("verifyin mapping is permutation");
	size_t numElements = instance.domainSize();
	return isPermutation(numElements,witness.permutation());
}

bool BREXWitnessChecker::verify_constraintsAssignment(const BREXInstance& instance, const BREXWitness& witness){
	TASK("verifying assignment constraints");
	size_t domainSize = instance.domainSize();
	const ConstraintSys& constraintsAssignment = instance.constraintsAssignment();
    for (size_t i=0; i < domainSize-1; i++){
        BREXWitness::color_t c1 = witness.get_color(i);
        BREXWitness::color_t c2 = witness.get_color(i+1);
        bool currSatisfied = isSatisfied(constraintsAssignment,c1,c2);
        if (!currSatisfied) {
            return false;
        }
    }

	return true;
}

bool BREXWitnessChecker::verify_constraintsPermutation(const BREXInstance& instance, const BREXWitness& witness){
	TASK("verifying permutation constraints");
	size_t domainSize = instance.domainSize();
	const ConstraintSys& constraintsPermutation = instance.constraintsPermutation();
	const BREXWitness::permutation_t& permutation = witness.permutation();
    for (size_t i=0; i < domainSize-1; i++){
        BREXWitness::color_t c1 = witness.get_color(i);
        size_t perm_next = permutation.getElementByIndex(i);
		BREXWitness::color_t c2 = witness.get_color(perm_next);
        bool currSatisfied = isSatisfied(constraintsPermutation,c1,c2);
        if (!currSatisfied) {
            return false;
        }
    }

	return true;
}

bool BREXWitnessChecker::verify_boundary(const BREXInstance& instance, const BREXWitness& witness){
	TASK("verifying boundary constraints");

	for(const auto& p : instance.boundaryConstraints()){
		const BREXInstance::point_t inputLocation = p.first;
		const size_t vector_index = inputLocation.first;
		const size_t element_index = inputLocation.second;

		const  BREXWitness::color_t color = witness.get_color(vector_index);
		assert(color.size() == instance.vectorsLen());
		if(!(p.second == color[element_index])) {
			return false;
		}
	}
	return true;
}

} // namespace PCP_Project
