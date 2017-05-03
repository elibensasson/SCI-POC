/**************************************** BREXtoACSP.cpp *****************************************/
/**
 * @file.
 *
 * Implements the BREXtoACSP reduction class.
 */
  /***********************************************************************************************/

#include "BREXtoACSP.hpp"
#include "common/Infrastructure/Infrastructure.hpp"
#include "common/Algebra/LinearSpace.hpp"
#include "common/Utils/TaskReporting.hpp"
#include "Details/witnessReduction.hpp"
#include "Details/constraints.hpp"
#include "Details/common.hpp"
#include "Details/instanceMappings.hpp"
#include "Details/neighborsConstructor.hpp"
#include "Details/spaces.hpp"
#include "Details/boundaryConstraints.hpp"
#include "Details/ACSPSummandsPolynomial.hpp"
#include <algebraLib/UnivariatePolynomialGeneral.hpp>
#include <PCP/PCP_common.hpp>
#include "common/Algebra/details/FiniteFields.hpp" //for "basisShiftPair"

#include <queue>

namespace PCP_Project {

using Algebra::FiniteField;
using Algebra::LinearSpace;
using Algebra::UnivariatePolynomialInterface;
using PCP_Project::BREXtoACSP::instanceMappings;
using PCP_Project::BREXtoACSP::common;
using PCP_Project::BREXtoACSP::AcspNeighbors;
using PCP_Project::BREXtoACSP::CS_testLocations;
using PCP_Project::BREXtoACSP::spaces;
using PCP_Project::BREXtoACSP::constraints;
using PCP_Project::BREXtoACSP::reduceBoundaryConstraints;

using std::unique_ptr;
using std::move;

unique_ptr<ACSPWitness> CBREXtoACSP::reduceWitness( const BREXInstance& instance, const BREXWitness& witness){

    return BREXtoACSP::witnessReduction::reduceWitness(instance,witness);
}


/**********************************************
 *             Instance Reduction
 *********************************************/

namespace{
    //get vanishing space properties
    Algebra::details::basisShiftPair getVanishingSpace(const spaces& spacesGenerator){
        Algebra::details::basisShiftPair result;

        const auto basis = spacesGenerator.getVanishingSpaceBasis();
        for(const auto& e : basis){
            //construct the standard basis
            result.basis.addElement((e));
        }

        //return the result
        return result;
    }
	

std::unique_ptr<UnivariatePolynomialInterface> compositionAlg(const ACSPInstance& instance, const UnivariatePolynomialInterface& witness){
    using std::vector;
    using std::min;
    using std::max;
    using std::unique_ptr;
    using Algebra::FieldElement;
    using Algebra::PolynomialDegree;
    using Algebra::zero;
    using Algebra::UnivariatePolynomialGeneral;
    using Algebra::mapIntegerToFieldElement;
    using Infrastructure::Log2;
    using BREXtoACSP::ACSPSummandsPolynomial;

    //the constraints poly
    const ACSPSummandsPolynomial& constraintsPoly = *(dynamic_cast<const ACSPSummandsPolynomial*>(&instance.constraintPoly()));
    
    //
    //get the composition degree bound
    //
    vector<PolynomialDegree> constraintsInputDegrees;

    // first input is "x" which has degree 1
    constraintsInputDegrees.push_back(PolynomialDegree(1));

    // rest are composition of neighbor with witness
    const auto witnessDegree = witness.getDegree();
    for (const auto& n : instance.neighborPolys()){
        constraintsInputDegrees.push_back(n->getDegreeBound(witnessDegree));
    }

    // get the composition degree bound
    const PolynomialDegree degBound = constraintsPoly.getDegreeBound_DividedByZH(constraintsInputDegrees);

    //
    //define interpolation space
    //
    const size_t degForPoly = degBound.isInteger()? ceil(log2(1+PolynomialDegree::integral_t(degBound))) : 0;
    const size_t interpolationBasisSize = min(instance.contextField().degree(), degForPoly);
    const size_t witnessEvalBasisSize = min(instance.contextField().degree(), degForPoly+2);
    const auto interpolationBasis = Algebra::details::buildStandardBasis(interpolationBasisSize);
    const auto witnessEvalBasis = Algebra::details::buildStandardBasis(witnessEvalBasisSize);
    const size_t vanishingSpaceDeg = Log2(instance.vanishingSet().size());
    const FieldElement testsShift = mapIntegerToFieldElement(interpolationBasis.asVector().size(),1,1);

    vector<FieldElement> evaluation(Infrastructure::POW2(interpolationBasis.asVector().size()));
    
    {
    //
    //Evaluate witness on relevant spaces for faster composition
    //
    vector<FieldElement> witnessEvaluation;
    {
	    TASK("Evaluate the witness polynomial over 'big enough' spaces");
		witnessEvaluation = witness.eval(witnessEvalBasis.asVector(),zero());
    }
	//question to Michael: This is not part of BREX->ACSP reduction
    //construct evaluation
    {
        TASK("Evaluate the composition polynomial divided by Z_H over a space of degree " + std::to_string(interpolationBasis.asVector().size()));
        const size_t LOG_BLOCK_SIZE = 2;
        const size_t BLOCK_SIZE = 1<<LOG_BLOCK_SIZE;
        const size_t NUM_EVALUATION_BLOCKS = evaluation.size() >> LOG_BLOCK_SIZE;
		
        #pragma omp parallel for schedule(guided) num_threads(SCIPR_NUM_THREADS)
        for(int i=0; i< NUM_EVALUATION_BLOCKS; i++){
        
            const size_t numThreads = omp_get_num_threads();
            const size_t currThreadId = omp_get_thread_num();
        
#ifdef PRINT_PROGRESS
            if (i % 4000 == 0){
                std::cout << i << "/" << NUM_EVALUATION_BLOCKS << ",";
                std::fflush(stdout);
            }
#endif //PRINT_PROGRESS
            
            std::queue<size_t> indexesToEval;
            vector<vector<FieldElement>> assignment;
            for(size_t j=0; j<BLOCK_SIZE; j++){
                //iteration constants
                const size_t curr_index = (i<<LOG_BLOCK_SIZE)+j;
                const FieldElement x = getSpaceElementByIndex(interpolationBasis.asVector(),testsShift,curr_index);
                
                indexesToEval.push(curr_index);

                //construct the assignment for the constraints poly
                vector<FieldElement> curr_assignment(1+instance.neighborPolys().size());
                curr_assignment[0] = x;
                for(size_t n=0; n < instance.neighborPolys().size(); n++){
                    const FieldElement neighborRes = instance.neighborPolys()[n]->eval(x);
                    const size_t neighborRes_indx = mapFieldElementToInteger(0,witnessEvalBasisSize,neighborRes);
                    curr_assignment[1+n] = witnessEvaluation[neighborRes_indx];
                }
                assignment.push_back(curr_assignment);
            }
            //evaluate and return
            if(!assignment.empty()){
#ifdef ZERO_COMPOSITION
				const vector<FieldElement> res(assignment.size(), zero());
#else
                const auto res = constraintsPoly.evalDividedByZH(assignment);

#endif
                for(const auto& r : res){
                    evaluation[indexesToEval.front()] = r;
                    indexesToEval.pop();
                }
            }
        }
    }
    }

    {
    TASK("Interpolating composition evaluation");
    Algebra::IFFT_inplace(&(evaluation[0]),interpolationBasis.asVector(),testsShift);
    }
    
    //
    //return result
    //
    return unique_ptr<UnivariatePolynomialInterface>(new UnivariatePolynomialGeneral(std::move(evaluation)));
}

}
/*
 * Construct an ACSP partial instance from a BREX full instance (top-level arithmetization function for partial instances)
 */
unique_ptr<ACSPInstance> CBREXtoACSP::reduceInstance(const BREXInstance&  instance){
    // get common deffintion
    const common commonDef(instance);
    
    // get the context field and set it as current context
    const FiniteField contextField = commonDef.contextField();
    contextField.setContext();

    // fetch more parameters
    const instanceMappings instanceMapping(commonDef);
    const spaces spacesGenerator(commonDef);
    const CS_testLocations testLocations(commonDef);
    const AcspNeighbors neighborsWithInterpretation(instance, commonDef, instanceMapping,testLocations);
    
    // construct the vanishing set
    unique_ptr<ACSPInstance::set> vanishingSet(new LinearSpace(getVanishingSpace(spacesGenerator)));
   
    // construct the neighbors vector
    ACSPInstance::polynomialsVec neighbors = neighborsWithInterpretation.getNeighborPolynomials();
    
    //construct the constraints polynomial
    unique_ptr<const ACSPInstance::polynomial> constraintsPoly(constraints::getConstraintsPoly(instance,neighborsWithInterpretation,commonDef,instanceMapping, spacesGenerator, testLocations));
    
    //construct the boundary constraints
    ACSPInstance::boundaryConstraints_t boundaryConstraints = reduceBoundaryConstraints(instance.boundaryConstraints(),commonDef);

    //construct and return
    unique_ptr<ACSPInstance> res(new ACSPInstance(contextField, move(vanishingSet), move(neighbors), move(constraintsPoly), commonDef.witnessDegreeBound(), boundaryConstraints, compositionAlg));

    return res;
}

} // namespace PCP_Project

