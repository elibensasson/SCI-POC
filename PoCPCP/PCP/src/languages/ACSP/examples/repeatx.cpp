#include "common/Algebra/AlgebraCommon.hpp"
#include "languages/ACSP/ACSPWitnessChecker.hpp"
#include "reductions/BREXtoACSP/BREXtoACSP.hpp"
#include "common/Algebra/details/FiniteFields.hpp"
#include "common/Algebra/details/Polynomials.hpp"
#include "common/Algebra/LinearSpace.hpp"
#include <algebraLib/UnivariatePolynomialGeneral.hpp>
#include "common/Infrastructure/Infrastructure.hpp"
#include "common/Algebra/MultiVarPoly.hpp"
#include "PCP/Tests/ACSP_PCP/PCPWorkflowTests.h"
#include "common/Algebra/LightUniPolyEval.hpp"
#include "common/Utils/TaskReporting.hpp"
//#include "common/lightCircLib/lightCircPoly.hpp"
//#include "common/lightCircLib/lightCircGate.hpp"
using namespace Algebra::details;
using namespace Algebra;
using namespace Infrastructure;
using namespace std;
using namespace PCP_Project;
//using namespace PCP_Project::LightCircLib;
///**
//* @class randomSequence
//* @brief A mapping of \f$\mathbb{N}\f$ into the field
//* such that the first \f$len\f$ integers are mapped into random field elements
//* and the rest are mapped to zeros
//*/
//class randSequence : public PCP_Project::Sequence<FieldElement> {
//public:
//	/**
//	* @brief   The constructor
//	* @param   n boundary of indexes that can be mapped to non-zero elements (namely \f$len\f$)
//	*/
//	randSequence(Sequence<FieldElement>::index_t n) :order(n){
//		for (unsigned int i = 0; i<order.size(); i++){
//			order[i] = generateRandom<FieldElement>();
//		}
//	}
//

#ifndef _MSC_VER//gcc doesn't recognize make_unique do need to define. MSVS does.
namespace std {
	template<typename T, typename... Args>
	std::unique_ptr<T> make_unique(Args&&... args)
	{
		return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
	}
}
#endif
    
std::unique_ptr<UnivariatePolynomialInterface> compositionAlg(const ACSPInstance& instance, const UnivariatePolynomialInterface& witness){
    using std::vector;
    using std::min;
    using std::unique_ptr;
    using Algebra::FieldElement;
    using Algebra::PolynomialDegree;
    using Algebra::zero;
    using Algebra::UnivariatePolynomialGeneral;

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
    const PolynomialDegree degBound = max(instance.constraintPoly().getDegreeBound(constraintsInputDegrees),witnessDegree);

    //
    //define interpolation space
    //
    const size_t degForPoly = degBound.isInteger()? ceil(log2(1+PolynomialDegree::integral_t(degBound))) : 0;
    const auto interpolationBasis = Algebra::details::buildStandardBasis(min(instance.contextField().degree(), degForPoly));

    //
    //Evaluate witness on relevant spaces for faster composition
    //
    vector<LightUniPolyEval> ANEval(1+instance.neighborPolys().size());
    {
	    TASK("Evaluate the witness compositions with neighbors, over 'big enough' spaces");
		#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
        for (int n=0; n<instance.neighborPolys().size();  n++){
			{
                std::cout<<n<<"/"<<instance.neighborPolys().size()<< ",";
                std::fflush(stdout);
            }
            const UnivariatePolynomialInterface& N = *(instance.neighborPolys()[n]);
            //for each neighbor N constructing the basis N(e1),..,N(e_d) , where d is the size of the interpolation basis
			Algebra::details::Basis currentBasis;
			//fullEval_ flag should be true only when all Neighbors have degree exactly one, as otherwise we are trying to apply FFT on something that may not be a basis of a subspace
			_COMMON_ASSERT(PolynomialDegree::integral_t(N.getDegree()) == 1, "flag should be true only when all Neighbors have degree exactly one");
			for (auto i = 0; i < interpolationBasis.asVector().size(); i++){
                const FieldElement e = interpolationBasis.asVector()[i];
				currentBasis.addElement(N.eval(e));//mess is because annoying need to go back and forth between FieldElement and its wrapper FieldElement
			}
			//evaluation of A(N_i)
			ANEval[n] = LightUniPolyEval(witness, currentBasis,zero());
			}
    }

    //construct evaluation
    vector<FieldElement> evaluation(Infrastructure::POW2(interpolationBasis.asVector().size()));
    {
        TASK("Evaluate the composition polynomial over 'big enough' spaces");
		#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
        for(int i=0; i< evaluation.size(); i++){
            if (i % 1000 == 0){
                std::cout << i << "/" << evaluation.size() << ",";
                std::fflush(stdout);
            }
            const FieldElement x = getSpaceElementByIndex(interpolationBasis.asVector(),zero(),i);
            
            //if x is in the vanishing space we assume the result is zero
            if(instance.vanishingSet().contains(x)){
                evaluation[i] = zero();
                continue;
            }
            
            //construct the assignment for the constraints poly
            vector<FieldElement> assignment(1+instance.neighborPolys().size());
            assignment[0] = x;
            for(size_t n=0; n < instance.neighborPolys().size(); n++){
                assignment[1+n] = ANEval[n].queryAtPoint(i);
            }

            //evaluate and return
            evaluation[i] = instance.constraintPoly().eval(assignment);
        }
    }

    //
    //return result
    //
    TASK("Interpolating composition evaluation");
    UnivariatePolynomialGeneral compositionPoly(evaluation,interpolationBasis.asVector(),zero());
    
    //build denominator and divide
    return instance.vanishingSet().vanishingPoly()->divideByMe(compositionPoly);
}

//#define SANITY_CHECK //mess up witness to see verifer won't take crap

static pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness>
generate_repeatx() {

	const size_t extensionDegree = 64;
	const size_t vanishingSpaceDim = 4;

	/**
	* define context field
	*/
	FiniteField contextField( extensionDegree);
	contextField.setContext();
	elementsSet_t fieldBasis = contextField.getBasis();

	/**
	* Construct vanishing set - affine
	*/
	Algebra::details::basisShiftPair BasisH;
	FieldElement ONE = one();

	FieldElement x = FieldElement(NTL::to_GF2E(NTL::GF2X(1, 1)));

//complicated way of getting generator g of F_{16}(doing it in this step by step way, cause NTL expectes long as power parameter)
//First getting generator z of F_{2^16} , z= x^{1+2^(16)+2^(32)+2^(48)}
	FieldElement z1 = x;
	FieldElement z2 = FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(z1), POW2(16)));
	FieldElement z3 = FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(z2), POW2(16)));
	FieldElement z4 = FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(z3), POW2(16)));
	FieldElement z = z1*z2*z3*z4;
//Now getting g = z^(1+2^4+2^8+2^12)
	FieldElement g = FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(z), 1 + POW2(4) + POW2(8) + POW2(12)));


	


	BasisH.basis.addElement(ONE);
	BasisH.basis.addElement(g);
	FieldElement t = g*g;
	BasisH.basis.addElement(t);
	BasisH.basis.addElement(g*t);
	/*
	NTL::GF2EInfo->power(t,2);
	BasisH.basis.addElement(t);
	t = g;
	NTL::GF2EInfo->power(t,3);
	BasisH.basis.addElement(t);
	*/

	unique_ptr<LinearSpace> vanishingSet(new LinearSpace(BasisH));

	/**
	* Construct neighbours
	*/
	ACSPInstance::polynomialsVec neighborPolys;

	
	//construct the vector of neighbor polynomials N= (X,g*X)

	unique_ptr<UnivariatePolynomialGeneral> tPoly(
		std::make_unique<UnivariatePolynomialGeneral>(
		UnivariatePolynomialGeneral()));
	tPoly->setCoefficient(1, ONE);
	neighborPolys.push_back(move(tPoly));

	tPoly =
		std::make_unique<UnivariatePolynomialGeneral>(
		UnivariatePolynomialGeneral());
	tPoly->setCoefficient(1, (g));
	neighborPolys.push_back(move(tPoly));
	/**
	* Construct constraints polynomial:
	P(z,N1z,N2z)=z*(g^15-z)*(N2z-x*N1z)
	=z*(1-z)*(N2z-x*N1z)
	=A(1+A)(C+xB)=A*(C+xB+AC+xAB)
	=AC+xAB+A^2C+xA^2B
	*/
unique_ptr<MultiVarPoly> ConstraintsPolynomial(
		make_unique<MultiVarPoly>());
	MultiVarMonomial M;
	M.coeff = (x);
	M.vars = { 0,0,1 };
	ConstraintsPolynomial->AddMonomial(M); //xA^2B
	ConstraintsPolynomial->AddMonomial({0,0,2 }); //A^2C
	M.vars = { { 0, 1 } };
	ConstraintsPolynomial->AddMonomial(M); //xAB
	ConstraintsPolynomial->AddMonomial({ 0,2 } ); //AC

	/**
	* Construct witness
	*/
	UniPolynomialEvaluation eval;
	vector<FieldElement> points;

	//   evaluation_t tMap(Algebra::compareFieldElements);
	t = g;
	FieldElement nPt =  x*x*x + x + ONE; //arbitrary initial value
	points.push_back(t);
#ifdef SANITY_CHECK
		eval.addPoint(t, nPt+ONE);
#else 
	eval.addPoint(t, nPt );
#endif
	for (int i = 2; i<POW2(4); ++i) {
		nPt = x * nPt;
		t = g * t;
		eval.addPoint(t, nPt);
		points.push_back(t);

	}

	//adding root at zero 
	/*FieldElement zero = zero();
	eval.addPoint(zero, ONE);
	points.push_back(zero);*/
	UnivariatePolynomial P;
	P.interpolation(eval, points);
	unique_ptr<UnivariatePolynomialGeneral> P2(new UnivariatePolynomialGeneral());
	for (unsigned long int i = 0; i<=P.getDegree(); ++i) {
		FieldElement coeff = (P.getCoeff(i));
		P2->setCoefficient(i, coeff);
	}

	ACSPWitness witness(move(P2));

	/** Construct the instance data */
	ACSPInstance instance(contextField, move(vanishingSet), move(neighborPolys), move(ConstraintsPolynomial),witness.assignmentPoly().getDegree(),ACSPInstance::boundaryConstraints_t(),compositionAlg);

	/** Return ACSP pair */
	return move(pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness>(move(instance), move(witness)));
}

static pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness>
generate_repeatx(int subfieldDim) {

	const size_t extensionDegree = 64;
	
	/**
	* define context field
	*/
	FiniteField contextField( extensionDegree);
	contextField.setContext();
	elementsSet_t fieldBasis = contextField.getBasis();

	/**
	* Construct vanishing set - affine
	*/
	Algebra::details::basisShiftPair BasisH;
	FieldElement ONE = one();

	FieldElement x = FieldElement(NTL::to_GF2E(NTL::GF2X(1, 1)));

	//complicated way of getting generator g of F_{2^(subfieldDim)}(doing it in this step by step way, cause NTL expectes long as power parameter)
	//First getting generator z of F_{2^16} , z= x^{1+2^(16)+2^(32)+2^(48)}
	FieldElement z1 = x;
	FieldElement z2 = FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(z1), POW2(16)));
	FieldElement z3 = FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(z2), POW2(16)));
	FieldElement z4 = FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(z3), POW2(16)));
	FieldElement z = z1*z2*z3*z4;

	//Now getting g = z^(1+2^(subfieldDim)+2^(2*subfieldDim)..)
	FieldElement g;
	if (subfieldDim >8)
		g = z;
	else {
		if (subfieldDim == 8)
			g = FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(z), 1 + POW2(8)));
		else
			if (subfieldDim == 4){
				g = FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(z), 1 + POW2(4) + POW2(8) + POW2(12)));
			}
			else
				_COMMON_ASSERT(false, "repeatx test works with parameter 4,8 or 16");


	}

	FieldElement current = ONE;
	for (int i = 0; i < subfieldDim; i++){
		BasisH.basis.addElement(current);
		current *= g;
	}
	unique_ptr<LinearSpace> vanishingSet(new LinearSpace(BasisH));

	/**
	* Construct neighbours
	*/
	ACSPInstance::polynomialsVec neighborPolys;


	//construct the vector of neighbor polynomials N= (X,g*X)

	unique_ptr<UnivariatePolynomialGeneral> tPoly(
		std::make_unique<UnivariatePolynomialGeneral>(
		UnivariatePolynomialGeneral()));
	tPoly->setCoefficient(1, ONE);
	neighborPolys.push_back(move(tPoly));

	tPoly =
		std::make_unique<UnivariatePolynomialGeneral>(
		UnivariatePolynomialGeneral());
	tPoly->setCoefficient(1, (g));
	neighborPolys.push_back(move(tPoly));
	
	FieldElement s= x*x*x + x + ONE; //arbitrary initial value
	/**
	* Construct constraints polynomial:
	P(z,N1z,N2z)=z*(N2z-x*N1z) + 
	=z*(1-z)*(N2z-x*N1z)
	=A(1+A)(C+xB)=A*(C+xB+AC+xAB)
	=AC+xAB+A^2C+xA^2B
	*/
	unique_ptr<MultiVarPoly> ConstraintsPolynomial(
		make_unique<MultiVarPoly>());
	MultiVarMonomial M;
	M.coeff = (x);
	M.vars = { 0, 0, 1 };
	ConstraintsPolynomial->AddMonomial(M); //xA^2B
	ConstraintsPolynomial->AddMonomial({ 0, 0, 2 }); //A^2C
	M.vars = { { 0, 1 } };
	ConstraintsPolynomial->AddMonomial(M); //xAB
	ConstraintsPolynomial->AddMonomial({ 0, 2 }); //AC

	/**
	* Construct witness
	*/
	UniPolynomialEvaluation eval;
	vector<FieldElement> points;

	FieldElement startPoint = x*x*x + x + ONE; //arbitrary initial value
	FieldElement endPoint = startPoint*FieldElement(NTL::GF2EInfo->power(::NTL::GF2E(x), POW2(subfieldDim) - 2));
	//   evaluation_t tMap(Algebra::compareFieldElements);
	FieldElement t = g;
	FieldElement nPt =startPoint;
	points.push_back(t);
#ifdef SANITY_CHECK
	eval.addPoint(t, nPt + ONE);
#else 
	eval.addPoint(t, nPt);
#endif
	cout << "filling evaluation" << endl;
	for (int i = 2; i<POW2(subfieldDim); ++i) {
		nPt = x * nPt;
		t = g * t;
		eval.addPoint(t, nPt);
		points.push_back(t);

	}

	//adding root at zero -not important it's root, just need to feel evaluation at 0 somehow, cause interpolation functions expects a subspace
	FieldElement ZERO = zero();
	eval.addPoint(ZERO, ONE);
	points.push_back(ZERO);
	UnivariatePolynomial P;
	cout << "interpolating" << endl;
	ExplicitLinearSubset H(BasisH.basis, BasisH.affineShift);
	P.interpolation(eval, BasisH.basis, H, BasisH.affineShift);
	unique_ptr<UnivariatePolynomialGeneral> P2(new UnivariatePolynomialGeneral());
	for (unsigned long int i = 0; i <= P.getDegree(); ++i) {//Ariel - had < rather than <=, but in any case better to use constructor as above
		FieldElement coeff = (P.getCoeff(i));
		P2->setCoefficient(i, coeff);
	}

	ACSPWitness witness(move(P2));

	/** We don't care about the input indicator, the input is empty, so it passes any way */
	ACSPInstance::boundaryConstraints_t boundaryConstraints;

	boundaryConstraints[FieldElement(g)] = FieldElement(startPoint);
	boundaryConstraints[Algebra::one()] = FieldElement(endPoint);

	/** Construct the instance data */
	ACSPInstance instance(contextField, move(vanishingSet), move(neighborPolys), move(ConstraintsPolynomial),witness.assignmentPoly().getDegree(), boundaryConstraints, compositionAlg);

	/** Return ACSP pair */
	return move(pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness>(move(instance), move(witness)));
}


TEST(ACSPWitnessChecker, repeatx){
	pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness> validPair = generate_repeatx();
	EXPECT_TRUE(PCP_Project::ACSPWitnessChecker::verify(validPair.first, validPair.second));
	
}

TEST(DISABLED_PCP4ACSP, repeatx){
	pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness> validPair = generate_repeatx();
	EXPECT_TRUE(PCP_Project::PCP_Prove_and_Verify_ACSP(validPair, PCP_Project::Prover::ProverAlgorithm(), PCP_Project::QueriesResFetcher::resultsFillingAlgorithm()));
}


TEST(ACSPWitnessChecker, repeatx8){
	pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness> validPair = generate_repeatx(8);
	EXPECT_TRUE(PCP_Project::ACSPWitnessChecker::verify(validPair.first, validPair.second));

}

TEST(DISABLED_PCP4ACSP, repeatx8){
	pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness> validPair = generate_repeatx(8);
	EXPECT_TRUE(PCP_Project::PCP_Prove_and_Verify_ACSP(validPair, PCP_Project::Prover::ProverAlgorithm(), PCP_Project::QueriesResFetcher::resultsFillingAlgorithm()));
}


TEST(ACSPWitnessChecker, DISABLED_repeatx16_VERY_SLOW){//disabled it cause it's extremeley slow right now
	pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness> validPair = generate_repeatx(16);
	EXPECT_TRUE(PCP_Project::ACSPWitnessChecker::verify(validPair.first, validPair.second));

}

TEST(DISABLED_PCP4ACSP, repeatx16){
	pair<PCP_Project::ACSPInstance, PCP_Project::ACSPWitness> validPair = generate_repeatx(16);
	EXPECT_TRUE(PCP_Project::PCP_Prove_and_Verify_ACSP(validPair, PCP_Project::Prover::ProverAlgorithm(), PCP_Project::QueriesResFetcher::resultsFillingAlgorithm()));
}
