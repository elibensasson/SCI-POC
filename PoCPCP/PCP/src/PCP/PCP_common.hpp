#ifndef PCP_COMMON_HPP_
#define PCP_COMMON_HPP_
//#define SANITY_CHECK //ruins witness poly in ACSP_Producer to check verifer won't accept corrupted proof
#define SET_RECURSION_DEPTH //set recursion depth so that all recursive rows have intersection with original poly rather than manually
#include "SoundnessParams.hpp"
#define NO_ZK // a flag for certain lines that might need change if we add Zero Knowledge
#define NO_RSLESS // `round up' degree to 2^t-1 and avoid using RS_LESS, which excutes 2 PCPs instead of 1. Eli claims this doesn't hurt soundness
//#define	 DUMMY_EVALS//Just evaluate the composition poly as all 0 for faster debugging - also comments out p0p1 consistency check in decision alg
//#define ONLY_INPUT_CONSISTENCY //omitting other parts to check bugs in input consistency part faster
#include "languages/ACSP/ACSPInstance.hpp"
#include "common/Algebra/details/FiniteFields.hpp"
#include "common/Algebra/AlgebraCommon.hpp"
namespace PCP_Project {
//#define IDENTITY_WITNESS //sometimes for debugging it's convenient to use the witness (after divided by boundary) to be A(X)=X
//#define ZERO_COMPOSITION //sometimes for debugging it's convenient to evaluate composition is identically zero//
//#define ARIEL_DEBUG // comments out various things that are inconvenient for Ariel when debugging

namespace PCP_common {
	

Algebra::PolynomialDegree composition_div_ZH_degreeBound(const ACSPInstance& src);

Algebra::PolynomialDegree witness_div_Z_Boundery_degreeBound(const ACSPInstance& src);

/**
 * Returns the basis over which the low degree tests are done
 * (aka RS PCPP)
 */
std::vector<Algebra::FieldElement> basisForPCPP(const ACSPInstance& src);

/**
 * Returns the basis over which the consistency of the witness and the composition
 * polynomial is proved
 */
std::vector<Algebra::FieldElement> basisForConsistency(const ACSPInstance& src);

//The three methods below are to return the relevant subbases of L in the PCPP recursion
vector<Algebra::FieldElement> getL0Basis(const vector<Algebra::FieldElement>& L);
vector<Algebra::FieldElement> getL0PrimeBasis(const vector<Algebra::FieldElement>& L);
vector<Algebra::FieldElement> getL1Basis(const vector<Algebra::FieldElement>& L);


//This method returns the first basis index of L, not in L0'. This is the index of the element needed to be added to LBeta, when Beta is already contained in L_0'
Algebra::details::BasisIndex firstIndexOutsideL0_Prime(const Algebra::details::BasisIndex k);
//the dim LBeta should have if dim(L)=k
Algebra::details::BasisIndex dimOfLBeta(const Algebra::details::BasisIndex k);
//the dim L0' should have if dim(L)=k
Algebra::details::BasisIndex dimOfL0_Prime(const Algebra::details::BasisIndex k);
//the dim L1 should have if dim(L)=k

//LBeta basis order needs to be changed before recursion, to ensure all recursive rows have intersection with original poly. See Remark in the paper in section about 
// BSS extensions of depth greater than one
vector<Algebra::FieldElement> changeRowBasisBeforeRecursion(const vector<Algebra::FieldElement>& LBeta);

//generates Lbeta. 
//Currently assumes that the intersectin of L'0 and L1 is in the first 2^mu indices of L1
vector<Algebra::FieldElement> getLBeta(const vector<Algebra::FieldElement>& BasisL, const Algebra::details::BasisIndex betaIndex);


Algebra::details::BasisIndex L0Dim(const Algebra::details::BasisIndex k);
inline Algebra::details::BasisIndex lastIndexOfL0(const Algebra::details::BasisIndex k){
	return L0Dim(k) - 1;
}
inline Algebra::details::BasisIndex dimOfL1(const Algebra::details::BasisIndex k){
	return k - L0Dim(k);
}
Algebra::details::BasisIndex shuffledRowIndexToOriginal(const Algebra::details::BasisIndex i, const Algebra::details::BasisIndex& rowLength);
vector<Algebra::FieldElement> changeRowValuesBeforeRecursion(const vector<FieldElement>& row);

} //namespace PCP_common

} //namespace PCP_Project

#endif //PCP_COMMON_HPP_
