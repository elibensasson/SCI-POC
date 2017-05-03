#include "PCP_common.hpp"
#include <NTL/GF2XFactoring.h>

namespace PCP_Project {

namespace PCP_common {
	
using Algebra::FiniteSetInterface;
using Algebra::FieldElement;
using Algebra::PolynomialDegree;
using Algebra::PolynomialInterface;
using Algebra::mapIntegerToFieldElement;
using Infrastructure::Log2;
using Algebra::details::BasisIndex;
using std::max;
using std::ceil;

PolynomialDegree composition_div_ZH_degreeBound(const ACSPInstance& src){
    //
    // calculate degree bound of composition polynomial divided by the vanishing space poly
    //
    const FiniteSetInterface& vanishingSet = src.vanishingSet();
	const ACSPInstance::polynomialsVec& neighborPolys = src.neighborPolys();
	const PolynomialInterface& constraintPoly = src.constraintPoly();

	vector< PolynomialDegree> inputDegrees;
	inputDegrees.push_back(PolynomialDegree(1));
	for (int i = 0; i < neighborPolys.size(); i++){
		inputDegrees.push_back(neighborPolys[i]->getDegreeBound(src.witnessDegreeBound()));
	}

	const auto compositionDegree = PolynomialDegree::integral_t(constraintPoly.getDegreeBound(inputDegrees)) - vanishingSet.size();
    
    return PolynomialDegree(compositionDegree);
}

PolynomialDegree witness_div_Z_Boundery_degreeBound(const ACSPInstance& src){
    return PolynomialDegree(PolynomialDegree::integral_t(src.witnessDegreeBound()) - src.boundaryConstraints().size());
}

/**
 * Returns the basis over which the low degree tests are done
 * (aka RS PCPP)
 */
vector<FieldElement> basisForPCPP(const ACSPInstance& src){
	
    //
    // find the space size, it should be at least \f$ 2^{\Eta} \f$ times bigger
    // than the maximal degree we will construct PCPP for
    //
    const PolynomialDegree maxDegree = max(max(composition_div_ZH_degreeBound(src),witness_div_Z_Boundery_degreeBound(src)),PolynomialDegree(2*Gamma));
    const short maxDegreeLog = ceil(Log2(PolynomialDegree::integral_t(maxDegree)));
    const short rankForPCPP_space = maxDegreeLog+Eta;

	//setting recursion depth according to dimension of PCPP space.
	//if dim is k, then recursion depth should be largest d s.t. 2d2^d<k. This will ensure every recursive row has intersection with original poly (see Remark 4.11 in paper)
#ifdef SET_RECURSION_DEPTH
	SoundnessParameters::setRecursionDepth(rankForPCPP_space);
#endif

    //
    //Construct and return the standard basis of needed size
    //
    vector<FieldElement> res;
    while(res.size() < rankForPCPP_space){
        res.push_back(mapIntegerToFieldElement(res.size(),1,1));
    }
    return res;
}

/**
 * Returns the basis over which the consistency of the witness and the composition
 * polynomial is proved
 */
std::vector<Algebra::FieldElement> basisForConsistency(const ACSPInstance& src){
    auto res = basisForPCPP(src);
    res.resize(res.size()-1);
    return res;
}


//Helper function for methods below
BasisIndex L0Dim(const BasisIndex k){
	double test = (double)k / (double)2 - Gamma;
	BasisIndex res = ceil(test);
	return res;
}

vector<Algebra::FieldElement> changeRowBasisBeforeRecursion(const vector<Algebra::FieldElement>& LBeta){
    vector<FieldElement> res(LBeta.begin(), LBeta.begin() + LBeta.size() - Mu - 1);
	for (auto i = 1; i <= Mu + 1; i++)
		res.insert(res.begin(),LBeta[LBeta.size() - i]);
	return res;

}

//translating i in changed row basis for recursion to it's value in origina LBeta basis
Algebra::details::BasisIndex shuffledRowIndexToOriginal(const Algebra::details::BasisIndex i, const BasisIndex& rowLength){
    auto d = Infrastructure::POW2(Mu + 1);
	//if (i == 0) return 0;
	auto j = (i%d) * rowLength / d;
	auto k = std::floor(i / d);
	return j+k;
}
vector<Algebra::FieldElement> changeRowValuesBeforeRecursion(const vector<FieldElement>& row){
	vector<FieldElement> res;
	for (auto i = 0; i < row.size(); i++)
		res.push_back(row[shuffledRowIndexToOriginal(i,row.size())]);
	return res;
}



//retrieves beta's index in L given its index in L1, assuming L's dimension is k
BasisIndex  L1IndexToLIndex(const Algebra::details::BasisIndex betaIndex, const BasisIndex k){
	BasisIndex index = Infrastructure::POW2(L0Dim(k))*betaIndex;
	return index;
}

//methods duplicated cause in Prover.cpp using vecotr of FieldElement rather than Basis object
vector<Algebra::FieldElement> getL0Basis(const vector<Algebra::FieldElement>& BasisL){
	return vector<FieldElement>(BasisL.begin(), BasisL.begin() + L0Dim(BasisL.size()));

}
vector<Algebra::FieldElement> getL0PrimeBasis(const vector<Algebra::FieldElement>& BasisL){
	auto d = L0Dim(BasisL.size());
	return vector<FieldElement>(BasisL.begin(), BasisL.begin() + d+Mu);


}
vector<Algebra::FieldElement> getL1Basis(const vector<Algebra::FieldElement>& BasisL){
	auto d = L0Dim(BasisL.size());
	return vector<FieldElement>(BasisL.begin() + d, BasisL.end());

}


//This method returns the first basis index of L, not in L0'. This is the index of the element needed to be added to LBeta, when Beta is already contained in L_0'
Algebra::details::BasisIndex firstIndexOutsideL0_Prime(const Algebra::details::BasisIndex k){
	return L0Dim(k) + Mu;//there is a cancelation of -1 cause indices start from 0, and +1 to get outside of L0'
}


//generates Lbeta. 
//Currently assumes that the intersectin of L'0 and L1 is in the first 2^mu indices of L1
vector<Algebra::FieldElement> getLBeta(const vector<Algebra::FieldElement>& BasisL, const Algebra::details::BasisIndex betaIndex){
	vector<FieldElement> res = getL0PrimeBasis(BasisL);
	BasisIndex betaIndexInL = L1IndexToLIndex(betaIndex, BasisL.size());
	if (betaIndex < Infrastructure::POW2(Mu))
		res.push_back(BasisL[firstIndexOutsideL0_Prime(BasisL.size())]);
	else{
		FieldElement beta = Algebra::getSpaceElementByIndex(BasisL, Algebra::zero(), betaIndexInL);
		res.push_back(beta);
	}
	return res;
}


//the dim LBeta should have if dim(L)=k
Algebra::details::BasisIndex dimOfLBeta(const Algebra::details::BasisIndex k){
	return L0Dim(k) + Mu + 1;
}
//the dim L0' should have if dim(L)=k
Algebra::details::BasisIndex dimOfL0_Prime(const Algebra::details::BasisIndex k){
	return L0Dim(k) + Mu;
}

} //namespace PCP_common

} //namespace PCP_Project
