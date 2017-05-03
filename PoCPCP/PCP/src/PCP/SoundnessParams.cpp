/******************************************* SoundenessParams.cpp *********************************************/
/**
 * @file.
 *
 * The file SoundenessParams.cpp contains the implementation of the parameters relevant to the soundness analysis, i.e. 
 * the parameters that are used to compute the amplified verifier's number of repetitions.
 * 
 * For more information - Read the documentation of the header file SoundenessParams.hpp.
 */
  /************************************************************************************************************/
#include "SoundnessParams.hpp"
#include "common/Infrastructure/Infrastructure.hpp"
#include <assert.h>
#include <algorithm>
#include <PCP/PCP_common.hpp>
using namespace std;
using namespace NTL;
using Infrastructure::POW2;
//#define FASTER_DEBUGGING //added by ARIEL  --> shallow recursion --> faster running when debugging
namespace PCP_Project {

/** Default initialization of the soundness parameters. */
int SoundnessParameters::mu;						///L0' has mu more dimensions than L0.
int SoundnessParameters::eta;						///The degree we test in RS= is |L|/(2^eta)-1, e.g. |L|/8-1 for eta=3.
int SoundnessParameters::gamma;						///L0 contains [k/2] basis elements minus gamma.
ProximityParam SoundnessParameters::deltaVRS;		///The proximity parameter for AMP-VRS
bool SoundnessParameters::useFeasible;
bool SoundnessParameters::useFixedDepth;
int SoundnessParameters::recursionDepth;
int SoundnessParameters::securityParam;
int SoundnessParameters::consistencyRepNum;
/** Sets the default values for the soundness parameters */
void SoundnessParameters::setupSoundnessParams() {
#ifndef FASTER_DEBUGGING

	gamma = 1;
    mu = 1;								///L0' has mu more dimensions than L0.
	eta = 3;							///An upper bound on The degree we test in RS is |L|/(2^eta)-1, e.g. |L|/8-1 for eta=3.
	deltaVRS = 1.0 / 8;					///The proximity parameter for AMP-VRS
	securityParam = 80; 						///log, base 0.5, of soundness error we want
	useFeasible =true;
	recursionDepth = 3;
	useFixedDepth = true;
	consistencyRepNum = std::ceil(securityParam/(double)eta);//each ACSP consistency check misses catching error w.p 2^{-eta}

#else
    gamma = 1;
	mu =2;								///L0' has mu more dimensions than L0.
	eta = 3;							///The degree we test in RS= is |L|/(2^eta)-1, e.g. |L|/8-1 for eta=3.
	deltaVRS = 1.0 / 8;						///The proximity parameter for AMP-VRS
	securityParam = 1; 						///log, base 0.5, of soundness error we want

#endif
}

/**
 * The function computes the number of repetitions done by AMP-VRS verifiers.
 * The formulas are according to the paper
 * @param k is the size of the basis.
 * @param delta is the parameter we want to test proximity to (need to get this as parameter as can depend on number of ACSP neighbors)
 */
unsigned long SoundnessParameters::computeRepetitions(const int k, const ProximityParam& delta) {
//#define SINGLE_QUERY
#ifdef SINGLE_QUERY
	return 1;
#else
	double alpha0 = (double)(pow(2, eta) - 1) / (2 * pow(2, eta) - 1);//the best probability to choose a column in first recursion level assuming row and column compliant attacks analysis on random functions
	double alpha1 = (double)(pow(2, mu+1) - 1) / (2 * pow(2, mu+1) - 1);//the best probability to choose a column in subsequent recursion levels assuming row and column compliant attacks analysis on random functions

	
	//Problem: There is no taking into account of the fact that delta of witness need to be multiplied by 1/q, q being # of neighb.
	//when using conjectured soundness from the paper, we assume a depth d 
	double rejProb; //the proven or conjectured rejection probability at each round
	unsigned long repNum;
	if (useFeasible)//Based on section 6.4 on paper, assuming  function is as bad as 1/x, ignoring rather negligible 1/|F| factor
		           //TODO: update explanation of this new computation
		            
		rejProb = alpha0 * pow(alpha1,RECURSION_DEPTH-1);

	
	else {
		if (RECURSION_DEPTH == 1)//based on lemma 4.12 - soundness of depth 1 test in the paper

			rejProb = std::min(delta, 0.5);
		if (RECURSION_DEPTH == 2)//based on lemma 4.15 Soundness of depth 2 test
			rejProb = std::min((delta / 20), 0.005);//TODO:this bound uses the somewhat smarter test using rows and extended rows not yet implemented
													//assumes eta>= 3
		if (RECURSION_DEPTH > 2)
			_COMMON_FATAL("no concrete analysis written for depth>2, not useful for any reasonable input length any way");
	}
	//computing number of reps needed to get acc prob down to 2^{-securityParam}
	repNum = std::ceil((double)securityParam* (std::log(0.5)/std::log(1 - rejProb) ));
	/*cout << "Number of repetitions is " << repNum << endl;
	cout << "alpha0 is " << alpha0 << endl;
	cout << "alpha1 is " << alpha1 << endl;
	cout << "rejprob is " << rejProb << endl;
*/
	return repNum;
#endif
}

//choose whether to recurse on row or column to maximize rejection probability (currently againt row\col compliant attackes on random functions)
queryType_t SoundnessParameters::chooseRowOrCol(int depth){
	_COMMON_ASSERT(FIXED_DEPTH == true, "PCP currently assumed to run with fixed recursion depth. Tell Ariel if this appears");
	size_t r = rand();
	if (depth == RECURSION_DEPTH){//on first recursion level, choose column with prob 2^eta-1/2*2^eta-1
		size_t a = (r % (2*POW2(eta)-1))+1;
		if (a <= POW2(eta) - 1)
			return COLUMN;
		else
			return ROW;
	}
	else
	{//on next levels, choose column with prob 2^(mu+1)-1/2*2^(mu+1)-1
		size_t a = (r % (2 *POW2(mu+1) - 1)) + 1;
		if (a <= POW2(mu+1) - 1)
			return COLUMN;
		else
			return ROW;
	}
		
}



//setting recursion depth according to dimension of PCPP space.
//if dim is k, then recursion depth should be largest d s.t. every recursive row has intersection with original poly (see Remark 4.11 in paper)
//This is computed by each time increasing d, and reducing dim to dim of L_Beta, until we are left with something smaller than 2d
//(See remark  in paper for rationale)
void SoundnessParameters::setRecursionDepth(const short& basisSize){
	_COMMON_ASSERT(FIXED_DEPTH == true, "PCP currently assumed to run with fixed recursion depth. Tell Ariel if this appears");
	_COMMON_ASSERT(Mu == 1, "current computation here assumes Mu =1. Tell Ariel if this appears");
	int d = 1;
	auto k = PCP_common::dimOfLBeta(basisSize);
	while (k > 2*d){
		k = PCP_common::dimOfLBeta(k);
		d++;
	}
	recursionDepth = d-1;
}

//Ariel:commented out version  below, which followed concrete paper. New version follows the paper
bool areSoundnessParametersValid(){

	return ((FIXED_DEPTH) && (Mu == 1) &&
		(USE_FEASIBLE || (FIXED_DEPTH && RECURSION_DEPTH <= 2))&&
		(Eta >= 3)
		);
	//bounds from paper for concrete section are stated only for depth 
	//1 and 2 - and would be better for 3, only for dimension about 80
	// The paper deals only with Mu =1 for simplicity, and as this seems to be always the best choice 
	// Simpler for now to run PCP to fixed depth, rather than down to certain base case dimension
	// so enforcing that
	//current setting of params, (only in concrete case), in ComputeRepetitions assume Eta>=3
}
}	//Of namespace PCP_Project
