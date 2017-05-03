/********************************************** Prover.cpp ************************************************/
/**
 * @file.
 *
 * The file Prover.cpp contains the implementation of PCP prover, starting from the RS= prover all the way up to proverLang.
 * 
 * For more information - Read the documentation of the header file Prover.hpp.
 */
  /*********************************************************************************************************/
#include "Prover.hpp"
#include "common/Utils/TaskReporting.hpp"
#include "PCP/PCP_common.hpp"
#include "common/Algebra/AlgebraTestUtils.hpp"
#include "common/Algebra/FFT/LDERunner.hpp"
#include "common/Algebra/LinearSpace.hpp"
#include "common/Algebra/ShiftedSubspacePolynomial.hpp"
#include <algebraLib/FFT.hpp>
#include <string>
#include <vector>
#include <omp.h>

using namespace std;
using namespace PCP_Project::PCP_common;
using Infrastructure::Log2;
using Infrastructure::Log2ceil;
using Infrastructure::POW2;
using Algebra::mapIntegerToFieldElement;
using Algebra::details::Basis;
using namespace Algebra::details;
using namespace Algebra;
using Algebra::LightUniPolyEval;
namespace PCP_Project {
namespace Prover{

/*****************************************************************************************/
/************************************ PRS= functions *************************************/
/*****************************************************************************************/

//anonymus namespace for BivariateExpansionProof
namespace{
class BivariateExpansionProof_basic : public BivariateExpansionProof{
public: 
    size_t getValIndex(const size_t rowId, const size_t colId)const{
        return colId + rowLen_*rowId;
    }
    
    size_t getIndexRow(const size_t valIndex)const{
        return valIndex / rowLen_;
    }
    
    size_t getIndexColumn(const size_t valIndex)const{
        return valIndex % rowLen_;
    }

    size_t size() const {
        return size_;
    }

    const MerkleTreeInterface& merkleTree()const {
        return *merkleTree_;
    }
    
    vector<FieldElement> getRow(const size_t rowId)const{
        vector<FieldElement> res;
        const size_t startIndex = getValIndex(rowId,0);
        for(size_t i=0; i< rowLen_; i++){
            res.push_back(evaluation_[startIndex+i]);
        }
        return res;
    }

protected:
    FieldElement const*const getProof()const{
        return evaluation_;
    }

    size_t rowLen_;
    size_t size_;
    
    //evaluation row by row
    FieldElement* evaluation_;
    unique_ptr<MerkleTree> merkleTree_;
};


//In next two classes is the actual heart of the construction is. The RS-PCPP recursion
//small and big differ in that in the big case not all proof is sto
//What does const * const mean?
class BivariateExpansionProof_small : public BivariateExpansionProof_basic{
public:
    BivariateExpansionProof_small(FieldElement const * const evaluation, const vector<FieldElement>& evaluationBasis, const FieldElement& affineShift){
        ///Splitting the received basis vectors into L0, L0' and L1:
       //const short splitingIndex = HALF_K(evaluationBasis.size()) - Gamma + 1;
		vector<FieldElement> L0_Basis(PCP_common::getL0Basis(evaluationBasis));
		vector<FieldElement> L0_Prime_Basis(PCP_common::getL0PrimeBasis(evaluationBasis));
		vector<FieldElement> L1_Basis(PCP_common::getL1Basis(evaluationBasis));

        ///Calculating the spaces spanned by L0' and L1:
        const size_t L0_size = POW2(L0_Basis.size());
        const size_t L0_Prime_size = POW2(L0_Prime_Basis.size());
        const size_t L1_size = POW2(L1_Basis.size());

        //resize the evaluation buffer to fit the proof
        rowLen_ = 2*L0_Prime_size;
        size_ = L1_size * rowLen_;
        evaluation_ = new FieldElement[size_];

        //construct proof
#ifndef _DEBUG
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
#endif
        for (long long i=0; i<L1_size; i++) {

			vector<FieldElement> L_Beta_Basis(PCP_common::getLBeta(evaluationBasis, i));
            const FieldElement beta = getSpaceElementByIndex(L1_Basis,affineShift,i);	//Note - beta is shifted by affineShift
            LDERunner(evaluation + i*L0_size, L0_Basis, beta, evaluation_ + (rowLen_*i), L_Beta_Basis, affineShift);
        }

        //merkle commitment tree construction
        merkleTree_ = unique_ptr<MerkleTree>(new MerkleTree(evaluation_,Log2(size_)));
    }

    ~BivariateExpansionProof_small(){
        delete[] evaluation_;
    }

};

//This method called for first level of recursion. small method called for next levels.
class BivariateExpansionProof_big : public BivariateExpansionProof_basic{
public:
    
    BivariateExpansionProof_big(const vector<FieldElement>& evaluationBasis, vector<FieldElement>&& evaluation) : 
    univariateEvaluation_(std::move(evaluation)), evaluationBasis_(evaluationBasis){
        ///Splitting the received basis vectors into L0, L0' and L1:
		vector<FieldElement> L0_Prime_Basis(PCP_common::getL0PrimeBasis(evaluationBasis));
		vector<FieldElement> L1_Basis(PCP_common::getL1Basis(evaluationBasis));

		///Calculating the spaces spanned by L0' and L1:
        const size_t L0_Prime_size = POW2(L0_Prime_Basis.size());
        const size_t L1_size = POW2(L1_Basis.size());

        //resize the evaluation buffer to fit the proof
        rowLen_ = 2*L0_Prime_size;
        size_ = L1_size * rowLen_;
        src_for_final_hash_.resize(2*L1_size);

        //calculate Merkle commitment
        addUnivariateSegmentToProof(&univariateEvaluation_[0],univariateEvaluation_.size(),0,evaluationBasis,Algebra::zero());//Question to Michael:Name Univariate confusing here - in this method you compute bivariate first level
        finishMerkleCommitment(evaluationBasis.size());
    }

    //segment size is expected to be a multiple of L0_size
    void addUnivariateSegmentToProof(FieldElement const * const segment, const size_t segmentLen, const short segmentIndex, const vector<FieldElement>& evaluationBasis, const FieldElement& affineShift){
		_COMMON_ASSERT(affineShift == Algebra::zero(), "new Code in PCPCommon assuming affineShift is always zero. Contact Ariel if this message appears");
        ///Splitting the received basis vectors into L0, L0' and L1:
        //const short splitingIndex = HALF_K(evaluationBasis.size()) - Gamma + 1;
		vector<FieldElement> L0_Basis(PCP_common::getL0Basis(evaluationBasis));
		vector<FieldElement> L0_Prime_Basis(PCP_common::getL0PrimeBasis(evaluationBasis));
		vector<FieldElement> L1_Basis(PCP_common::getL1Basis(evaluationBasis));

        ///Calculating the spaces spanned by L0' and L1:
        const size_t L0_size = POW2(L0_Basis.size());
        const size_t L0_Prime_size = POW2(L0_Prime_Basis.size());
        const size_t L1_size = POW2(L1_Basis.size());

        //segment specific constants
        const size_t rowsPerSegment = segmentLen/L0_size;
        const size_t startIndex = rowsPerSegment*segmentIndex;
#ifndef _DEBUG
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
#endif
        for (long long i=0; i<rowsPerSegment; i++) {
            const size_t globalIndex = i+startIndex;

			vector<FieldElement> L_Beta_Basis = PCP_common::getLBeta(evaluationBasis,globalIndex);
            const FieldElement beta = getSpaceElementByIndex(L1_Basis,affineShift,globalIndex);	//Note - beta is shifted by affineShift. globalIndex is beta's index in L1


    //        if (globalIndex < POW2(Mu)) {
    //            L_Beta_Basis_tmp.insert(L_Beta_Basis_tmp.begin(),evaluationBasis[firstIndexOutsideL0_Prime(evaluationBasis.size())]);
    //        }
    //        else {
				//L_Beta_Basis_tmp.insert(L_Beta_Basis_tmp.begin(),beta + affineShift);	//beta is already shifted, so we cancel by adding affineShift again.beta added in front to make sure all rows in next recursion level have intersection with original function
    //        }
    //       

            vector<FieldElement> rowData(rowLen_);
            LDERunner(segment + i*L0_size, L0_Basis, beta, &rowData[0], L_Beta_Basis, affineShift);
            MerkleTree currRowTree(&rowData[0],L_Beta_Basis.size());
            src_for_final_hash_[globalIndex*2] = currRowTree.tree_[2];
            src_for_final_hash_[globalIndex*2+1] = currRowTree.tree_[3];
        }
    }

	//here prover computes and commits 
    void finishMerkleCommitment(const short evaluationBasisSize){
        short logSize = Infrastructure::Log2(src_for_final_hash_.size()*sizeof(CryptoCommitment::hashDigest_t)/sizeof(FieldElement));
        merkleTree_ = unique_ptr<MerkleTree>(new MerkleTree((FieldElement*)&(src_for_final_hash_[0]) , logSize));
        src_for_final_hash_.resize(0);
    }
    
    vector<FieldElement> getRow(const size_t rowId)const{
        const FieldElement affineShift = Algebra::zero();
        ///Splitting the received basis vectors into L0, L0' and L1:
        //const short spliti = HALF_K(evaluationBasis_.size()) - Gamma + 1;
		vector<FieldElement> L0_Basis(PCP_common::getL0Basis(evaluationBasis_));
		vector<FieldElement> L0_Prime_Basis(PCP_common::getL0PrimeBasis(evaluationBasis_));
		vector<FieldElement> L1_Basis(PCP_common::getL1Basis(evaluationBasis_));

        ///Calculating the spaces spanned by L0' and L1:
        const size_t L0_size = POW2(L0_Basis.size());
        const size_t L0_Prime_size = POW2(L0_Prime_Basis.size());
        const size_t L1_size = POW2(L1_Basis.size());
        
        //calculate row
        vector<FieldElement> res(POW2(L0_Prime_Basis.size()+1));
        {
			vector<FieldElement> L_Beta_Basis = PCP_common::getLBeta(evaluationBasis_, rowId);
			const FieldElement beta = getSpaceElementByIndex(L1_Basis, affineShift, rowId);	//Note - beta is shifted by affineShift

            //vector<FieldElement> L_Beta_Basis_tmp(L0_Prime_Basis);//Question to Michael:we're copying L0_Prime_Basis many times. Could be avoided if L_Beta_Basis_tmp defined outside of {}
            
            //if (rowId < POW2(Mu)) {//this condition checks whether beta is already in L0_Prime, in which case we add a different element than beta
            //    L_Beta_Basis_tmp.insert(L_Beta_Basis_tmp.begin(),evaluationBasis_[firstIndexOutsideL0_Prime(evaluationBasis_.size())]);
            //}
            //else {
            //    L_Beta_Basis_tmp.insert(L_Beta_Basis_tmp.begin(),beta + affineShift);	//beta is already shifted, so we cancel by adding affineShift again.beta added in front to make sure all rows in next recursion level have intersection with original function
            //}
           LDERunner(&univariateEvaluation_[rowId*L0_size], L0_Basis, beta, &res[0], L_Beta_Basis, affineShift);
        }

        return res;

    }

    const vector<FieldElement>& getUnivariateEvaluation()const{
        return univariateEvaluation_;
    }

    ~BivariateExpansionProof_big(){
    }

private:
    vector<CryptoCommitment::hashDigest_t> src_for_final_hash_;
    vector<FieldElement> univariateEvaluation_;
    vector<FieldElement> evaluationBasis_;
};
}

/**
 * A single step of a RS PCPP proof
 */
unique_ptr<BivariateExpansionProof> ProverAlgorithm::expandUnivariateToBivariate(FieldElement const * const evaluation, const vector<FieldElement>& evaluationBasis, const FieldElement& affineShift)const{
    return unique_ptr<BivariateExpansionProof>(new BivariateExpansionProof_small(evaluation,evaluationBasis,affineShift));
}
/**
 * Generation of low degree (aka RS) proof for a row of a bivariate expansion
 */
unique_ptr<BivariateExpansionProof> ProverAlgorithm::proofRow_lowDegree(const BivariateExpansionProof& src, const size_t row_index ,const vector<FieldElement>& rowBasis, const FieldElement& affineShift)const{
    const size_t rowLen = POW2(rowBasis.size());
    vector<FieldElement> row(rowLen);
#ifndef _DEBUG
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
#endif
    for(long long i=0; i< rowLen; i++){
        row[i] = src.readValue(row_index,PCP_common::shuffledRowIndexToOriginal(i,rowLen));//as recursive rows basis needs to be shuffled when constructing subproof , but src holds row according to basis before shuffling need to translate indicies between bases
    }

    return expandUnivariateToBivariate(&row[0],rowBasis,affineShift);
}

/**
 * Generation of low degree (aka RS) proof for a column of a bivariate expansion
 */
unique_ptr<BivariateExpansionProof> ProverAlgorithm::proofColumn_lowDegree(const BivariateExpansionProof& src, const size_t column_index ,const vector<FieldElement>& columnBasis, const FieldElement& affineShift)const{
    const size_t columnLen = POW2(columnBasis.size());
    vector<FieldElement> column(columnLen);
#ifndef _DEBUG
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
#endif
    for(long long i=0; i< columnLen; i++){
        column[i] = src.readValue(i,column_index);
    }

    return expandUnivariateToBivariate(&column[0],columnBasis,affineShift);
}

class cachedPolynomial : public UnivariatePolynomialInterface{
public:
    cachedPolynomial(const UnivariatePolynomialInterface& poly, const BivariateExpansionProof_big& cache, const vector<FieldElement>& proofBasis, const FieldElement& proofShift, const ACSPInstance& instance):
        poly_(poly), cache_(cache), proofBasis_(proofBasis), proofShift_(proofShift){
    
        using Algebra::details::UnivariatePolynomial;
        using Algebra::details::UniPolynomialEvaluation;

        UnivariatePolynomial Ax;
        UnivariatePolynomial Z_X;
        {
            vector<FieldElement> S_x;
            UniPolynomialEvaluation f_eval;
            for (const auto& p : instance.boundaryConstraints()) {
                f_eval.addPoint(p.first, p.second);
                S_x.push_back(p.first);
            }
            Ax.interpolation(f_eval, S_x);
            Z_X.findNonSparseVanishingPoly(S_x);
        }

        boundaryPoly = UnivariatePolynomialGeneral(Ax.getCoeffsVec());
        boundaryVanishingPoly = UnivariatePolynomialGeneral(Z_X.getCoeffsVec());
    }
    
    FieldElement eval(const FieldElement& x)const{
        return poly_.eval(x);
    }
	
    PolynomialDegree getDegree()const{
        return poly_.getDegree();
    }
    
    unique_ptr<PolynomialInterface> clone()const{
        return poly_.clone();
    }
	
    FieldElement getCoefficient(const unsigned int index)const{
        return poly_.getCoefficient(index);
    }
    
    vector<FieldElement> eval(const vector<FieldElement>& orderedBasis, const FieldElement& shift)const{
        if(shift != proofShift_){
            return poly_.eval(orderedBasis,shift);
        }
        
        for(short i=0; i<orderedBasis.size();i++){
            if( orderedBasis[i] != proofBasis_[i]){
                return poly_.eval(orderedBasis,shift);
            }
        }


        //return from catch
        {
        /*const short cosetSizeLog = HALF_K(proofBasis_.size()) - Gamma + 1;
        const size_t numCommonColumns = POW2(cosetSizeLog + Mu);*/
		const size_t numCommonColumns = POW2(dimOfL0_Prime(proofBasis_.size()));
        vector<FieldElement> res(POW2(orderedBasis.size()));
#ifndef DEBUG
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
#endif
        for(long long i=0; i<res.size(); i++){
            const FieldElement srcPoint = getSpaceElementByIndex(orderedBasis,shift,i);
            const FieldElement val_cache =  cache_.getUnivariateEvaluation()[i];
            const FieldElement polyVal = val_cache*boundaryVanishingPoly.eval(srcPoint) + boundaryPoly.eval(srcPoint);
            res[i] = polyVal;
        }
        return res;
        }
    }
private:
    const UnivariatePolynomialInterface& poly_;
    const BivariateExpansionProof_big& cache_;
    const vector<FieldElement> proofBasis_;
    const FieldElement proofShift_;
    UnivariatePolynomialGeneral boundaryPoly;
    UnivariatePolynomialGeneral boundaryVanishingPoly;
};

/*****************************************************************************************/
/********************************* P_ACSP functions **************************************/
/*****************************************************************************************/
proof_t ACSP_Producer(const std::pair<ACSPInstance, ACSPWitness>& ACSPPair) {

    using Algebra::zero;
    using std::unique_ptr;
    using std::cout;
    using std::endl;

    proof_t result;

	/**************************************/
	/** Step 1 - Parameter Instantiation: */
	/**************************************/

	const ACSPInstance& acspInstance = ACSPPair.first;
	const Algebra::UnivariatePolynomialInterface& witnessPoly = ACSPPair.second.assignmentPoly();

    //
    //Construct ordered basis for proofs
    //
    const vector<FieldElement> basisPCPP(PCP_common::basisForPCPP(acspInstance));
	

	/**************************************************************/
	/*** Step 2 - Prove witness A is low degree and satisfies boundary constraints (by *single* VRS PCPP)
	/**************************************************************/

    {
    TASK("Constructing boundary constraints proof");
	cout << "Consistency with the instance: |x| = " << acspInstance.boundaryConstraints().size() << ", |S_x| = " << acspInstance.boundaryConstraints().size() << endl;
    
    unique_ptr<Algebra::UnivariatePolynomialInterface> polyToTest;
    std::vector<FieldElement> S_x;
    {
        TASK("Constructing formal sum (division of polynomial)");
        ///Compute the evaluation f:Sm->{0,1} such that forall i=1..|x| f(alpha_i)=x_i
        ///where alpha_i is the i'th element of Sm and x_i is the i'th bit of the instance x.
        ///And writing the input location to the vector S_x
        Algebra::details::UniPolynomialEvaluation f_eval;
        for (const auto& p : acspInstance.boundaryConstraints()) {
            f_eval.addPoint(p.first, p.second);
            S_x.push_back(p.first);
        }

        ///Generating Ax - the low-degree extension of f over S_x, and calling the VRS prover:
        Algebra::details::UnivariatePolynomial Ax;
        Ax.interpolation(f_eval, S_x);
        
        Algebra::details::UnivariatePolynomial Z_X;
        Z_X.findNonSparseVanishingPoly(S_x);
        
        const auto deltaPoly = Algebra::UnivariatePolynomialGeneral(Ax.getPoly()) - Algebra::UnivariatePolynomialGeneral(witnessPoly);

#ifdef IDENTITY_WITNESS
		polyToTest.reset(new Algebra::UnivariatePolynomialGeneral(std::vector<FieldElement>({ Algebra::zero(), Algebra::one() })));
#else
		polyToTest = Algebra::UnivariatePolynomialGeneral(Z_X.getPoly()).divideByMe(deltaPoly);	
#endif //IDENTITY_WITNESS

    }

	///Step 1 - Produce the RS proof PI for A-Ax:
	cout << "Calling Input Consistency VRS proof" << endl;
    {
        TASK("Constructing proof");
        BivariateExpansionProof_big* proof = new BivariateExpansionProof_big(basisPCPP,polyToTest->eval(basisPCPP,zero()));
        result.boundaryRS_proof = unique_ptr<BivariateExpansionProof>(proof);//Interesting. how does he allow you to do this with unique pointer?
    }
    }


	/***************************************************************/
	/** Step 3 - Proof that A satisfies the constraint polynomial:  */
	/***************************************************************/

	//The proof works by computing composition of A with constraint polynomial, dividing by Z_H and using PCPP to prove low degreeness
	//question to Michael: Where did the consistnecy queries go?
    {
	TASK("Computing B and p1. Size of p0 = " + to_string(POW2(basisPCPP.size())));
    const cachedPolynomial cachedWitness(witnessPoly,*((BivariateExpansionProof_big*)result.boundaryRS_proof.get()),basisPCPP,zero(),acspInstance);
	//computing composition with witness
    unique_ptr<UnivariatePolynomialInterface> polyToTest = acspInstance.composeWithWitness_and_divideByVanishingSpacePoly(cachedWitness);
    
    cout<<"deg P(A(N))/Z_H =  "<< PolynomialDegree::integral_t(polyToTest->getDegree())<<endl;
    {
        TASK("Constructing proof");
        BivariateExpansionProof_big* proof = new BivariateExpansionProof_big(basisPCPP,polyToTest->eval(basisPCPP,zero()));
        result.compositionRS_proof = unique_ptr<BivariateExpansionProof>(proof);
    }
        
    }

    return move(result);
}


/*****************************************************************************************/
/********************************* P_Lang functions **************************************/
/*****************************************************************************************/

proof_t ProverAlgorithm::PCP_Lang_Producer(const std::pair<ACSPInstance,ACSPWitness>& ACSPPair)const {

	///Using witness reduction to generate the univariate witness polynomial A:
	cout << "Reduced witness A is of degree: " << PolynomialDegree::integral_t(ACSPPair.second.assignmentPoly().getDegree()) << endl;

	///Calling ACSP prover:
	return ACSP_Producer(ACSPPair);
}


} //namespace Prover
}	//of namespace PCP_Project
