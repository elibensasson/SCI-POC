#include "Prover_Dummy.hpp"
#include "PCP/PCP_common.hpp"

namespace PCP_Project {
namespace Prover{

using CryptoCommitment::hashDigest_t;
using CryptoCommitment::path_t;
using Algebra::FieldElement;
using Algebra::generateRandom;
using Infrastructure::POW2;
using Infrastructure::Log2;
using std::vector;
using std::unique_ptr;
using std::move;

namespace{

 
class MerkleDummy : public MerkleTreeInterface{
public:
MerkleDummy(short src_log_size):logSizeBytes_(src_log_size+Log2(sizeof(FieldElement))), c_(rand()), root_(getHVal(1)){};
const CryptoCommitment::hashDigest_t& getRoot()const{
    return root_;
}

path_t getPathToBlock(const size_t blockIndex)const{
    
    path_t result;
    const short firstRow_logLen = logSizeBytes_ - CryptoCommitment::logBytesPerHash - 1;
    size_t curr_offset = (1UL<<firstRow_logLen) + (blockIndex>>1);
    
    while(curr_offset > 1UL){
        result.push_back(getHVal(curr_offset ^ 1UL));
        curr_offset >>= 1;
    }
    
    return result;

}
private:
size_t c_;
short logSizeBytes_;
hashDigest_t hd_;
hashDigest_t root_;

hashDigest_t getHVal(size_t indx)const{
    hashDigest_t res = hd_;
    ((size_t*)(&res))[0] = c_;
    ((size_t*)(&res))[1] = indx;
    
    return res;
}
};



class BivariateDummy : public BivariateExpansionProof{
public:
    
    BivariateDummy(const vector<FieldElement>& evaluationBasis):
        rowLen_(POW2(1+PCP_common::getL0PrimeBasis(evaluationBasis).size())),
        size_(POW2(1+PCP_common::getL0PrimeBasis(evaluationBasis).size() + PCP_common::getL1Basis(evaluationBasis).size())),
        merkle_(1+PCP_common::getL0PrimeBasis(evaluationBasis).size() + PCP_common::getL1Basis(evaluationBasis).size()){};

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
    
    FieldElement readValue(const size_t idx)const{
        return generateRandom();
    }
    
    vector<FieldElement> getRow(const size_t rowId)const{
        vector<FieldElement> res;
        for(size_t i=0; i< rowLen_; i++){
            res.push_back(generateRandom());
        }
        return res;
    }
    
    const MerkleTreeInterface& merkleTree()const{
        return merkle_;
    }
protected:
    FieldElement const*const getProof()const{
        _COMMON_FATAL("Not implemented");
    }
private:

const MerkleDummy merkle_;
const size_t size_;
const size_t rowLen_;
};

} //anonymous namespace

proof_t Prover_Dummy::PCP_Lang_Producer(const std::pair<ACSPInstance,ACSPWitness>& ACSPPair)const {
    
    proof_t result;
    const vector<FieldElement> basisPCPP(PCP_common::basisForPCPP(ACSPPair.first));
    result.boundaryRS_proof = unique_ptr<BivariateExpansionProof>(new BivariateDummy(basisPCPP));
    result.compositionRS_proof = unique_ptr<BivariateExpansionProof>(new BivariateDummy(basisPCPP));

    return move(result);
}

unique_ptr<BivariateExpansionProof> Prover_Dummy::expandUnivariateToBivariate(FieldElement const * const evaluation, const vector<FieldElement>& evaluationBasis, const FieldElement& affineShift)const{
    return unique_ptr<BivariateExpansionProof>(new BivariateDummy(evaluationBasis));
}
    
unique_ptr<BivariateExpansionProof> Prover_Dummy::proofRow_lowDegree(const BivariateExpansionProof& src, const size_t row_index ,const vector<FieldElement>& rowBasis, const FieldElement& affineShift)const{
    return unique_ptr<BivariateExpansionProof>(new BivariateDummy(rowBasis));
}
    
unique_ptr<BivariateExpansionProof> Prover_Dummy::proofColumn_lowDegree(const BivariateExpansionProof& src, const size_t column_index ,const vector<FieldElement>& columnBasis, const FieldElement& affineShift)const{
    return unique_ptr<BivariateExpansionProof>(new BivariateDummy(columnBasis));
}

//}

} //namespace PCP_Project
} //namespace Prover
