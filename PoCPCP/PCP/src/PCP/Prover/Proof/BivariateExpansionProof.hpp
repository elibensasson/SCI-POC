#ifndef BIVARIATEEXPANSIONPROOF_HPP__
#define BIVARIATEEXPANSIONPROOF_HPP__

#include <algebraLib/FieldElement.hpp>
#include <vector>
#include "PCP/CryptoCommitment/MerkleCommitment.hpp"

namespace PCP_Project{

class MerkleTreeInterface{
public:
    static short getBlockSize();
    static short getDualBlockSize();
    static size_t getBlockIndex(const size_t elementIndex);
    static short getOffsetInDualBlock(const size_t index);
    virtual const CryptoCommitment::hashDigest_t& getRoot()const = 0;
    virtual CryptoCommitment::path_t getPathToBlock(const size_t blockIndex)const = 0;
};

class MerkleTree : public MerkleTreeInterface{
public:
    MerkleTree(Algebra::FieldElement const*const src, short src_log_size);
    MerkleTree(short src_log_size);
    void constructSubtree(Algebra::FieldElement const*const src, const size_t sigment_logLen, const size_t sigment_index);
    void finishTreeAfterSegments(const size_t sigment_logLen);
    const CryptoCommitment::hashDigest_t& getRoot()const;
    CryptoCommitment::path_t getPathToBlock(const size_t blockIndex)const;
    bool verifyPathToBlock(const CryptoCommitment::path_t& path, Algebra::FieldElement const*const blockData, const size_t blockIndex)const;
    ~MerkleTree();
    CryptoCommitment::hashDigest_t* tree_;
private:
    short logSizeBytes_;
    CryptoCommitment::hashDigest_t root_;
};

class BivariateExpansionProof{
public:
    //size in FieldElements
    virtual size_t size()const = 0;
    virtual size_t getValIndex(const size_t rowId, const size_t colId)const = 0;
    virtual size_t getIndexRow(const size_t valIndex)const = 0;
    virtual size_t getIndexColumn(const size_t valIndex)const = 0;
    Algebra::FieldElement readValue(const size_t rowId, const size_t colId)const{
        return readValue(getValIndex(rowId,colId));
    }
    virtual Algebra::FieldElement readValue(const size_t idx)const{
        return getProof()[idx];
    }
    virtual std::vector<Algebra::FieldElement> getRow(const size_t rowId)const = 0;

    virtual const MerkleTreeInterface& merkleTree()const = 0;
    
    virtual ~BivariateExpansionProof(){};

protected:
    virtual Algebra::FieldElement const*const getProof()const = 0;
};

} //PCP_Project

#endif //#ifndef BIVARIATEEXPANSIONPROOF_HPP__
