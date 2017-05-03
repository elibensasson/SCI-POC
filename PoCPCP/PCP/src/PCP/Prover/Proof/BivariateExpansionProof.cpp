#include "BivariateExpansionProof.hpp"
#include "common/Infrastructure/Infrastructure.hpp"

namespace PCP_Project{

using CryptoCommitment::hashDigest_t;
using CryptoCommitment::path_t;
using CryptoCommitment::constructMerkleTree;
using CryptoCommitment::logBytesPerHash;
using Infrastructure::Log2;
using Infrastructure::POW2;
using Algebra::FieldElement;

MerkleTree::MerkleTree(FieldElement const*const src, short src_log_size): logSizeBytes_(src_log_size + Log2(sizeof(FieldElement))){
    
    tree_ = new hashDigest_t[POW2(src_log_size)];
    root_ = constructMerkleTree(src,logSizeBytes_,tree_);
}

MerkleTree::MerkleTree(short src_log_size): logSizeBytes_(src_log_size+Log2(sizeof(FieldElement))){
    
    tree_ = new hashDigest_t[POW2(src_log_size)];
}

void MerkleTree::constructSubtree(FieldElement const*const src, const size_t sigment_logLen, const size_t sigment_index){
    constructMerkleSubTree(src,logSizeBytes_,sigment_logLen + Log2(sizeof(FieldElement)),sigment_index,tree_);
}

void MerkleTree::finishTreeAfterSegments(const size_t sigment_logLen){
    short currSrcLogLen = logSizeBytes_ - (sigment_logLen + Log2(sizeof(FieldElement)));
    root_ = constructMerkleTree(tree_ + POW2(currSrcLogLen), currSrcLogLen + logBytesPerHash, tree_);
}

const hashDigest_t& MerkleTree::getRoot()const{
    return root_;
}
short MerkleTreeInterface::getBlockSize(){
    return (1<<logBytesPerHash)/sizeof(FieldElement);
}

short MerkleTreeInterface::getDualBlockSize(){
    return 2*getBlockSize();
}

size_t MerkleTreeInterface::getBlockIndex(const size_t elementIndex){
    return elementIndex/getBlockSize();
}

short MerkleTreeInterface::getOffsetInDualBlock(const size_t index){
    return index - (getBlockIndex(index) & ~1UL)*getBlockSize();
}

path_t MerkleTree::getPathToBlock(const size_t blockIndex)const{
    return CryptoCommitment::getPathToBlock(tree_,logSizeBytes_,blockIndex);
}

bool MerkleTree::verifyPathToBlock(const path_t& path, FieldElement const*const blockData, const size_t blockIndex)const{
    return CryptoCommitment::verifyPathToBlock(blockData, root_, path, blockIndex);
}
MerkleTree::~MerkleTree(){
    delete[] tree_;
}
} //PCP_Project
