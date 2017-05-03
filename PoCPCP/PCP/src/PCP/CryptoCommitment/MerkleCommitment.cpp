#include "MerkleCommitment.hpp"
#include "common/Infrastructure/Infrastructure.hpp"
#include <omp.h>
#include <string>
//#include <sys/mman.h>

namespace PCP_Project{
namespace CryptoCommitment{

//hashes 64 bytes from src into 32 bytes in dst
void hash(void const* const src, void * const dst){
    
    static const short hash_src_len = 2*sizeof(hashDigest_t);
    
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256,src,hash_src_len);
    SHA256_Final((unsigned char*)dst,&sha256);
}

hashDigest_t hash(void const* const src){
    hashDigest_t res;
    hash(src,&res);
    return res;
}

bool operator==(const hashDigest_t& a, const hashDigest_t& b){
    return 0 == std::memcmp(&a,&b,sizeof(hashDigest_t));
}

bool operator!=(const hashDigest_t& a, const hashDigest_t& b){
    return !(a==b);
}

bool operator<(const hashDigest_t& a, const hashDigest_t& b){
    return 0 < std::memcmp(&a,&b,sizeof(hashDigest_t));
}



//
// Constructs a Merkle tree for the src buffer (srcLen expected in bytes)
// The tree is written to dst, and its root is returned.
// It is expected src_logLen is in bytes.
// It is expected the size of dst in bytes is at least srcLen.
//
hashDigest_t constructMerkleTree(void const* const src, const short src_logLen, void * const dst){

    using Infrastructure::POW2;

    short curr_dst_logLen = src_logLen - logBytesPerHash - 1;
    hashDigest_t* curr_src = (hashDigest_t*)src;
    _COMMON_ASSERT(src_logLen > logBytesPerHash,"It is assumed the source length contains a power of 2 blocks of 256 bits, src_logLen = " + std::to_string(src_logLen));

    while (curr_dst_logLen >= 0){
        hashDigest_t* curr_dst = ((hashDigest_t*)dst) + (1UL<<curr_dst_logLen);
        const size_t curr_dst_len = POW2(curr_dst_logLen);
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
        for(long long i=0; i<curr_dst_len; i++){
            hash(curr_src+(i<<1UL),curr_dst+i);
        }
        //msync(curr_src,curr_dst_len<<(1+logBytesPerHash),MS_ASYNC);

        curr_src = curr_dst;
        curr_dst_logLen--;
    }

        //return tree root
        return ((hashDigest_t*)dst)[1];
}

//
// Constructs a Merkle sub-tree for a sigment in the src buffer (srcLen expected in bytes)
// The sub - tree is written to dst
// It is expected src_logLen is in bytes.
// It is expected the size of dst in bytes is at least srcLen.
//
void constructMerkleSubTree(void const* const src, const short src_logLen, const size_t sigment_logLen, const size_t sigment_index, void * const dst){

    using Infrastructure::POW2;

    short curr_dst_row_logLen = src_logLen - logBytesPerHash - 1;
    short curr_sigment_logLen = sigment_logLen - logBytesPerHash - 1;
    hashDigest_t* curr_src_shifted = (hashDigest_t*)src + POW2(sigment_logLen - logBytesPerHash)*sigment_index;
    
    _COMMON_ASSERT(sigment_logLen > logBytesPerHash,"It is assumed the source length contains a power of 2 blocks of 256 bits");

    while (curr_sigment_logLen >=0){
        hashDigest_t* curr_dst_row = ((hashDigest_t*)dst) + POW2(curr_dst_row_logLen);
        hashDigest_t* curr_dst_shifted = curr_dst_row + POW2(curr_sigment_logLen)*sigment_index;
        const size_t curr_sigment_len = POW2(curr_sigment_logLen);
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
        for(long long i=0; i<curr_sigment_len; i++){
            hash(curr_src_shifted+(i<<1UL),curr_dst_shifted+i);
        }
        //msync(curr_src_shifted,curr_sigment_len<<(1+logBytesPerHash),MS_ASYNC);

        curr_src_shifted = curr_dst_shifted;
        curr_dst_row_logLen--;
        curr_sigment_logLen--;
    }
}

path_t getPathToBlock(void const*const tree, const short src_logLen, const size_t blockIndex){
    path_t result;
    hashDigest_t const*const tree_ = (hashDigest_t const*const)tree;
    const short firstRow_logLen = src_logLen - logBytesPerHash - 1;
    size_t curr_offset = (1UL<<firstRow_logLen) + (blockIndex>>1);
    
    while(curr_offset > 1UL){
        result.push_back(tree_[curr_offset ^ 1UL]);
        curr_offset >>= 1;
    }
    
    return result;
}

bool verifyPathToBlock(void const*const blockData, const hashDigest_t& root, const path_t& path, const size_t blockIndex){
    const short firstRow_logLen = path.size();
    size_t curr_offset = (1UL<<firstRow_logLen) + (blockIndex>>1);

    auto currHash = hash(blockData);

    for(size_t i=0; i < path.size(); i++){
        hashDigest_t hash_src[2];
        hash_src[(curr_offset&1UL) ^ 1UL] = path[i];
        hash_src[(curr_offset&1UL)] = currHash;

        currHash = hash(&hash_src);

        curr_offset >>= 1;
    }
    
    //verify root
    return currHash == root;
}

} //namespace CryptoCommitment
} //namespace PCP_Project
