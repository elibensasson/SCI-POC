#ifndef MERKLECOMMITMENT_HPP__
#define MERKLECOMMITMENT_HPP__

#include <openssl/sha.h>
#include <vector>
#include <array>
#include <cstring>

namespace PCP_Project{
namespace CryptoCommitment{

typedef struct hashDigest_t{
    char buffer[SHA256_DIGEST_LENGTH];
};

bool operator==(const hashDigest_t& a, const hashDigest_t& b);
bool operator!=(const hashDigest_t& a, const hashDigest_t& b);
bool operator<(const hashDigest_t& a, const hashDigest_t& b);

typedef std::vector<hashDigest_t> path_t;

//hashes 64 bytes from src into 32 bytes in dst
void hash(void const* const src, void * const dst); 
hashDigest_t hash(void const* const src); 

const short logBytesPerHash = 5;

//
// Constructs a Merkle tree for the src buffer (srcLen expected in bytes)
// The tree is written to dst, and its root is returned.
// It is expected src_logLen is in bytes.
// It is expected the size of dst in bytes is at least srcLen.
//
hashDigest_t constructMerkleTree(void const* const src, const short src_logLen, void * const dst);

//
// Constructs a Merkle sub-tree for a sigment in the src buffer (srcLen expected in bytes)
// The sub - tree is written to dst
// It is expected src_logLen is in bytes.
// It is expected the size of dst in bytes is at least srcLen.
//
void constructMerkleSubTree(void const* const src, const short src_logLen, const size_t sigment_logLen, const size_t sigment_index, void * const dst);

path_t getPathToBlock(void const*const tree, const short src_logLen, const size_t blockIndex);

bool verifyPathToBlock(void const*const blockData, const hashDigest_t& root, const path_t& path, const size_t blockIndex);

} //namespace CryptoCommitment
} //namespace PCP_Project

#endif //#ifndef MERKLECOMMITMENT_HPP__
