#include "MerkleCommitment.hpp"
#include <gtest/gtest.h>
#include <vector>
#include <stdint.h>

namespace{

using PCP_Project::CryptoCommitment::hashDigest_t;
using PCP_Project::CryptoCommitment::path_t;
using PCP_Project::CryptoCommitment::constructMerkleTree;
using PCP_Project::CryptoCommitment::getPathToBlock;
using PCP_Project::CryptoCommitment::verifyPathToBlock;
using PCP_Project::CryptoCommitment::hash;
using PCP_Project::CryptoCommitment::logBytesPerHash;
using std::vector;

void printHash(void const*const src){
    std::cout<<((int64_t*)src)[0]<<",";
    std::cout<<((int64_t*)src)[1]<<",";
    std::cout<<((int64_t*)src)[2]<<",";
    std::cout<<((int64_t*)src)[3];
}

TEST(MerkleCommitment,verifyTree){
    const short NUM_BLOCKS_LOG = 15;
    const short NUM_BYTES_LOG = NUM_BLOCKS_LOG+logBytesPerHash;
    vector<unsigned char> data(1UL<<NUM_BYTES_LOG);
    vector<hashDigest_t> merkleTree(1UL<<NUM_BLOCKS_LOG);

    //fill data
    for(size_t i=0; i<data.size(); i++){
        data[i] = (i<2? 1 : data[i-1] + data[i-2]);
    }

    //construct Merkle commitment
    const hashDigest_t root = constructMerkleTree(&data[0], NUM_BYTES_LOG, &merkleTree[0]);

    //verify tree
    for(size_t i=1; i < merkleTree.size()/2; i++){
        EXPECT_EQ(merkleTree[i],hash(&merkleTree[2*i]));
    }
}

TEST(MerkleCommitment,verifySubTree){
    const short NUM_BLOCKS_LOG = 15;
    const short NUM_BYTES_LOG = NUM_BLOCKS_LOG+logBytesPerHash;
    const short NUM_BLOCKS_PER_SIGMENT_LOG = NUM_BLOCKS_LOG/2;
    const short NUM_BYTES_PER_SIGMENT_LOG = NUM_BLOCKS_PER_SIGMENT_LOG+logBytesPerHash;
    
    //construct data
    vector<unsigned char> data(1UL<<NUM_BYTES_LOG);
    for(size_t i=0; i<data.size(); i++){
        data[i] = (i<2? 1 : data[i-1] + data[i-2]);
    }
    
    //construct referance
    vector<hashDigest_t> merkleTree_ref(1UL<<NUM_BLOCKS_LOG);
    const hashDigest_t root_ref = constructMerkleTree(&data[0], NUM_BYTES_LOG, &merkleTree_ref[0]);

    //construct tree using segment subtrees
    vector<hashDigest_t> merkleTree(1UL<<NUM_BLOCKS_LOG);
    for(size_t i=0; i< (1UL<<(NUM_BYTES_LOG - NUM_BYTES_PER_SIGMENT_LOG)); i++){
        constructMerkleSubTree(&data[0],NUM_BYTES_LOG,NUM_BYTES_PER_SIGMENT_LOG,i,&merkleTree[0]);
    }
    const short left_src_log_len = NUM_BYTES_LOG - NUM_BYTES_PER_SIGMENT_LOG;
    const hashDigest_t root = constructMerkleTree(&merkleTree[1UL<<left_src_log_len], left_src_log_len+logBytesPerHash, &merkleTree[0]);

    //verify results
    for(size_t i=0; i< merkleTree.size(); i++){
        EXPECT_EQ(merkleTree[i],merkleTree_ref[i]);
    }

}

TEST(MerkleCommitment,verifyPaths){
    const short NUM_BLOCKS_LOG = 6;
    const short NUM_BYTES_LOG = NUM_BLOCKS_LOG+logBytesPerHash;
    vector<unsigned char> data(1UL<<NUM_BYTES_LOG);
    vector<unsigned char> merkleTree(data.size());

    //fill data
    for(size_t i=0; i<data.size(); i++){
        data[i] = (i<2? 1 : data[i-1] + data[i-2]);
    }

    //construct Merkle commitment
    const hashDigest_t root = constructMerkleTree(&data[0], NUM_BYTES_LOG, &merkleTree[0]);

    //verify all paths
    for(size_t i=0; i< (1UL<<NUM_BLOCKS_LOG); i++){
        hashDigest_t const*const currData = ((hashDigest_t*)&data[0])+(i & ~1UL);
        const path_t path = getPathToBlock(&merkleTree[0],NUM_BYTES_LOG,i);
        
        EXPECT_TRUE(verifyPathToBlock(currData,root,path,i));
    }
}

}
