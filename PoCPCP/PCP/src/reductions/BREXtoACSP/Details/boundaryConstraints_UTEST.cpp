#include "languages/BREX/BrexWitnessChecker_UTEST.hpp"
#include "../BREXtoACSP.hpp"
#include "languages/ACSP/ACSPWitnessChecker.hpp"

#include <gtest/gtest.h>

namespace {

using PCP_Project::CBREXtoACSP;
using PCP_Project::BREXWitness;
using PCP_Project::BREXInstance;
using PCP_Project::ACSPWitnessChecker;
using std::pair;

typedef pair<BREXInstance,BREXWitness> BrexPair;

/******************************************
 *     Input consistency tests
 ******************************************/
TEST(BREXtoACSP_boundaryConstraints,complitness){
	const BrexPair src = PCP_UTESTS::generate_valid_boundary();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = CBREXtoACSP::reduceWitness(src.first,src.second);
    
    EXPECT_TRUE(ACSPWitnessChecker::verify(*instance,*witness));
}

TEST(BREXtoACSP_boundaryConstraints,soundness){
	const BrexPair src = PCP_UTESTS::generate_invalid_boundary();
    
    const auto instance = CBREXtoACSP::reduceInstance(src.first);
    const auto witness = CBREXtoACSP::reduceWitness(src.first,src.second);
    
    EXPECT_FALSE(ACSPWitnessChecker::verify(*instance,*witness));
}

}
