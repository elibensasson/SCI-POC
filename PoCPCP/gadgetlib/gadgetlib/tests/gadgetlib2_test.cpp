#include <NTL/GF2XFactoring.h>
#include <NTL/GF2E.h>
#include <gtest/gtest.h>

#define IRR_DEGREE 63

void initNTL() {
    ::NTL::GF2X Irr;
    ::NTL::BuildIrred(Irr, IRR_DEGREE);
    ::NTL::GF2E::init(Irr);
}

int main(int argc, char **argv) {
    initNTL();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}