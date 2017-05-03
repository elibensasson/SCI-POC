/** @file
*****************************************************************************
Unit tests for gadgetlib2 - main() for running all tests
*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#include <gtest/gtest.h>

#ifndef REMOVE_LD2_CODE

#include <NTL/GF2XFactoring.h>
#include <NTL/GF2E.h>
#define GADGETLIB_TESTS_IRR_DEGREE 63

void initNTL() {
	::NTL::GF2X Irr;
	::NTL::BuildIrred(Irr, GADGETLIB_TESTS_IRR_DEGREE);
	::NTL::GF2E::init(Irr);
}

#endif // REMOVE_LD2_CODE

int main(int argc, char **argv) {
#   ifndef REMOVE_LD2_CODE
	initNTL();
#   endif // REMOVE_LD2_CODE
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}