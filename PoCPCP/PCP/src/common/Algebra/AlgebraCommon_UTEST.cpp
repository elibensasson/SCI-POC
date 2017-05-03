
#include <gtest/gtest.h>
#include "AlgebraCommon.hpp"
#include "AlgebraTestUtils.hpp"


/* The purpose of this class is to test the methods from the AlgebraCommon file.
See explanation there about the methods.
*/


using namespace Algebra;
using Algebra::details::Basis;
using Algebra::FieldElement;




#define ARRAY_SIZE 10
//Testing conversion between FieldElement Array and FieldElement vector
TEST(Algebra, FieldElementArrayToElementVector){
	FieldElement c[ARRAY_SIZE];
	randomFieldElementArray(ARRAY_SIZE, c);
	vector<FieldElement> v = fieldPointArrayToElementVector(c,ARRAY_SIZE);
	for (int i = 0; i < ARRAY_SIZE; i++){
		EXPECT_EQ(c[i], (v[i]));

	}

}
