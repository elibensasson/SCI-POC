#include <gtest/gtest.h>
#include <algebraLib/FiniteField.hpp>

int main(int argc, char *argv[]) {
    
    Algebra::FiniteField(64).setContext();
	testing::InitGoogleTest(&argc, argv);
	size_t res = RUN_ALL_TESTS();
    Algebra::FiniteField::releaseContextField();
    return res;
}
