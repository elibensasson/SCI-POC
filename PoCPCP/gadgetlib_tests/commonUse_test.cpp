/** @file
*****************************************************************************
Unit tests for gadgetlib2 - main() for running all tests
*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#include <gtest/gtest.h>
#include <gadgetlib/common_use.hpp>
#include <algebraLib/CircuitPolynomial.hpp>

using namespace gadgetlib;
using namespace Algebra;

namespace{
	TEST(commmonUse, getSelectorTEST){
		int programLength = 20;
		Algebra::UnpackedWord wannaBePc(3,"wannaBePc");
		CircuitPolynomial selector;
		EXPECT_ANY_THROW(selector = getSelector(programLength, 2, wannaBePc));
		wannaBePc = UnpackedWord(10, "wannaBePc");
		selector = getSelector(programLength, 3, wannaBePc);
		VariableAssignment assignment;
		for (int i = 0; i < wannaBePc.size(); ++i){
			assignment.insert(std::pair<Variable, FElem>(wannaBePc[i], Algebra::zero()));
		}
		assignment[wannaBePc[0]] = Algebra::one();
		assignment[wannaBePc[1]] = Algebra::one();
		EXPECT_EQ(selector.eval(assignment), Algebra::one());
		EXPECT_FALSE(selector.isSatisfied(assignment));
		assignment[wannaBePc[0]] = Algebra::zero();
		EXPECT_TRUE(selector.isSatisfied(assignment));
		assignment[wannaBePc[0]] = Algebra::one();
		assignment[wannaBePc[2]] = Algebra::one(); 
		EXPECT_TRUE(selector.isSatisfied(assignment));
	}

} // namespace