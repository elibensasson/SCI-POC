#include <gadgetlib/protoboard.hpp>
#include <gadgetlib/common_use.hpp>
#include <algebraLib/variable.hpp>
#include <algebraLib/variable_operators.hpp>
#include <gtest/gtest.h>
#include <memory>

using namespace gadgetlib;
using namespace Algebra;

namespace{
	TEST(protoboard, boolianityCheck){
		initGF2EwithPrimitivePoly();
		ProtoboardPtr pb = Protoboard::create();
		Algebra::Variable x("x");
		pb->enforceBooleanity(x,Opcode::NONE);
		pb->val(x) = Algebra::zero();
		EXPECT_TRUE(pb->isSatisfied(Opcode::NONE));
		pb->val(x) = Algebra::one();
		EXPECT_TRUE(pb->isSatisfied(Opcode::NONE));
		pb->val(x) = Algebra::FElem(getGF2E_X());
		EXPECT_FALSE(pb->isSatisfied(Opcode::NONE));
	}
	TEST(protoboard, getUsedvariablesCheck){
		initGF2EwithPrimitivePoly();
		ProtoboardPtr pb = Protoboard::create();
		Algebra::Variable::set usedVariables = pb->getUsedVariables();
		EXPECT_TRUE(usedVariables.size() == 0);
		Algebra::Variable x("x");
		Algebra::Variable y("y");
		Algebra::Variable z("z");
		CircuitPolynomial p1(x+y);
		pb->addGeneralConstraint(p1,"x + y", Opcode::NONE);
		usedVariables = pb->getUsedVariables();
		EXPECT_TRUE(usedVariables.size() == 2);
		CircuitPolynomial p2(x + y + FElem(getGF2E_X()));
		pb->addGeneralConstraint(p2, "x + y + FElem(x)", Opcode::NONE);
		usedVariables = pb->getUsedVariables();
		EXPECT_TRUE(usedVariables.size() == 2);
		CircuitPolynomial p3(z + y);
		pb->addGeneralConstraint(p3, "z + y", Opcode::NONE);
		usedVariables = pb->getUsedVariables();
		EXPECT_TRUE(usedVariables.size() == 3);
		std::vector<CircuitPolynomial> input = { p1, p2, p3};
		CircuitPolynomial p4(PolynomialOperation::MUL, input);
		pb->addGeneralConstraint(p4, "(x + y) * (x + y + FElem(x)) * (z + y)", Opcode::NONE);
		usedVariables = pb->getUsedVariables();
		EXPECT_TRUE(usedVariables.size() == 3);
		pb->addGeneralConstraint(p4, "(x + y) * (x + y + FElem(x)) * (z + y)", Opcode::ADD);
		usedVariables = pb->getUsedVariables();
		EXPECT_TRUE(usedVariables.size() == 3);
	}

	TEST(protoboard, memory){
		initGF2EwithPrimitivePoly();
		ProtoboardPtr pb = Protoboard::create();
		const FElem g = Algebra::FElem(getGF2E_X());
		std::vector<FElem> elements(3);
		elements[0] = g;
		for (int i = 1; i < elements.size(); i++){
			elements[i] = elements[i - 1] * g;
		}
		EXPECT_ANY_THROW(pb->loadValue(g));
		pb->storeValue(g, elements[1]);
		EXPECT_EQ(pb->loadValue(g), elements[1]);
		pb->storeValue(g, elements[2]);
		EXPECT_EQ(pb->loadValue(g), elements[2]);
		EXPECT_NE(pb->loadValue(g), elements[1]);

	}
}
