///**
//*       @file  AffinePolynomialNTL_UTEST.cpp
//*      @brief  A file for unit testing the class- AffinePolynomialNTL_UTEST
//*      Unit tests use GTEST library
//*      This file can be linked to external main function
//*      (for example for wide system unit-testing)
//*
//*
//*	Additionally uses libraries NTL and GTEST.
//*
//*     @author  Ariel Gabizon, ariel.gabizon@gmail.com
//* =====================================================================================
//*/
//
//#include <iostream>
//#include "SubspacePolynomial.hpp"
//#include "AffinePolynomialNTL.hpp"
//#include <algebraLib/FiniteField.hpp>
//#include "AlgebraCommon.hpp"
//#include "AlgebraTestUtils.hpp"
//#include <gtest/gtest.h>
//
//
//
//using namespace Algebra;
//using NTL::GF2E;
//
//using std::vector;
//
//	namespace{
//		const int NUM_OF_ELEMENTS = 20;
//		const int NUM_OF_COEFF = 15;
//		const int NUM_OF_REPETITIONS = 20;
//		const int NUM_OF_SPEEDREPETITIONS = 100000;
//
//		/******************************************************************************************************/
//		/***************************************** Helper functions *******************************************/
//		/******************************************************************************************************/
//
//		/* The function evaluates an Affine polynomial at a given point and returns the result.
//		It is used in the tests here to check that AffinePolynomialNTL.eval is correct
//		*/
//		FieldElement evalNoMat(const FieldElement& x, const vector<FieldElement>& coefficients, const FieldElement& constantFactor)  {
//
//			unsigned long expo = 1;
//			FieldElement res = zero();
//
//			/*we take advantage of the fact squaring is very efficient in GF2E,
//			and that x^{2^i} = (x^{2^{i-1})^2. Each iteration we square the current power
//			of x- power_x, multiply it by the relevant coefficient, and add to the sum,
//			resulting in sum_{i=0}^{size-1} coeffcieints[i]*x^{2^i} */
//
//			FieldElement currPower, power_x = x;
//			for (unsigned long i = 0; i < coefficients.size(); i++) {
//				currPower = coefficients[i] * power_x;
//
//				res += currPower;
//				power_x = FieldElement::sqr(power_x);
//			}
//
//			return	res + constantFactor;
//		}
//
//
//
//
//		/******************************************************************************************************/
//		/***************************************** Tests *******************************************/
//		/******************************************************************************************************/
//
//		//Testing constructor of AffinePolynomial 
//		TEST(Algebra, AffinePolynomialNTLConstructor){
//			const FieldElement ZERO = zero();
//			const FieldElement ONE = one();
//
//			vector<FieldElement> coeff = randomFieldElementVector(NUM_OF_COEFF);
//			FieldElement constFact = randomFieldElement();
//			AffinePolynomialNTL p(coeff, constFact);
//			for (int i = 0; i < NUM_OF_COEFF; i++){
//				if (coeff[i] != ZERO)
//					EXPECT_EQ(coeff[i], p.getCoefficient(Infrastructure::POW2(i)));
//
//			}
//			EXPECT_EQ(constFact, p.getCoefficient(0));
//
//
//		}
//
//		//testing AffinePolynomialNTL evaluates correctly
//		TEST(Algebra, AffinePolynomialNTLEval){
//			const FieldElement ZERO = zero();
//			const FieldElement ONE = one();
//
//			vector<FieldElement> coeff = randomFieldElementVector(NUM_OF_COEFF);
//			FieldElement constFact = randomFieldElement();
//			AffinePolynomialNTL p(coeff, constFact);
//			FieldElement x;
//			for (int i = 0; i < NUM_OF_REPETITIONS; i++){
//				x = randomFieldElement();
//				EXPECT_EQ(evalNoMat(x, coeff, constFact), p.eval(x));
//
//			}
//
//		}
//
//
//#define CHECK_SPEED
//#ifdef CHECK_SPEED
//
//		//testing how fast AffinePolynomialNTL evaluates
//		TEST(Algebra, AffinePolynomialNTLSpeed){
//			const FieldElement ZERO = zero();
//			const FieldElement ONE = one();
//
//			vector<FieldElement> coeff = randomFieldElementVector(NUM_OF_COEFF);
//			FieldElement constFact = randomFieldElement();
//			AffinePolynomialNTL p(coeff, constFact);
//			FieldElement x, y;
//			for (int i = 0; i < NUM_OF_SPEEDREPETITIONS; i++){
//				x = randomFieldElement();
//				p.eval(x);
//
//			}
//
//		}
//
//#endif
//		//
//		//Tests the AffinePolynomial.eval function - mainly it's speed compared to evalNoMat
//		//TEST(Algebra, affineEvalSpeed){
//
//
//		//FieldElement zero = generateConstant<FieldElement, 0>();
//		//FieldElement one = generateConstant<FieldElement, 1>();
//		//FieldElement one2 = generateConstant<FieldElement, 1>();
//		//FieldElement coeff[1] = { one2 };
//		//FieldElement* coefficients = coeff;
//		//AffinePolynomial a;
//
//		//FieldElement x = mapIntegerToFieldElement(0, a.fieldDeg, 2);
//		//FieldElement y = FieldElement(NTL::GF2EInfo->random_GF2E());
//		////AffinePolynomial a(coefficients, 1, one);//TODO - this constructor seems problematic
//		//a.allocateSize(35);
//		//for (int i = 0; i < 35; i++){
//		//	a.setCoefficient(i, FieldElement(NTL::GF2EInfo->random_GF2E()));
//		//}
//		//
//
//		//y = FieldElement(NTL::GF2EInfo->random_GF2E());
//		//a.computeMat();
//		////a.setCoefficient(1, x);
//		//EXPECT_EQ(a.evalNoMat(y), a.eval(y));
//		//for (int i = 0; i < 200000; i++){
//		//	a.eval(y);
//		//}
//		//std::cout << "y:" << y << "a(y)" << a.eval(y);
//
//	}
//	//Tests the AffinePolynomial.eval function - mainly it's speed compared to evalNoMat
//	TEST(Algebra, affineEvalSpeed2){
//
//
//		//FieldElement zero = generateConstant<FieldElement, 0>();
//		//FieldElement one = generateConstant<FieldElement, 1>();
//		//FieldElement one2 = generateConstant<FieldElement, 1>();
//		//FieldElement coeff[1] = { one2 };
//		//FieldElement* coefficients = coeff;
//		//AffinePolynomial a;
//
//		//FieldElement x = mapIntegerToFieldElement(0, a.fieldDeg, 2);
//		//FieldElement y = FieldElement(NTL::GF2EInfo->random_GF2E());
//		////AffinePolynomial a(coefficients, 1, one);//TODO - this constructor seems problematic
//		//a.allocateSize(35);
//		//for (int i = 0; i < 35; i++){
//		//	a.setCoefficient(i, FieldElement(NTL::GF2EInfo->random_GF2E()));
//		//}
//
//
//		//y = FieldElement(NTL::GF2EInfo->random_GF2E());
//		//a.computeMat();
//		////a.setCoefficient(1, x);
//		//EXPECT_EQ(a.evalNoMat(y), a.eval(y));
//		//for (int i = 0; i < 200000; i++){
//		//	a.evalNoMat(mapIntegerToFieldElement(0, a.fieldDeg, i));
//		//}
//		//std::cout << "y:" << y << "a(y)" << a.eval(y);
//
//		//	}
//
//		////checks that the circuit containing the selector poly, now constructed with the new subspacepolynomial.eval, works well
//		//TEST(Algebra, selectPoly){
//		//	for (int i = 0; i < 7; i++){
//		//	a.(mapIntegerToFieldElement(0, a.fieldDeg, i));
//		//	}
//
//
//
//	} //anonimus namespace
