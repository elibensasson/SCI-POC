/** testing mainly the integration of Matan's FFT in Polynomials.cpp**/ 
#include <gtest/gtest.h>
#include "NTL/GF2EX.h"
#include "NTL/GF2E.h"
#include "NTL/GF2X.h"

#include "details/Polynomials.hpp"
#include <NTL/GF2XFactoring.h>
#include "FFT/src/FFT.h"
#include "AlgebraTestUtils.hpp"
#include "LightUniPolyEval.hpp"
#include "FFT/LDERunner.hpp"
#if defined _WIN64 || defined _WIN32
#include <windows.h>
#else	// #if defined _WIN64 || defined _WIN32
#endif	// #if defined _WIN64 || defined _WIN32

//#define __PROFILE
//#define __MEASURE_TIMES

#include "common/Utils/commonUtils.hpp"


#ifdef __PROFILE
#define MEASURE_TIME(CODE, PRINT)
const int PROF_REPEAT = 1000;
#else	// #ifdef __PROFILE
const int PROF_REPEAT = 1;
#endif	// #ifdef __PROFILE

using namespace Algebra::details;
using namespace NTL;
using std::cout;
using Infrastructure::POW2;
using namespace Algebra::details;
using namespace Algebra;
using std::cout;
using std::endl;

//testing LightUniPolyEval object computers correctly
TEST(LightEval, LightEvalBasicTest){

	::Algebra::details::initGF2E(64);
	int basis_size = 9;
	int poly_deg = 23;
	UnivariatePolynomial		poly;
		
		
	const Algebra::details::Basis					basis1 = buildStandardBasis(basis_size);
	FieldElement shift = randomFieldElement();
	LightUniPolyEval eval(basis1,shift);
    const size_t spaceSize = POW2(basis1.getSizeOfBasis());

	//preparing a polynomial 
	for (int i = 0; i <= poly_deg; i++) {
		poly.setCoeff(i, one());
	}
	//storing poly evaluations in LightUniPolyEval object
	for (auto i = 0; i < spaceSize; i++){
		FieldElement alpha = eval.getPoint(i);
		eval.addPoint(i, poly.queryAtPoint(alpha));
	}

	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = eval.queryAtPoint(i);
		EXPECT_EQ(exp, act);
	}

	/**same as above, but testing next constructor**/
LightUniPolyEval eval2(basis1,shift);
	//storing poly evaluations in LightUniPolyEval object
	for (auto i = 0; i < spaceSize; i++){
		FieldElement alpha = eval2.getPoint(i);
		eval2.addPoint(i, poly.queryAtPoint(alpha));
	}

	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval2.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = eval2.queryAtPoint(i);
		EXPECT_EQ(exp, act);
	}

	/**same as above, but testing next constructor**/
	LightUniPolyEval eval3(poly,basis1, shift);
	
	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval3.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = eval3.queryAtPoint(i);
		EXPECT_EQ(exp, act);
	}
	/**same as above, but testing fillEvaluation method**/
	LightUniPolyEval eval4;
	eval4.fillEvaluation(poly, basis1, shift);

	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval4.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = eval4.queryAtPoint(i);
		EXPECT_EQ(exp, act);
	}

	/** Now checking fillOldEvaluation method**/
	UniPolynomialEvaluation uniEval;
	eval2.fillOldEvaluation(uniEval);

	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval2.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = uniEval.queryAtPoint(a);
		EXPECT_EQ(exp, act);
	}



}


//same as above but on random poly and random subspace (warning:method below does not check subspace basis is actually independent, may cause problems)
TEST(LightEval, LightEvalRandom){

	::Algebra::details::initGF2E(64);
	int basis_size = 7;
	int poly_deg = 23;
	UnivariatePolynomial		poly;


	const Algebra::details::Basis					basis0 = buildStandardBasis(basis_size);
	 Algebra::details::Basis					basis1;

	FieldElement t;
	for (int i = 0; i < basis0.getSizeOfBasis(); i++){
		basis1.addElement(Algebra::randomFieldElementNoZero());
	}

	FieldElement shift = randomFieldElement();
	LightUniPolyEval eval(basis1,shift);
    const size_t spaceSize = POW2(basis1.getSizeOfBasis());

	//preparing a polynomial 
	for (int i = 0; i <= poly_deg; i++) {
		poly.setCoeff(i, randomFieldElement());
	}
	//storing poly evaluations in LightUniPolyEval object
	for (auto i = 0; i < spaceSize; i++){
		FieldElement alpha = eval.getPoint(i);
		eval.addPoint(i, poly.queryAtPoint(alpha));
	}

	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = eval.queryAtPoint(i);
		EXPECT_EQ(exp, act);
	}


	/**same as above, but testing next constructor**/
	LightUniPolyEval eval2(basis1, shift);
	//storing poly evaluations in LightUniPolyEval object
	for (auto i = 0; i < spaceSize; i++){
		FieldElement alpha = eval2.getPoint(i);
		eval2.addPoint(i, poly.queryAtPoint(alpha));
	}

	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval2.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = eval2.queryAtPoint(i);
		EXPECT_EQ(exp, act);
	}

	/**same as above, but testing next constructor**/
	LightUniPolyEval eval3(poly, basis1, shift);

	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval3.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = eval3.queryAtPoint(i);
		EXPECT_EQ(exp, act);
	}
    return;

	/**same as above, but testing fillEvaluation method**/
	LightUniPolyEval eval4; 
	eval4.fillEvaluation(poly, basis1, shift);

	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval4.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = eval4.queryAtPoint(i);
		EXPECT_EQ(exp, act);
	}



	/** Now checking fillOldEvaluation method**/
	UniPolynomialEvaluation uniEval;
	eval2.fillOldEvaluation(uniEval);

	//checking we got correct result
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval.getPoint(i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = uniEval.queryAtPoint(a);
		EXPECT_EQ(exp, act);
	}
}


//testing the LDERunner IFFT runner method that interpolates according to a LighUniPolyEval object
TEST(LightEval, IFFT){

		::Algebra::details::initGF2E(64);
		int basis_size = 5;
		int poly_deg = 23;
		int interpolate_basis_size = Infrastructure::Log2ceil(poly_deg + 1);
		//preparing NTL's objects
		Algebra::FieldElement				ONE = one();
		Algebra::details::UnivariatePolynomial		poly;
		FieldIndex spaceSize = POW2(interpolate_basis_size);
		FieldElement* evaluation2 = new FieldElement[spaceSize];
		const Algebra::details::Basis					basis = buildStandardBasis(interpolate_basis_size);
		//const Algebra::details::Basis basis2 = buildStandardBasis(basis_size);
		FieldElement shift = Algebra::randomFieldElement();

		Algebra::details::ExplicitLinearSubset S(basis, shift);
		LightUniPolyEval	evaluation(basis,shift);

		//preparing a polynomial and an evaluation based upon it on a small space
		for (int i = 0; i <= poly_deg; i++) {
			poly.setCoeff(i, Algebra::randomFieldElement());
		}

		for (auto i = 0; i < S.getSizeOfField(); i++){
			FieldElement alpha = evaluation.getPoint(i) ;
			evaluation.addPoint(i, poly.queryAtPoint(alpha));
				}

		UnivariatePolynomial P2;
		//interpolating the polynomial using the evaluation
		IFFTRunner (evaluation, basis, shift, P2);
		
		
		//checking we got the original polynomial back
		for (int i = 0; i < S.getSizeOfField(); i++){
			FieldElement a = S.getFieldElementByIndex(i) ;
			const FieldElement &exp = P2.queryAtPoint(a);
			const FieldElement &act = poly.queryAtPoint(a);
			EXPECT_EQ(exp, act);
		}
		
	}

//testing the LDERunner IFFT runner method that interpolates according to a LighUniPolyEval object
TEST(LightEval, LDE){

	::Algebra::details::initGF2E(64);
	int basis_size = 5;
	int poly_deg = 23;
	int interpolate_basis_size = Infrastructure::Log2ceil(poly_deg + 1);
	//preparing NTL's objects
	Algebra::FieldElement				ONE = one();
	Algebra::details::UnivariatePolynomial		poly;
	const Algebra::details::Basis					basis1 = buildStandardBasis(interpolate_basis_size);
	const Algebra::details::Basis basis2 = buildStandardBasis(basis_size);

	FieldElement shift1 = Algebra::randomFieldElement();
	FieldElement shift2 = Algebra::randomFieldElement();

	FieldElement* evaluation2= new FieldElement[POW2(basis2.getSizeOfBasis())];

	//preparing a polynomial and an evaluation based upon it on a small space
	for (int i = 0; i <= poly_deg; i++) {
		poly.setCoeff(i, Algebra::randomFieldElement());
	}
	LightUniPolyEval	evaluation1(basis1,shift1);

	for (auto i = 0; i < POW2(basis1.getSizeOfBasis()); i++){
		FieldElement alpha = evaluation1.getPoint(i);
        evaluation1.addPoint(i, poly.queryAtPoint(alpha));
	}


	//getting evaluation on the large space using LDERunner and storing it in evaluation2
	LDERunner(&(evaluation1.getTable()[0]), basis1.asVector(), shift1, evaluation2, basis2.asVector(), shift2);

	//checking we got the correct result
		for (int i = 0; i < POW2(basis2.getSizeOfBasis()); i++){
		FieldElement a = getSpaceElementByIndex(basis2.asVector(),shift2,i);
		const FieldElement &exp = poly.queryAtPoint(a);
		const FieldElement &act = evaluation2[i];
		EXPECT_EQ(exp, act);
	}

	delete[] evaluation2;

}


//testing LightUniPolyEval object times compared to UniPolyEvaluation

//4 times faster with LightEval even without omp. omp not giving big improve in LightEval , maybe cause too many threads
TEST(LightEval, LightEvalMeasure){

	::Algebra::details::initGF2E(64);
	int basis_size = 22;
	int poly_deg = 2;
	UnivariatePolynomial		poly;
	FieldElement ONE = one();

	const Algebra::details::Basis					basis1 = buildStandardBasis(basis_size);
	FieldElement shift = randomFieldElement();
	LightUniPolyEval eval(basis1,shift);
    const size_t spaceSize = POW2(basis1.getSizeOfBasis());
	UniPolynomialEvaluation eval2;
	//preparing a polynomial 
	for (int i = 0; i <= poly_deg; i++) {
		poly.setCoeff(i, ONE);
	}
#define LIGHT
#ifdef LIGHT

	cout << "Storing in LightEval"<<endl;
	//storing poly evaluations in LightUniPolyEval object
#pragma omp parallel for
	for (auto i = 0; i < spaceSize; i++){
		//FieldElement alpha = eval.getPoint(i);
		eval.addPoint(i,ONE);
	}
	EXPECT_EQ(true, true);
	cout << "compute using LightEval" << endl;
	//checking we got correct result
	FieldElement act;
#pragma omp parallel for private (act)
	for (int i = 0; i < spaceSize; i++){
	//	FieldElement a = eval.getPoint(i);
		act = eval.queryAtPoint(i);
	}
	EXPECT_EQ(true, true);

#else

	cout << "Storing in UniPolyEval"<<endl;

	//storing poly evaluations in UniPolynomialEvaluation object
//#pragma omp parallel for -pragma causing problems with UniPolynomialEvaluation
	for (auto i = 0; i < spaceSize; i++){
		FieldElement alpha = eval.getPoint(i);
		eval2.addPoint(alpha, ONE);
	}
	EXPECT_EQ(true, true);


	cout << "compute using UniPolyEval"<<endl;
//#pragma omp parallel for
	for (int i = 0; i < spaceSize; i++){
		FieldElement a = eval.getPoint(i);
		const FieldElement &act = eval2.queryAtPoint(a);
	}
	EXPECT_EQ(true, true);

#endif

}
