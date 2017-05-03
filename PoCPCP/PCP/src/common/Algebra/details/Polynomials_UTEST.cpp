//** testing mainly the integration of Matan's FFT in Polynomials.cpp**//
#include <gtest/gtest.h>
#include "NTL/GF2EX.h"
#include "NTL/GF2E.h"
#include "NTL/GF2X.h"

#include "Polynomials.hpp"
#include <NTL/GF2XFactoring.h>
#include "FFT/src/FFT.h"
#include "../AlgebraTestUtils.hpp"
#include "common/Algebra/FFT/LDERunner.hpp"
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
using namespace FFF;
using Infrastructure::POW2;
using Algebra::FieldElement;
using Algebra::zero;
using Algebra::one;

void generateBasis(FFF::Element* e, len_t l) {
	for (unsigned int i = 0; i < l; ++i) {
		Element::c_setZero(e + i);
		e[i].c[i >> Element::log_bits_in_cell] ^= (1ULL << (i & andMask(Element::log_bits_in_cell)));
	}
}

cell_t generateByte() {
	return rand() % (256);
}

cell_t generateCell() {
	cell_t res = 0;
	for (unsigned int i = 0; i < sizeof(cell_t); ++i)
		res ^= (generateByte() << (i * bits_in_byte));
	return res;
}

void generateElement(Element& e) {
	for (unsigned int i = 0; i < Element::element_len; ++i)
		e.c[i] = generateCell();
}

void generatePolynomial(Element* p, len_t p_len) {
	for (unsigned int i = 0; i < p_len; ++i)
		generateElement(p[i]);
}

TEST(Polynomials, MatanFftTest) {
	//FFF::Tests::testFFT(3, 5);
	int basis_size = 3;
	int poly_len = 5;
	Element* fft_basis = (Element*)alloca((sizeof(Element)*(basis_size + 1)));
	Element* fft_poly = Malloc(Element, poly_len); // TODO: AlgFFT frees the allocated Element - this is so WRONG!!
	generateBasis(fft_basis, basis_size + 1);
	generatePolynomial(fft_poly, poly_len);
	Element shift = fft_basis[basis_size];
	FFF::Basis b(fft_basis, basis_size, shift);
	FFF::FFT f(b,FFF::FFT_OP);
	f.AlgFFT(&fft_poly, poly_len);
	free(fft_poly);
}
