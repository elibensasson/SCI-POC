/*
 * Tests.h
 *
 *  Created on: Jul 2, 2014
 *      Author: matan
 */

#ifndef TESTS_H_
#define TESTS_H_
//#ifndef _MSC_VER //Ariel:not recognizing these on windows
#include "Definitions.h"
#include "Basis.h"
#include "Element.h"
#include "Chunk.cuh"
#include "Polynomials.h"
#include "FFT.h"
//#else
//#include "Definitions.h"
//#include "Basis.h"
//#include "Element.h"
//#include "FFT.h"
//#include "Polynomials.h"
//#include ""
//#endif
#ifdef __GPU
#include "FFT/src/GPU_FFT.cuh"
#endif
//#define __TEST_CORRECTNESS
namespace FFF{

class Tests {
public:
	Tests();
	virtual ~Tests();

	/**** Auxilliary Functions ****/
	static void generateBasis(Element* e, len_t l);
	static cell_t generateByte();
	static cell_t generateCell();
	static void generateElement(Element& e);
	static void generatePolynomial(Element* p, len_t p_len);
	static void generateElementsChunkArray(Elements_Chunk* a, len_t l);
	static void generateRandomElementsChunk(Elements_Chunk* t);
	static void generateChunkArray(Chunk* a,len_t l);
	static void generateRandomChunk(Chunk* t);
	static void extractElementFromChunk(idx_t i, const Chunk& c, Element&  e);
#ifdef __GPU
	static void generateChunkPolynomial(Chunk** p, len_t l);
	static Chunk* copyToGPU(Chunk* a, len_t l);
	static Chunk* copyFromGPU(Chunk* d_a,len_t l);
	static FFT* generateFFT(len_t dim);
#endif
template <typename T>
	static T generateRandom()
	{
		T ans;
		for(unsigned int i = 0 ; i< sizeof(T) ; ++i)
		{
				ans^=(((T)generateByte())<<(i*8));
		}
		return ans;

	}
/**** Tests ****/
	static bool test_generateSubspace(Basis& b);
	static bool test_cMul(unsigned int muls);
	static double measure_cMul();
	static bool testFFT(len_t basis_size, len_t poly_len);
	static bool testiFFT(len_t basis_size);
	static bool measureFFT(	unsigned int subspaceSize_low, unsigned int subapceSize_high,
							unsigned int* threads, len_t threads_len);
	static bool measureiFFT(	unsigned int subspaceSize_low, unsigned int subapceSize_high,
							unsigned int* threads, len_t threads_len);
	static bool chunkTestCellToChunk(Elements_Chunk* a, Chunk* c);
	static bool testCellToChunk(len_t chunks);
	static bool testChunkToCell(len_t chunks);
	static bool testMul(len_t chunks);
	static bool testMulByChunk(len_t chunks);
	static bool chunkTestMul(Chunk* c, Element* e, Chunk* c_copy);
	static bool chunkTestMulByChunk(Chunk* c, Chunk& b, Chunk* c_copy);
#ifdef __GPU
    static bool TestWFromUV_GPU(int dim, int WDim);
	static bool testGPUMultiExponentiate(int dim, int len_multiExp);
	static bool testmultiExpMult(int dim, int len_multiExp);
	static bool testTaylorExpansion_GPU(int dim, int taylor_dim);
	static bool testPartition_GPU(int dim, int partitionDim);
	static bool testLinearEvaluation_GPU(int dim);
#endif

};
} /* namespace FFF */
#endif /* TESTS_H_ */
