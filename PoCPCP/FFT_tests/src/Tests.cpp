/*
 * Tests.cpp
 *
 *  Created on: Jul 2, 2014
 *      Author: matan
 */

#include "Tests.h"
#include <cstdlib>
#include <cstring>
#include <ctime>
#ifdef __GNUC__
#include <sys/time.h>
#endif	// #ifdef __GNUC__
#include <omp.h>
#include <iostream>
#include <cstdio>
namespace FFF {

void Tests::generateBasis(Element* e, len_t l) {
	for (unsigned int i = 0; i < l; ++i) {
		Element::c_setZero(e + i);
		e[i].c[i >> Element::log_bits_in_cell] ^= (1 << (i & andMask(Element::log_bits_in_cell)));
	}
}

cell_t Tests::generateByte() {
	return rand() % (256);
}

cell_t Tests::generateCell() {
	cell_t res = 0;
	for (unsigned int i = 0; i < sizeof(cell_t); ++i)
		res ^= (generateByte() << (i * bits_in_byte));
	return res;
}

void Tests::generateElement(Element& e) {
	for (unsigned int i = 0; i < Element::element_len; ++i)
		e.c[i] = generateCell();
}

void Tests::generatePolynomial(Element* p, len_t p_len) {
	for (unsigned int i = 0; i < p_len; ++i)
		generateElement(p[i]);
}

bool Tests::testFFT(len_t basis_size, len_t poly_len) {
	double begin,end;
	bool flag = true;
	Element* basis = Malloc(Element,basis_size+1) ;
	Element* poly = Malloc(Element,poly_len);
	generateBasis(basis, basis_size + 1);
	generatePolynomial(poly, poly_len);
	Element shift = basis[basis_size];
	Basis b(basis, basis_size, shift);
	FFT f(b);

#ifdef __TEST_CORRECTNESS
	Element* poly_eval = Malloc(Element,1<<basis_size);
	Element t;
	for(unsigned int i = 0; i < (1<<basis_size); ++i ) {
		Basis::getElement(b,t,i);
		Polynomials::evaluate(poly,poly_len,t,poly_eval[i]);
	}
#endif
#ifdef __MEASURE
#ifndef __GPU
	begin = omp_get_wtime();
//	clock_t s,e;
//	s = clock();
#else	// #ifndef __GPU
#endif	// #ifndef __GPU
#endif	// #ifdef __MEASURE
	f.AlgFFT(&poly, poly_len);
#ifdef __MEASURE
#ifndef __GPU
	end = omp_get_wtime();
	double elapsed_secs = double(end - begin);
	std::cout << elapsed_secs << std::endl;
#else	// #ifndef __GPU
#endif	// #ifndef __GPU
#endif	// #ifdef __MEASURE
#ifdef __TEST_CORRECTNESS
	for(unsigned int i = 0; i < (1<<basis_size); ++i) {
		if(!Element::equals(poly_eval[i],poly[i])) {
			flag= false;
			std::cout << i << std::endl;
			Element::printElement(poly_eval[i]);
			std::cout << std::endl;
			Element::printElement(poly[i]);
			std::cout << std::endl;
			break;
		}
	}
	free(poly_eval);
#endif
	free(basis);
	free(poly);
	return flag;
}

bool Tests::testiFFT(len_t basis_size) {
	bool flag = true;
	len_t poly_len = 1 << basis_size;
	Element* basis = Malloc(Element,basis_size+1)
	;
	Element* poly = Malloc(Element,poly_len)
	;
	generateBasis(basis, basis_size + 1);
	generatePolynomial(poly, poly_len);
	Element shift = basis[basis_size];
	Basis b(basis, basis_size, shift);
	FFT f(b);
#ifdef __TEST_CORRECTNESS
	Element* poly_eval = Malloc(Element,1<<basis_size);
	memcpy(poly_eval,poly,poly_len*sizeof(Element));
#endif
#ifndef GPU
	double begin = omp_get_wtime();
#else
#endif
	f.AlgIFFT(poly);
#ifndef GPU
	double end = omp_get_wtime();
	double elapsed_secs = double(end - begin);
	std::cout << elapsed_secs << std::endl;
#else
#endif
#ifdef __TEST_CORRECTNESS
	Element t;
	Element res;
	for(unsigned int i = 0; i < (1<<basis_size); ++i ) {
		Basis::getElement(b,t,i);
		Polynomials::evaluate(poly,poly_len,t,res);
		if(!Element::equals(res,poly_eval[i])) {
			flag = false;
			break;
		}
	}
	free(poly_eval);
	if(!flag) {
		for(unsigned int i = 0; i < (1<<basis_size); ++i ) {
			Basis::getElement(b,t,i);
			Polynomials::evaluate(poly,poly_len,t,res);
			Element::printElement(res);
			std::cout << std::endl;
		}
	}
#endif
	free(basis);
	free(poly);
	return flag;
}
bool testSingleElement() {
	Element a, b, c, d;
	Tests::generateElement(a);
	Tests::generateElement(b);
	Element::c_mul(a, b, c);
	Element::naiveMul(a.c, b.c, d.c);
	if (!Element::equals(c, d)) {
		Element::printElement(a);
		std::cout << std::endl;
		Element::printElement(b);
		std::cout << std::endl;
		Element::printElement(c);
		std::cout << std::endl;
		Element::printElement(d);
		return false;
	}
	return true;
}
bool Tests::test_cMul(unsigned int muls) {
	for (unsigned int i = 0; i < muls; ++i)
		if (!testSingleElement())
			return false;
	return true;
}
#ifdef __MEASURE
double Tests::measure_cMul() {
	Element a, b, c;
	Tests::generateElement(a);
	Tests::generateElement(b);
#ifdef __GNUC__
	timespec start;
	timespec end;
	clock_gettime(CLOCK_REALTIME, &start);
#endif	//#ifdef __GNUC__
	Element::c_mul(a, b, c);
#ifdef __GNUC__
	clock_gettime(CLOCK_REALTIME, &end);
	return end.tv_nsec - start.tv_nsec;
#else
	return 0;
#endif	//#ifndef __GNUC__
}
#endif	// #ifdef __MEASURE
bool Tests::measureFFT(unsigned int subspaceSize_low,
		unsigned int subspaceSize_high, unsigned int* threads,
		len_t threads_len) {
	for (unsigned int t = 0; t < threads_len; ++t) {
		omp_max_threads = threads[t];
		for (unsigned int size = subspaceSize_low; size <= subspaceSize_high;
				++size) {
			if (!Tests::testFFT(size, 1 << size)) {
				return false;
			}
		}
	}
	return true;
}
bool Tests::measureiFFT(unsigned int subspaceSize_low,
		unsigned int subspaceSize_high, unsigned int* threads,
		len_t threads_len) {
	for (unsigned int t = 0; t < threads_len; ++t) {
		omp_max_threads = threads[t];
		for (unsigned int size = subspaceSize_low; size <= subspaceSize_high;
				++size) {
			if (!Tests::testiFFT(size)) {
				return false;
			}
		}
	}
	return true;
}
bool Tests::test_generateSubspace(Basis& b) {
	len_t b_size = b.getSize();
	bool ans = true;
	Element t;
	Element* es = b.getBasis();
	Element* res = (Element*) malloc(sizeof(Element) * (1 << b_size));
	FFT::generateSubspaceElements_cpu(es, b.getShift(), b_size, res);
	for (idx_t i = 0; i < (1 << b_size); ++i) {
		Basis::getElement(b, t, i);
		if (!Element::equals(t, res[i])) {
			ans = false;
			break;
		}
	}
	free(res);
	return ans;

}
#ifdef __GPU
bool Tests::chunkTestCellToChunk(Elements_Chunk* a, Chunk* c) {
	unsigned int j;
	for (unsigned int i = 0; i < Chunk::elements_in_chunk; ++i) {
		for (unsigned int j = 0; j < Element::element_len; ++j) {
			cell_t t = 0;
			for (unsigned int l = 0; l < Element::bits_in_cell; ++l)
				t ^= (((cell_t) (((c->v[l + j * Element::bits_in_cell]) >> i)
						& 1)) << l);
			if (t != a->e[i].c[j])
				return false;
		}
	}
	return true;
}
void Tests::generateRandomElementsChunk(Elements_Chunk* t) {
	unsigned int j;
	for (unsigned int i = 0; i < Chunk::elements_in_chunk; ++i)
		generateElement(t->e[i]);
}
void Tests::generateElementsChunkArray(Elements_Chunk* a, len_t l) {
	for (unsigned int i = 0; i < l; ++i)
		Tests::generateRandomElementsChunk(&(a[i]));
}
void Tests::generateRandomChunk(Chunk* t) {
	for (unsigned int i = 0; i < Chunk::cells_in_chunk; ++i) {
		t->v[i] = generateRandom<chunk_cell_t>();
	}
}
void Tests::generateChunkArray(Chunk* a, len_t l) {
	for (unsigned int i = 0; i < l; ++i)
		Tests::generateRandomChunk(&(a[i]));
}
bool Tests::testCellToChunk(len_t chunks) {
	Elements_Chunk *a = (Elements_Chunk*) malloc(
			sizeof(Elements_Chunk) * chunks);
	Chunk *c = (Chunk*) malloc(sizeof(Chunk) * chunks);
	generateElementsChunkArray(a, chunks);
	Chunk::normalToChunk(a, c, chunks, true);
	for (unsigned int ck = 0; ck < chunks; ++ck)
		if (!chunkTestCellToChunk(&(a[ck]), &(c[ck])))
			return false;
	return true;
}
bool Tests::testChunkToCell(len_t chunks) {
	//Allocate memory
	Elements_Chunk *a = (Elements_Chunk*) malloc(
			sizeof(Elements_Chunk) * chunks);
	Chunk *c = (Chunk*) malloc(sizeof(Chunk) * chunks);
	generateChunkArray(c, chunks);
	Chunk::chunkToNormal(c, a, chunks);
	for (unsigned int ck = 0; ck < chunks; ++ck)
		if (!chunkTestCellToChunk(&(a[ck]), &(c[ck])))
			return false;
	return true;
}
/*
 * C - is the output of the original algorithm.
 * c_copy - is a copy of the original input to the algorithm.
 * b - is the chunks that c_copy was multiplied by.
 * Everything should work of multiplying each chunk of c_copy in b will give the
 * 	corresponding chunk in c.
 */
bool Tests::chunkTestMulByChunk(Chunk* c, Chunk& b, Chunk* c_copy) {
	Element e_t;
	Element e;
	Element e_p;
	for (unsigned int i = 0; i < Chunk::elements_in_chunk; ++i) {
		extractElementFromChunk(i, *c_copy, (e_t));
		extractElementFromChunk(i, *c, (e_p));
		extractElementFromChunk(i, b, e);
		Element::naiveMul(e_t, e, e_t);
		if (!Element::equals((e_t), (e_p))) {
			std::cout << i;
			std::cout << std::endl << std::endl << std::endl;
			Element::printElement(e);
			std::cout << std::endl << std::endl << std::endl;
			c_copy->print();
			std::cout << std::endl << std::endl << std::endl;
			Element::printElement((e_t));
			std::cout << std::endl << std::endl << std::endl;
			Element::printElement((e_p));
			return false;
		}
	}
	return true;
}
bool Tests::chunkTestMul(Chunk* c, Element* e, Chunk* c_copy) {
	Element e_t;
	Element e_p;
	for (unsigned int i = 0; i < Chunk::elements_in_chunk; ++i) {
		extractElementFromChunk(i, *c_copy, (e_t));
		extractElementFromChunk(i, *c, (e_p));
		Element::naiveMul(e_t, *e, e_t);
		if (!Element::equals((e_t), (e_p))) {
			std::cout << i;
			std::cout << std::endl << std::endl << std::endl;
			Element::printElement(*e);
			std::cout << std::endl << std::endl << std::endl;
			c_copy->print();
			std::cout << std::endl << std::endl << std::endl;
			Element::printElement((e_t));
			std::cout << std::endl << std::endl << std::endl;
			Element::printElement((e_p));
			return false;
		}
	}
	return true;
}
bool Tests::testMulByChunk(len_t chunks) {
	//Decalre Variables & Allocate Memory
	bool flag = true;
	Chunk cm;
	Chunk* c = (Chunk*) malloc(sizeof(Chunk) * chunks);
	Chunk* c_copy = (Chunk*) malloc(sizeof(Chunk) * chunks);
	//Initiazlie Variables
	generateChunkArray(c, chunks);
	for (unsigned int i = 0; i < Element::ord; ++i)
		c->v[i] = 0;
	c->v[0] = -1;
	generateRandomChunk(&cm);

	memcpy(c_copy, c, sizeof(Chunk) * chunks);
	//Test
	Chunk::chunk_mul(c, &cm, chunks, c);
#ifdef __TEST_CORRECTNESS
//		Inspect Results
	for(unsigned int i = 0; i < chunks; ++i)
	{
		if(!chunkTestMulByChunk(&(c[i]),cm,&(c_copy[i]))) {
			flag = false;
			break;
		}
	}
//		Return answer
#endif
	free(c);
	free(c_copy);
	return true;
}
bool Tests::testMul(len_t chunks) {
	//Decalre Variables & Allocate Memory
	Element* e = (Element*) malloc(sizeof(Element));
	Chunk* c = (Chunk*) malloc(sizeof(Chunk) * chunks);
	Chunk* c_copy = (Chunk*) malloc(sizeof(Chunk) * chunks);
	//Initiazlie Variables
	generateChunkArray(c, chunks);
	for (unsigned int i = 0; i < Element::ord; ++i)
		c->v[i] = 0;
	c->v[0] = -1;
	for (unsigned int i = 0; i < Element::element_len; ++i)
		e->c[i] = generateRandom<cell_t>();
//		e->c[0]=-1;
//		e->c[1]=0;
	memcpy(c_copy, c, sizeof(Chunk) * chunks);
	//Test
	Chunk::mul(c, e, chunks, c);
#ifdef __TEST_CORRECTNESS
//		Inspect Results
	for(unsigned int i = 0; i < chunks; ++i)
	{
		if(!chunkTestMul(&(c[i]),e,&(c_copy[i])))
		return false;
	}
//		Return answer
#endif
	free(e);
	free(c);
	free(c_copy);
	return true;
}
void Tests::extractElementFromChunk(idx_t i, const Chunk& c, Element& e) {
	unsigned int j;
	Element::c_setZero(e);
	for (unsigned int j = 0; j < Element::ord; ++j) {
		e.c[j >> Element::log_bits_in_cell] ^= (((cell_t) (((c.v[j]) >> i) & 1))
				<< (j & andMask(Element::log_bits_in_cell)));
	}
}

/*
 * Number of elements in the polynomial.
 */
void Tests::generateChunkPolynomial(Chunk** p, len_t l) {
	len_t p_l = sizeCiel(l, Chunk::elements_in_chunk);
	chunk_cell_t mask = (1 << (l & andMask(Chunk::log_elements_in_chunk))) - 1;
	*p = (Chunk*) malloc(sizeof(Chunk) * p_l);
	generateChunkArray(*p, p_l);
	if(mask!=0){
	for (unsigned int i = 0; i < Chunk::cells_in_chunk; ++i) {
		(*p)[p_l - 1].v[i] &= mask;
	}
	}
	return;
}
Chunk* Tests::copyToGPU(Chunk* a, len_t l) {
	Chunk* d_a;
	cudaMalloc(&d_a, sizeof(Chunk) * l);
	cudaMemcpy(d_a, a, sizeof(Chunk) * l, cudaMemcpyHostToDevice);
	return d_a;
}
Chunk* Tests::copyFromGPU(Chunk* d_a, len_t l) {
	Chunk* a = (Chunk*) malloc(sizeof(Chunk) * l);
	cudaMemcpy(a, d_a, sizeof(Chunk) * l, cudaMemcpyDeviceToHost);
	return a;
}

FFT* Tests::generateFFT(len_t dim) {
	Element* e = (Element*) malloc(sizeof(Element) * (dim + 1));
	generateBasis(e, dim + 1);
	Basis b(e, dim, e[dim]);
	FFT* ans = new FFT(b);
	free(e);
	return ans;
}

bool Tests::testmultiExpMult(int dim, int len_multiExp) {
	bool flag = true;
	Element* e_p = (Element*) malloc(
			sizeof(Element) * (1 << (dim + Chunk::log_elements_in_chunk)));
	Element* e_s = (Element*) malloc(
			sizeof(Element)
					* (1 << (len_multiExp + Chunk::log_elements_in_chunk)));
	Chunk* p;
	Chunk* s;
	generateChunkPolynomial(&p, 1 << (dim+Chunk::log_elements_in_chunk));
	generateChunkPolynomial(&s, 1 << (len_multiExp+Chunk::log_elements_in_chunk));
	Chunk::chunkToNormal(p, (Elements_Chunk*) e_p, 1 << dim, true);
	Chunk::chunkToNormal(s, (Elements_Chunk*) e_s, 1 << len_multiExp, true);
	for (unsigned int i = 0; i < 1 << (dim + Chunk::log_elements_in_chunk); i +=
			(1 << (len_multiExp + Chunk::log_elements_in_chunk))) {
		for (unsigned int j = 0;
				j < 1 << (len_multiExp + Chunk::log_elements_in_chunk); ++j) {
			Element::c_mul(e_p[i + j], e_s[j], e_p[i + j]);
		}
	}
	Chunk* d_p = copyToGPU(p, 1 << dim);
	Chunk* d_s = copyToGPU(s, 1 << len_multiExp);
	GPU_FFT::multiExp_mult(1 << dim, d_p, d_s, 1 << len_multiExp);
	Chunk * h_res = copyFromGPU(d_p,1<<dim);
	Element * res = (Element*)malloc(sizeof(Element)*(1<<(dim+Chunk::log_elements_in_chunk)));
	Chunk::chunkToNormal(h_res,(Elements_Chunk*)res,1<<dim,true);
	for(unsigned int i = 0 ; i < 1<<(dim+Chunk::log_elements_in_chunk) ; ++i){
		if(!Element::equals(res[i],e_p[i])){
			flag = false;
			break;
		}
	}
	free(p);
	free(s);
	free(e_p);
	free(e_s);
	free(h_res);
	free(res);
	cudaFree(d_p);
	cudaFree(d_s);
	return flag;
}
bool Tests::testGPUMultiExponentiate(int dim, int len_multiExp) {
	bool flag = true;
	FFT* fft = generateFFT(dim);
	Chunk* p;
	generateChunkPolynomial(&p, 1 << dim);
	len_t p_len = 1;
	len_t len_sub = 1;
	if (dim > Chunk::log_elements_in_chunk) {
		p_len = 1 << (dim - Chunk::log_elements_in_chunk);
	}
	if (len_multiExp > Chunk::log_elements_in_chunk) {
		len_sub = 1 << (len_multiExp - Chunk::log_elements_in_chunk);
	}
	Element* e = (Element*) malloc(
			sizeof(Element) * (p_len << Chunk::log_elements_in_chunk));
	Chunk::chunkToNormal(p, (Elements_Chunk*) e, p_len, true);
	for (unsigned int i = 0; i < (1 << len_multiExp); ++i) {
		for (unsigned int j = i; j < (1 << dim); j += (1 << len_multiExp)) {
			Element::c_mul(e[j], fft->exps[dim - len_multiExp][i], e[j]);
		}
	}
	Chunk* h_cpu = (Chunk*) malloc(sizeof(Chunk) * p_len);
	Chunk::normalToChunk((Elements_Chunk*) e, h_cpu, p_len, true);
	Chunk* d_p = copyToGPU(p, p_len);
	Chunk* d_sub;
	cudaMalloc(&d_sub, sizeof(Chunk) * len_sub);
	GPU_FFT::multiExponentiate_gpu(fft, d_p, p_len, len_multiExp, d_sub);
	Chunk* h_p = copyFromGPU(d_p, p_len);
	for (unsigned int i = 0; i < p_len; ++i) {
		for (unsigned int j = 0; j < Chunk::cells_in_chunk; ++j) {
			if (h_cpu[i].v[j] != h_p[i].v[j]) {
				flag = false;
			}
		}
	}
	cudaFree(d_sub);
	cudaFree(d_p);
	free(h_cpu);
	free(h_p);
	free(p);
	free(e);
	delete fft;
	return flag;
}

bool Tests::testTaylorExpansion_GPU(int dim, int taylorDim){
	bool flag = true;
	FFT* fft = generateFFT(dim);
	GPU_FFT::setUpConstantMemory(fft);
	Chunk* p;
	generateChunkPolynomial(&p,1<<dim  );
	len_t p_len = 1;
	len_t len_sub = 1;
	if (dim > Chunk::log_elements_in_chunk) {
		p_len = 1 << (dim - Chunk::log_elements_in_chunk);
	}
	if (taylorDim > Chunk::log_elements_in_chunk) {
		len_sub = 1 << (taylorDim - Chunk::log_elements_in_chunk);
	}

	/*
	 * Calculate on CPU
	 */
	Element* e = (Element*) malloc(
			sizeof(Element) * (p_len << Chunk::log_elements_in_chunk));
	Chunk::chunkToNormal(p, (Elements_Chunk*) e, p_len, true);
	for(unsigned int i =0  ; i < (1<<dim) ; i+=(1<<taylorDim)){
		Polynomials::taylorExpansion(e+i,taylorDim);
	}


	/*
	 * Change the calculated elements to Chunk representation
	 */
	Chunk* h_cpu = (Chunk*) malloc(sizeof(Chunk) * p_len);
	Chunk::normalToChunk((Elements_Chunk*) e, h_cpu, p_len, true);


	/*
	 * Calculate using GPU
	 */
	Chunk* d_p = copyToGPU(p, p_len);
	GPU_FFT::taylorExpansion_gpu(fft,d_p,p_len,taylorDim);
	Chunk* h_p = copyFromGPU(d_p, p_len);

	for (unsigned int i = 0; i < p_len; ++i) {
		for (unsigned int j = 0; j < Chunk::cells_in_chunk; ++j) {
			if (h_cpu[i].v[j] != h_p[i].v[j]) {
				flag = false;
			}
		}
	}
	cudaFree(d_p);
	free(h_cpu);
	free(h_p);
	free(p);
	free(e);
	delete fft;
	return flag;
}
bool Tests::testPartition_GPU(int dim, int partitionDim){
	bool flag = true;
	FFT* fft = generateFFT(dim);
	GPU_FFT::setUpConstantMemory(fft);
	Chunk* p;
	generateChunkPolynomial(&p,1<<dim  );
	len_t p_len = 1;
	len_t len_sub = 1;
	if (dim > Chunk::log_elements_in_chunk) {
		p_len = 1 << (dim - Chunk::log_elements_in_chunk);
	}
	if (partitionDim > Chunk::log_elements_in_chunk) {
		len_sub = 1 << (partitionDim - Chunk::log_elements_in_chunk);
	}

	/*
	 * Calculate on CPU
	 */
	Element* tmp = (Element*) malloc(
			sizeof(Element) * (p_len << Chunk::log_elements_in_chunk));
	Element* e = (Element*) malloc(
			sizeof(Element) * (p_len << Chunk::log_elements_in_chunk));
	Chunk::chunkToNormal(p, (Elements_Chunk*) tmp, p_len, true);
	for(unsigned int i =0  ; i < (1<<dim) ; i+=(1<<partitionDim)){
		fft->GPartition_cpu(tmp+i,partitionDim,e+i);
	}
	free(tmp);

	/*
	 * Change the calculated elements to Chunk representation
	 */
	Chunk* h_cpu = (Chunk*) malloc(sizeof(Chunk) * p_len);
	Chunk::normalToChunk((Elements_Chunk*) e, h_cpu, p_len, true);


	/*
	 * Calculate using GPU
	 */
	Chunk* d_p = copyToGPU(p, p_len);
	Chunk* d_sub;
	cudaMalloc(&d_sub, sizeof(Chunk) * p_len);
	Chunk* h_p;
	if(GPU_FFT::partition(d_p,d_sub,p_len,partitionDim)){
		h_p = copyFromGPU(d_sub,p_len);
	} else {
	    h_p = copyFromGPU(d_p, p_len);
	}

	for (unsigned int i = 0; i < p_len; ++i) {
		for (unsigned int j = 0; j < Chunk::cells_in_chunk; ++j) {
			if (h_cpu[i].v[j] != h_p[i].v[j]) {
				flag = false;
			}
		}
	}
	cudaFree(d_sub);
	cudaFree(d_p);
	free(h_cpu);
	free(h_p);
	free(p);
	free(e);
	delete fft;
	return flag;
}
bool Tests::TestWFromUV_GPU(int dim, int UDim){
	bool flag = true;
	FFT* fft = generateFFT(dim);
	GPU_FFT::setUpConstantMemory(fft);
	Chunk* p;
	generateChunkPolynomial(&p,1<<dim  );
	len_t p_len = 1;
	len_t len_sub = 1;
	if (dim > Chunk::log_elements_in_chunk) {
		p_len = 1 << (dim - Chunk::log_elements_in_chunk);
	}
	if (UDim > Chunk::log_elements_in_chunk) {
		len_sub = 1 << (UDim - Chunk::log_elements_in_chunk);
	}

	/*
	 * Calculate on CPU
	 */
	Element* e = (Element*) malloc(
			sizeof(Element) * (p_len << Chunk::log_elements_in_chunk));
	Chunk::chunkToNormal(p, (Elements_Chunk*) e, p_len, true);
	fft->WFromUV_cpu_serial(e,UDim+1,MAX(dim,Chunk::log_elements_in_chunk));

	/*
	 * Change the calculated elements to Chunk representation
	 */
	Chunk* h_cpu = (Chunk*) malloc(sizeof(Chunk) * p_len);
	Chunk::normalToChunk((Elements_Chunk*) e, h_cpu, p_len, true);


	/*
	 * Calculate using GPU
	 */
	Chunk* d_p = copyToGPU(p, p_len);
	Chunk* d_sub;
	cudaMalloc(&d_sub, sizeof(Chunk) * len_sub);
	cudaMemcpy(d_sub,fft->gpu_subspace[fft->basis.getSize()-UDim-1],sizeof(Chunk)*len_sub,cudaMemcpyHostToDevice);
	Chunk* h_p;
	GPU_FFT::WFromUV(d_p,p_len,d_sub,UDim);
	h_p = copyFromGPU(d_p,p_len);
	for (unsigned int i = 0; i < p_len; ++i) {
		for (unsigned int j = 0; j < Chunk::cells_in_chunk; ++j) {
			if (h_cpu[i].v[j] != h_p[i].v[j]) {
				flag = false;
			}
		}
	}
	cudaFree(d_sub);
	cudaFree(d_p);
	free(h_cpu);
	free(h_p);
	free(p);
	free(e);
	delete fft;
	return flag;
}
bool Tests::testLinearEvaluation_GPU(int dim){
	bool flag = true;
	FFT* fft = generateFFT(dim);
	GPU_FFT::setUpConstantMemory(fft);
	Chunk* p;
	generateChunkPolynomial(&p,1<<dim  );
	len_t p_len = 1;
	if (dim > Chunk::log_elements_in_chunk) {
		p_len = 1 << (dim - Chunk::log_elements_in_chunk);
	}

	/*
	 * Calculate on CPU
	 */
	Element* e = (Element*) malloc(
			sizeof(Element) * (p_len << Chunk::log_elements_in_chunk));
	Chunk::chunkToNormal(p, (Elements_Chunk*) e, p_len, true);
	Element t;
	for(unsigned int i = 0 ; i < (p_len<<Chunk::log_elements_in_chunk) ; i+=2){
		Element::c_mul(e[i+1],fft->lastShift,t);
		Element::c_mul(e[i+1],fft->lastD,e[i+1]);
		Element::c_add(e[i+1],e[i],e[i+1]);
		Element::c_add(e[i+1],t,e[i+1]);
		Element::c_add(e[i],t,e[i]);
	}

	/*
	 * Change the calculated elements to Chunk representation
	 */
	Chunk* h_cpu = (Chunk*) malloc(sizeof(Chunk) * p_len);
	Chunk::normalToChunk((Elements_Chunk*) e, h_cpu, p_len, true);


	/*
	 * Calculate using GPU
	 */
	Chunk* d_p = copyToGPU(p, p_len);
	Chunk* d_sub;
	cudaMalloc(&d_sub, sizeof(Chunk) * p_len);
	Chunk* h_p;
#ifdef __GNUC__
	timespec start, end;
	clock_gettime(CLOCK_REALTIME, &start);
#endif	// #ifdef __GNUC__
	GPU_FFT::linearEvaluation(d_p,d_sub,p_len);
	h_p = copyFromGPU(d_p,p_len);
	for (unsigned int i = 0; i < p_len; ++i) {
		for (unsigned int j = 0; j < Chunk::cells_in_chunk; ++j) {
			if (h_cpu[i].v[j] != h_p[i].v[j]) {
				flag = false;
			}
		}
	}
	cudaFree(d_sub);
	cudaFree(d_p);
	free(h_cpu);
	free(h_p);
	free(p);
	free(e);
	delete fft;
	return flag;
}
#endif
} /* namespace FFF */
