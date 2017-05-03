//#include "Element.h"
//#include "Tests.h"
//#include <iostream>
//#include <omp.h>
//using namespace FFF;
//using namespace std;
//#include <climits>
//#include <cstring>
//#include <cstdio>
//#include <cstdlib>
////int measureFFT() {
////	unsigned int threads[] = { 1, 2, 4, 8 };
////	int threads_len = 4;
////	Tests::measureFFT(3, 22, threads, threads_len);
////}
////int measureiFFT() 
////	unsigned int threads[] = { 1, 2, 4, 8 };
////	int threads_len = 4;
////	Tests::measureiFFT(3, 30, threads, threads_len);
////}
//int main(int argc, char** argv) {
//#ifdef __GPU
//	cudaSetDevice(1);
//	Chunk::setMod();
//#endif	// #ifdef __GPU
////	measureiFFT();
//	if (argc == 4) {
//		int size = atoi(argv[3]);
//		omp_max_threads = atoi(argv[2]);
//		log_omp_max_threads = FFF::log[omp_max_threads];
//		if (!strcmp("i", argv[1]))
//			cout << Tests::testiFFT(size) << endl;
//		else
//			cout << Tests::testFFT(size, 1 << size) << endl;
//		return 0;
//	}
//	omp_max_threads = 8;
//	cout << Tests::testFFT(7,1<<7);
////	std::cout << Tests::testGPUMultiExponentiate(24,12);
////	for (unsigned int i = 6; i < 30; ++i) {
////		cout << i << endl;
////		Tests::testFFT(i, 1 << i);
////	}
////	unsigned int k = 15;
////	for(unsigned int i = k-1 ; i >= 2 ; --i){
////		std::cout << Tests::testGPUMultiExponentiate(k,i) << std::endl;
////		std::cout << Tests::testPartition_GPU(k,i) << std::endl;
////		std::cout << Tests::testTaylorExpansion_GPU(k,i) << std::endl;
////		std::cout << Tests::testLinearEvaluation_GPU(i) << std::endl;
////		std::cout << Tests::TestWFromUV_GPU(k,i) << std::endl;
////	}
////	len_t b_size = 3;
////	Element* es = (Element*)malloc(sizeof(Element)*b_size);
////	Tests::generateBasis(es,b_size);
////	Element l;
////	Element::c_setZero(l);
////	Basis b(es,b_size,l);
////	cout << Tests::test_generateSubspace(b);
////	double nSecSum=0;
////	for(unsigned int i =0  ; i < (1<<20) ; ++i)
////		nSecSum +=Tests::measure_cMul();
////	cout << nSecSum/(1<<20);
////	cout << Tests::test_cMul(1000000);
//}
