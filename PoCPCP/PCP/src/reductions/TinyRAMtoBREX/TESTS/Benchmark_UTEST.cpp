
/**
*      @file  Benchmark_UTEST.cpp
*      @brief  benchmark simple operation to compare different
*              operations on different data constructs.
*
*	It expected to be linked with an implementation for
*	FieldElement and Element class.
*
*	Additionally uses libraries NTL, FFT and GTEST.
*
*     @author  Lior Greenblatt, lior.greenblatt@gmail.com
* =====================================================================================
*/

#include "common/Configuration.hpp"
#include "algebraLib/FieldElement.hpp"
#include "algebraLib/FiniteField.hpp"
#include <NTL/GF2XFactoring.h>
#include <gtest/gtest.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/format.hpp>
#include "common/Utils/Timing.hpp"


namespace{
using FFF::Element;
using std::vector;

#ifdef _MSC_VER // Visual C++
	/* Workaround (shaul): This #ifdef is due to (yet another) violation of the OpenMP standard in the VC++ compiler.
	The threadprivate clause in VC++ is implemented directly as the application of the thread __declspec attribute */
	extern "C" __declspec(thread) unsigned int ntlGlobalThreadNum;  /**< Used to pass a unique thread number to NTL for memory management purposes
																	Used by GF2XRegister in NTL/GF2X1.cpp
																	Defined as extern "C" to avoid name mangling when used with C++.  */
#else
	extern "C" unsigned int ntlGlobalThreadNum;                     /**< Used to pass a unique thread number to NTL for memory management purposes
																	Used by GF2XRegister in NTL/GF2X1.cpp
																	Defined as extern "C" to avoid name mangling when used with C++.  */
#pragma omp threadprivate(ntlGlobalThreadNum)
#endif

	extern unsigned int globalFFTNumthreads;     /**< should be set to 0 (automatically sets to number of processors),
													 change to a positive integer to force a number of threads in
													 GathenGerhardInterpolation() and GathenGerhardEvaluation(). For
													 debug and analysis purposes. */

	extern "C" const unsigned int ntlMaxNumThreads;  /**< defines the maximum number of threads supported by multi-threaded version of NTL (change this in GF2X1.cpp) */


typedef double usec_t;
std::map<std::string, usec_t> measuredTimes;

const int BMARK_ITERATIONS = 1000000;
// Use globals so compiler will not optimize out

template<class T=int>
std::string getTestName()
{
    const ::testing::TestInfo* const test_info =
		::testing::UnitTest::GetInstance()->current_test_info();
	std::string testName = test_info->test_case_name();
	testName += ".";
	testName += test_info->name();
	return testName;
}

template<class T>
void measure_addition(const vector<T>& vec, int tnum)
{
	vector<T> res(vec.size()-1);
	std::string testName = getTestName();
	{
		PCP_Project::Timer t;
	    {
#pragma omp parallel for num_threads(tnum)
            for(int i=0; i< vec.size()-1; i++) {
                const size_t currThreadId = omp_get_thread_num();
                ntlGlobalThreadNum = currThreadId;
                res[i] = vec[i]+vec[i+1];
            }
        }
		measuredTimes[testName] = t.getElapsed();
	}
}


template<class T>
void measure_multiplication(const vector<T>& vec, int tnum)
{
	vector<T> res(vec.size()-1);
	std::string testName = getTestName();
	{
		PCP_Project::Timer t;
	    {
#pragma omp parallel for num_threads(tnum)
            for(int i=0; i< vec.size()-1; i++) {
                const size_t currThreadId = omp_get_thread_num();
                ntlGlobalThreadNum = currThreadId;
                res[i] = vec[i]*vec[i+1];
            }
        }
		measuredTimes[testName] = t.getElapsed();
	}
}

#define NTL_BENCHMARK(TNUM, OPSTR)										\
TEST(Benchmark, ntl_ ## OPSTR ## _ ## TNUM ## _thread) {							\
	std::vector<NTL::GF2E> vec(BMARK_ITERATIONS+1);                                 \
    for(auto& e : vec) e= NTL::GF2EInfo->random_GF2E();						\
	measure_ ## OPSTR (vec, TNUM);	\
}

#define NTL_BENCHMARK_GROUP(TNUM)				\
	NTL_BENCHMARK(TNUM, addition);		\
	NTL_BENCHMARK(TNUM, multiplication);

NTL_BENCHMARK_GROUP(1);
NTL_BENCHMARK_GROUP(2);
NTL_BENCHMARK_GROUP(4);
NTL_BENCHMARK_GROUP(8);
NTL_BENCHMARK_GROUP(16);
NTL_BENCHMARK_GROUP(32);
NTL_BENCHMARK_GROUP(64);
NTL_BENCHMARK_GROUP(128);

#define FFT_BENCHMARK(TNUM, OPSTR)							\
	TEST(Benchmark, fft_ ## OPSTR ## _ ## TNUM ## _thread) {		\
		std::vector<Element> vec(BMARK_ITERATIONS+1);               \
        for(auto& e : vec)e.c[0] = rand();                                       \
        measure_ ## OPSTR (vec, TNUM);	\
	}

#define FFT_BENCHMARK_GROUP(TNUM)		\
FFT_BENCHMARK(TNUM, addition);		\
FFT_BENCHMARK(TNUM, multiplication);

FFT_BENCHMARK_GROUP(1);
FFT_BENCHMARK_GROUP(2);
FFT_BENCHMARK_GROUP(4);
FFT_BENCHMARK_GROUP(8);
FFT_BENCHMARK_GROUP(16);
FFT_BENCHMARK_GROUP(32);
FFT_BENCHMARK_GROUP(64);
FFT_BENCHMARK_GROUP(128);

TEST(Benchmark, print_results) {
	std::map<std::string, usec_t> t = measuredTimes;
	std::cout << "Iterations number: " << BMARK_ITERATIONS << std::endl;
	std::cout << boost::format("%-15s") % "operation" << boost::format("%-9d") % "threads" <<
		boost::format("%-13d") % "ntl (secs)" << boost::format("%-13d") % "fft (secs)" << std::endl;
	for (std::map<std::string, usec_t>::iterator it_ntl = measuredTimes.begin(); it_ntl != measuredTimes.end(); ++it_ntl) {
		for (std::map<std::string, usec_t>::iterator it_fft = measuredTimes.begin(); it_fft != measuredTimes.end(); ++it_fft) {
			std::string ntl_key = it_ntl->first;
			std::string fft_key = it_fft->first;
			std::vector<std::string> ntl_tokens;
			std::vector<std::string> fft_tokens;
			boost::split(ntl_tokens, ntl_key, boost::is_any_of("_"));
			boost::split(fft_tokens, fft_key, boost::is_any_of("_"));
			if (ntl_tokens[0] == "Benchmark.ntl" && fft_tokens[0] == "Benchmark.fft" &&
				ntl_tokens[1] == fft_tokens[1] && ntl_tokens[2] == fft_tokens[2]) {
				usec_t ntl_value = it_ntl->second;
				usec_t fft_value = it_fft->second;
				std::cout << boost::format("%-15s") % fft_tokens[1] << boost::format("%-9d") % fft_tokens[2] <<
					boost::format("%-13d") % ntl_value << boost::format("%-13d") % fft_value << std::endl;
			}
		}
	}
}
}
