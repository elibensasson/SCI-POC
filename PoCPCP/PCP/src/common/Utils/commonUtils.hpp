/********************************************** commonUtils.hpp ************************************************/
/** @file.
 *
 * The file commonUtils.hpp contains various utility functions and date type definitions used by the different 
 * classes in the project.
 * 
 * Main classes defined in the file:
 * ElementsVector - Represents a vector of field elements, which can store a large number of elements (larger than 2^32).
 * DataString - Represents a std::string of data used as the instance input for Prover_Lang and Verifier-Lang.
 * ErrorHandling - Implements the functionality of displaying the content of error messages and exiting the program.
 * 
 * Main date types defined in the file:
 * _COMMON_Long - A 64-bit number used defined for the _COMMON code.
 * BasisIndex - An index of a basis element.
 * FieldIndex - An index of a field element.
 * FieldPoint - An element in an extension field of GF(2).
 * FieldPointPair - A two-dimensional point (x,y) in the field.
 * GF2Poly - Univariate polynomial over GF(2) type.
 *
 * Main functions defined in the file:
 * Basic Math Functions - Power macro, log with basis 2, extract a bit from an integer etc.
 * GF2E Functionality - Convert an integer to a vector of GF2 bits and vice versa, compare two GF2 extension field elements etc.
 * DataString Functions - Create a std::string, append two std::strings, extract a specific bit etc.
 * Memory and Performance - Print information regarding the global memory usage and the memory usage of the running process.
 * ErrorHandling Functions - The function PCP_Error which prints an error message and exits.
 */
  /************************************************************************************************************/

#ifndef _COMMON_UTILS_HPP_
#define _COMMON_UTILS_HPP_

#include <cstdio>  
#include <fstream>
#include <cstdlib>
#include <vector>
#include <list>
#include <cassert>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>
#include <ctime>

#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/lexical_cast.hpp>

// TODO remove the following includes to get rid of dependencies
#include "FileHandling.hpp"
#include "common/EnvParams.hpp"
#include "Logging.hpp"

#ifdef __GLIBC__
#include <execinfo.h> // backtraces
#endif


#ifndef _MSC_VER // emulate the MSVC-specific sprintf_s using the standard snprintf
#define sprintf_s snprintf //TODO: sprintf_s!=snprintf (http://blog.verg.es/2008/09/sprintfs-is-not-snprintf.html)
#endif

#if defined _WIN64 || defined _WIN32
// To use this macro you would have to include windows.h before
// including this file
#ifdef __MEASURE_TIMES
#define MEASURE_TIME(CODE, PRINT) do {									\
	LARGE_INTEGER counterStart;											\
	LARGE_INTEGER counterEnd;											\
	QueryPerformanceCounter(&counterStart);								\
	{ CODE; }															\
	QueryPerformanceCounter(&counterEnd);								\
	LARGE_INTEGER total;												\
	total.QuadPart = counterEnd.QuadPart - counterStart.QuadPart;		\
	if (PRINT) {														\
		std::cout														\
			<< "*** Code measured: "									\
			<< #CODE													\
			<< "; -> "													\
			<< "threads: "												\
			<< FFF::omp_max_threads										\
			<< " -> "													\
			<< (total.QuadPart)											\
			<< " micro seconds"											\
			<< std::endl;												\
	}																	\
} while (0)
#else	//#ifdef __MEASURE_TIMES
#define MEASURE_TIME(CODE, PRINT) { CODE; }
#endif	//#ifdef __MEASURE_TIMES
#else	// #if defined _WIN64 || defined _WIN32
#define MEASURE_TIME(CODE, PRINT) { CODE; }
#endif	// #if defined _WIN64 || defined _WIN32

namespace PCP_Project {

/********************************************************/
/******************* Error Handling *********************/
/********************************************************/

// declare a function as never returning, to quiet down "control reaches end of non-void function" warnings
#if defined(_MSC_VER) // VisualC++
  #define __noreturn _declspec(noreturn)
#elif defined(__GNUC__)
  #define __noreturn __attribute__((noreturn))
#else
  #define __noreturn
#endif



/**
 * The ErrorHandling class containimplements the functionality of displaying the content of error
 * messages (including content of call stack when error happened), and exiting the program.
 */
class ErrorHandling {
public:
    static void __noreturn fatalError(std::string msg);
    static void __noreturn fatalError(const std::stringstream& msg);
	static void warning(std::string msg);
	static void warning(const std::stringstream& msg);
	static void info(std::string msg);
	static void info(const std::stringstream& msg);
	static void printStacktrace();

};
    


} // of namespace

#define ASSERT(COND, ...)								\
do {													\
    if (!(COND)) {										\
        std::stringstream msg;							\
        msg << __VA_ARGS__ << " (" << #COND	<< " at " << __FILE__ << ":" << __LINE__ << " " << ")" << std::ends;	\
        PCP_Project::ErrorHandling::fatalError(msg);	\
    }													\
} while (0)

#endif	//_COMMON_UTILS_HPP_
