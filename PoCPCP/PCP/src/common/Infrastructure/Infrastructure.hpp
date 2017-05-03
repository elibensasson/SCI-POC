#ifndef  __Infrastructure_HPP
#define  __Infrastructure_HPP

#include <cstdint>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdarg>

#ifndef _MSC_VER // emulate the MSVC-specific sprintf_s using the standard snprintf
#define sprintf_s snprintf //TODO: sprintf_s!=snprintf (http://blog.verg.es/2008/09/sprintfs-is-not-snprintf.html)
#endif

#ifdef _DEBUG // MSVC Debug build
#define DEBUG // gcc Debug flag
#endif

/********************************************************/
/**************** Class Writing Helpers *****************/
/********************************************************/
// A macro to disallow any non-defined constructors
// This should be used in the private: declarations for a class
#define DISALLOW_CONSTRUCTION(TypeName) \
  TypeName();               

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

/********************************************************/
/*************** Debug String Formatting ****************/
/********************************************************/

namespace PCP_Project {
// someday, if/when MSVC supports C++0x variadic templates, change FMT in release version to the
// following in order to increase efficiency:
// #define FMT(...) ""
::std::string FMT(const char* format, ...);
} // namespace PCP_Project

namespace Infrastructure {
/** Safely converts 64-bit types to 32-bit, or from unsigned to signed */
long safeConvert(const int64_t num);

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
            static void __noreturn fatalError(const ::std::string& msg);
            static void __noreturn fatalError(const std::stringstream& msg);
            static void printStacktrace();

    };

#define _COMMON_DEBUG_MSG(msg) do {  \
            std::cerr << msg << " (In file " << __FILE__ << " line " << __LINE__ << ".)"; \
        } while (0)


#define _COMMON_FATAL(msg) do {  \
            ::std::stringstream msgStream; \
            msgStream << msg << " (In file " << __FILE__ << " line " << __LINE__ << ".)"; \
            ::Infrastructure::ErrorHandling::fatalError(msgStream.str()); \
        } while (0)

// TODO change _COMMON_ASSERT to not run in debug
#define _COMMON_ASSERT(predicate, msg) if(!(bool(predicate))) _COMMON_FATAL(msg);

/********************************************************/
/****************** Basic Math **************************/
/********************************************************/

//Calculates ::Infrastructure::Log2 of a number.
double Log2(double n);

//Calculates  upper bound of Log2 of a number (number of bits needed to represent value)
unsigned int Log2ceil(uint64_t i);

//Returns true iff the given number is a power of 2.
bool IsPower2(const long x);


//Returns a^b when a can be a and b are INTEGERS.
//constexpr int64_t POW(int64_t base, int exponent) {
//	return (int64_t) powl((long double)base, (long double)exponent);
//}
//#define POW(a,b) ((int64_t)(pow((float)(a),(int)(b))))

// Returns 2^exponent
/*constexpr*/ inline int64_t POW2(int exponent) {
    //assert(exponent>=0);
    return ((int64_t)1) << exponent;
}

//Returns the ceiling of a when a is of type double.
/*constexpr*/ inline int64_t CEIL(double a) {
    return (int64_t)ceil(a);
}
//#define CEIL(a)  ((int64_t)ceil((double)(a)))

} // namespace Infrastructure 

#endif   // __Infrastructure_HPP
