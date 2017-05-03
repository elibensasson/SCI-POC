/**
 *       @file  CXX11_macros.hpp
 *      @brief  This file contains some mocros
 *      for compatibility between gcc C++11 support and MS VS2010 tr1 support.
 *      This file is written at 2013.08.22, and should be temporal untill
 *      C++11 will be supported widely enough by some released version of
 *      MS Visual Studio.
 *      Writing macros for such things is not the best solution, but can be used
 *      if no easy way to overcome some gap exists.
 *      Each macro should be well documented and explain exactly what is the
 *      unsupported feature it is waiting for.
 *
 *     @author  Michael Riabzev, RiabzevMichael@gmail.com
 * =====================================================================================
 */

/**
 * C++11 presents the 'deleted constructor' syntax.
 * It used to notify the compiler that some constructor
 * that has a default implementation should not exist.
 * The to simulate something like that in C++97 would
 * implement that specific constructor as private.
 * Such a solution is not complete because it may be still
 * called by friend or member methods, this is why to be
 * completely sure it is never called one would implement it
 * to throw some exception.
 * 
 * An example for such a usage is when one wants his class
 * instances only to be move (using the move constructor,
 * which is enabled in tr1 on VS2010) and never copied.
 * In such case he would define a default move constructor
 * for his class, and delete its copy constructor.
 *
 * In C++11 such a case would have been implemented this way:
 *
 * class A {
 * public:
 * 	A() = default;
 * 	A(A&&) = default;
 * 	A(const A&) = deleted;
 * };
 *
 * In tr1 only supported compiler (as VS2010) such a case
 * would have been implemented this way:
 *
 * class A {
 * public:
 * 	A(){};
 * 	A(A&& src){};
 * private:
 *	A(const A& src){throw "deleted function";}
 * };
 *
 * I tried to implement that for tr1, but it seems that
 * g++ with C++11 support handles the fact that the
 * copy constructor exists but private (and not deleted)
 * by throwing an error in case there are places it
 * thinks suits for it to be used, instead of just
 * using the move constructor instead. This seems a
 * little weird an I did not check to much on that,
 * but this macro just fixes this.
 *
 * Usage example:
 * class A {
 * private:
 * 	A(const A& src) _COMMON_CXX11_DELETED
 * };
 */

 #include "common/Infrastructure/Infrastructure.hpp"

#if __cplusplus > 199711L
#define _COMMON_CXX11_DELETED = delete;
#else
#define _COMMON_CXX11_DELETED ; // {_COMMON_FATAL("Attempted to call deleted function.");} // Error
    // message commented out because it causes runtime check instead of static (linker) check. 
#endif

/**
 * explicit casting (conversion operator) is allowed only from C++11.
 * I do not know any alternative for VS2010.
 * It is still a good practice to use it, unfortunately
 * the macro has no affect on VS2010 compilation.
 */
#if __cplusplus > 199711L
#define _COMMON_CXX11_EXPLICIT explicit
#else
#define _COMMON_CXX11_EXPLICIT
#endif


#ifndef  __CXX11_macros_HPP
#define  __CXX11_macros_HPP

#endif   // __CXX11_macros_HPP
