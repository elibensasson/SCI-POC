
#ifndef NTL_WordVectorStatic__H
#define NTL_WordVectorStatic__H

/************************************************************************************************************/
/**
 * @Documentation - Partial.
 *
 * The WordVectorStatic class is a static implementation of NTL's WordVector class, made by the PCPCD team.
 * The static behavior should improve the efficiency of the class, by avoiding memory allocations and deallocations.
 * Instead of using a dynamic array of unsigned long, the array will be static, thus setting a limit to the 
 * polynomial's degree.
 * 
 *   A WordVector is functionally similar to
 *   a  generic NTL vector of _ntl_ulong.  
 * 
 * Be careful! the MaxLength() function does not return the max length ever set, but rather the max space allocated,
 * which *may* be more.
 * 
 * The FixLength() facility is not available.
 * 
 * The reason for special-casing is efficiency (of course).
 *
 */
  /************************************************************************************************************/

#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

NTL_OPEN_NNS


#ifndef NTL_RANGE_CHECK
#define NTL_WV_RANGE_CHECK_CODE 
#else
#define NTL_WV_RANGE_CHECK_CODE if (i < 0 || !dynamicRep || i >= long(dynamicRep[-1])) RangeError(i);
#endif

//#ifndef NTL_WV_STATIC_RANGE_CHECK_CODE
//#define NTL_WV_STATIC_RANGE_CHECK_CODE 
//#else
#define NTL_WV_STATIC_RANGE_CHECK_CODE if (i < 0 || !statRep || i >= STATIC_VECTOR_LENGTH) RangeError(i);
//#endif

// vectors are allocated in chunks of this size

#ifndef NTL_WordVectorMinAlloc
#define NTL_WordVectorMinAlloc (4)
#endif

// vectors are always expanded by at least this ratio

#ifndef NTL_WordVectorExpansionRatio
#define NTL_WordVectorExpansionRatio (1.2)
#endif

// controls initialization during input

#ifndef NTL_WordVectorInputBlock
#define NTL_WordVectorInputBlock 50
#endif

// controls release functionality

#define NTL_RELEASE_THRESH (10000)
// #define NTL_RELEASE_THRESH (0)

//We define the number of words in the vector. 
//Since each long is 32 bit, the max degree of the polynomial is limited to 32*STATIC_VECTOR_LENGTH - 1.
#define STATIC_VECTOR_LENGTH 2

/**
 * The class WordVectorStatic implements a statically allocated vector of words.
 * We use this vector to represent coefficients of polynomials over GF(2), specifically, each
 * bit of a word is a coefficient. 
 * The class is used by the GF2X type.
 */
class WordVectorStatic {  

public: 

	_ntl_ulong *dynamicRep;						///A dynamic array of unsigned longs. Each bit is a GF2 element. dynamicRep[-1] = array length
	_ntl_ulong statRep[STATIC_VECTOR_LENGTH];	///A static array of unsigned longs. Each bit is a GF2 element.
	_ntl_ulong staticMode;						///ntl_ulong used as a flag: 1 iff using static array, 0 otherwise.

	void RangeError(long i) const;				///Prints an error of sort "index out of range in vector: i (len)"

	
	/***********************************/
	/******** Class Constructors *******/
	/***********************************/

	///Default C'tor
	WordVectorStatic() { 
		staticMode = 1;
		for (int i=0; i<STATIC_VECTOR_LENGTH; i++)
			statRep[i] = 0;
		
		dynamicRep = 0; 
	}	

	//constructs a new object out of x.
	WordVectorStatic(WordVectorStatic& x, INIT_TRANS_TYPE) { 
		
		this->staticMode = x.staticMode;
		if(staticMode) {
			for (int i=0; i<STATIC_VECTOR_LENGTH; i++)
				statRep[i] = x.statRep[i];
			dynamicRep = 0;
		}

		else {
			dynamicRep = x.dynamicRep; 
			x.dynamicRep = 0; 
		}
	}	

	///Constructs a vector of n elements.
	WordVectorStatic(INIT_SIZE_TYPE, long n) { 
		if (n > STATIC_VECTOR_LENGTH) {
			dynamicRep = 0; 
			DoSetLength(n); 
			staticMode = 0;	
		}
		else {
			for (int i=0; i<STATIC_VECTOR_LENGTH; i++)
				statRep[i] = 0;
			staticMode = 1;
			dynamicRep = 0;
		}
	}	

	///Copy Constructor
	WordVectorStatic(const WordVectorStatic& a) { 
		this->staticMode = a.staticMode;
		
		if (staticMode) {
			for (int i=0; i<STATIC_VECTOR_LENGTH; i++)
				statRep[i] = a.statRep[i];		
			dynamicRep = 0;
		}

		else {
			dynamicRep = 0;
			*this = a; 
		}
	}     

	///Assignment operator
	WordVectorStatic& operator=(const WordVectorStatic& a);  

	//Class Destructor
	~WordVectorStatic();  
	void kill(); 

   ///Release the memory. In the static case - nothing ot do.
   void release() { 
	   //if (!staticMode)
			if (MaxLength() > NTL_RELEASE_THRESH) 
				kill(); 
   }
   // this conditinally kills the vector, if its size is excessive

	/*************************************/
	/******** Set Length Fucntions *******/
	/*************************************/

   ///Sets the vector's length. In the static case - nothing to do but parameters check...
   void DoSetLength(long n);
  
   ///Sets the vector's length. In the static case - nothing to do but parameters check...
   void SetLength(long n);

   ///Resets the vector
   void ZeroLength();

   	void SetMaxLength(long n); 

	void QuickSetLength(long n) { 
		if (!staticMode)
			dynamicRep[-1] = _ntl_ulong(n); 
	} 

	///Return the vector length (in the static case - constant).
	long length() const { 
		
		//Static mode - compute length:
		if (staticMode) {
			long res = STATIC_VECTOR_LENGTH;
			while (statRep[res-1] == 0) {
				if (res==1) return 0;
				res--;
			}
			return res;
		}
		
		//Dynamic mode - use dynamicRep[-1]
		return (!dynamicRep) ?  0 : long(dynamicRep[-1]); 
	}  

	///Return the vector max length (in the static case - constant).
	long MaxLength() const { 
		//Static mode:
		if (staticMode)
			return STATIC_VECTOR_LENGTH;
		
		//Dynamic mode:
		return (!dynamicRep) ?  0 : long(dynamicRep[-2] >> 1); 
	} 

	/** Converts the object to a dynamic one, keeping the value as the original. */
	void convertToDynamic();

	/***********************************/
	/******** Operators and Misc *******/
	/***********************************/

	///Access the relevant word in the vector.
	_ntl_ulong& operator[](long i)   
	{  
		//Static mode:
		if (staticMode) {
			NTL_WV_STATIC_RANGE_CHECK_CODE  ///Check the validity of i.
			return statRep[i];			///Return the i'th word.
		}

		//Dynamic mode:
		NTL_WV_RANGE_CHECK_CODE  
		return dynamicRep[i];  
	}  

	///Access the relevant word in the vector (read only).
	const _ntl_ulong& operator[](long i) const 
	{  
		//Static mode:
		if (staticMode) {
			NTL_WV_STATIC_RANGE_CHECK_CODE  ///Check the validity of i.
			return statRep[i];			///Return the i'th word.
		}

		//Dynamic mode:
		NTL_WV_RANGE_CHECK_CODE  
		return dynamicRep[i];  
	}  

	_ntl_ulong& operator()(long i) { return (*this)[i-1]; }  
	const _ntl_ulong& operator()(long i) const { return (*this)[i-1]; } 

	const _ntl_ulong* elts() const { 
		if (staticMode)		//Static case
			return statRep;
		
		return dynamicRep;			//Dynamic case
	}  

	_ntl_ulong* elts() { 
		if (staticMode)		//Static case
			return statRep;

		return dynamicRep;			//Dynamic case
	}

	static void swap_impl(WordVectorStatic& x, WordVectorStatic& y);  
	static void append_impl(WordVectorStatic& v, _ntl_ulong a); 
	static void append_impl(WordVectorStatic& v, const WordVectorStatic& w); 
}; 

	/********************************/
	/******** Swap & InAppend *******/
	/********************************/

inline void swap(WordVectorStatic& x, WordVectorStatic& y) 
{ WordVectorStatic::swap_impl(x, y); }

inline void append(WordVectorStatic& v, _ntl_ulong a)
{ WordVectorStatic::append_impl(v, a); }

inline void append(WordVectorStatic& v, const WordVectorStatic& w)
{ WordVectorStatic::append_impl(v, w); }

	/********************************/
	/******** Input \ Output ********/
	/********************************/

NTL_SNS istream& operator>>(NTL_SNS istream&, WordVectorStatic&);  
NTL_SNS ostream& operator<<(NTL_SNS ostream&, const WordVectorStatic&);  


long operator==(const WordVectorStatic& a, const WordVectorStatic& b);  
long operator!=(const WordVectorStatic& a, const WordVectorStatic& b);


long InnerProduct(const WordVectorStatic& a, const WordVectorStatic& b);

//void ShiftAdd(_ntl_ulong *cp, const _ntl_ulong* ap, long sa, long n);
// cp = cp + (a << n)

/******************************************************************/
/***** Block Constructors \ Destructors. Used by NTL vectors ******/
/******************************************************************/

long WV_BlockConstructAlloc(WordVectorStatic& x, long d, long n);

void WV_BlockConstructSet(WordVectorStatic& x, WordVectorStatic& y, long i);

long WV_BlockDestroy(WordVectorStatic& x);

//long WV_storage(long d);





NTL_CLOSE_NNS

#endif	//NTL_WordVectorStatic__H
