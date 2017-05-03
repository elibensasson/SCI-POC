

#ifndef NTL_ZZ_p__H
#define NTL_ZZ_p__H

#include <NTL/ZZ.h>
#include <NTL/ZZVec.h>

NTL_OPEN_NNS


// representation:  each ZZ_p is represented by a ZZ in the range 0..p-1.

// The constructor for a ZZ_p pre-allocates space for the underlying ZZ,
// and initializes it to zero.

const int MAX_ZZ_p_TEMPS = 16;

class ZZ_p;

class ZZ_pInfoT {
private:
   ZZ_pInfoT();                       // disabled
   ZZ_pInfoT(const ZZ_pInfoT&);   // disabled
   void operator=(const ZZ_pInfoT&);  // disabled
public:
   ZZ_pInfoT(const ZZ& NewP);
   ~ZZ_pInfoT();

  long ref_count;		// reference count for gargabge collection

  ZZ p;				// the modulus
  long size;			// p.size()
  long ExtendedModulusSize;

  // the following implement a "lazy" initialization strategy
  long initialized;		// flag if initialization really was done
   void init();
   void check() { if (!initialized) init(); }

  long NumPrimes;

  long MaxRoot;

  long QuickCRT;

  ZZ MinusMModP;		//  -M mod p, M = product of primes

// mmsize from multi-thread version of NTL, probably not necessary //shaul
//  long mmsize;			// size pre-computed for MultiMul

// From original NTL version, not in multi-threaded version from LinBox //shaul
   void *crt_struct;
   void *rem_struct;


  // the following arrays are indexed 0..NumPrimes-1
  // q = FFTPrime[i]


//From multi-thread version of NTL
  ZZVec CoeffModP;		// coeff mod p

  double *x;			// u/q, where u = (M/q)^{-1} mod q
  long *u;			// u, as above

  double **tbl;			// table used for MultiRem; only with NTL_SINGLE_MUL

  long **tbl1;			// table used for MultiRem; only with NTL_TBL_REM

  ZZ_p *temps[MAX_ZZ_p_TEMPS];
  long temps_top;

    // Methods to perform various operations on ZZ_p's with this field def
    // These allow operations to be performed in multiple fields without having
    // to save and restore the context, as well as in different fields in 
    // parallel.
    
    void conv (ZZ_p & x, const ZZ & a);
    void conv (ZZ_p & x, long a);
    void add (ZZ_p & x, const ZZ_p & a, const ZZ_p & b); // x = a + b
    void add (ZZ_p & x, const ZZ_p & a, long b);
    inline void add (ZZ_p & x, long a, const ZZ_p & b) { add (x, b, a); }
    void sub (ZZ_p & x, const ZZ_p & a, const ZZ_p & b); // x = a - b
    void sub (ZZ_p & x, const ZZ_p & a, long b);
    void sub (ZZ_p & x, long a, const ZZ_p & b);
    void negate (ZZ_p & x, const ZZ_p & a); // x = -a
    void mul (ZZ_p & x, const ZZ_p & a, const ZZ_p & b); // x = a*b
    void mul (ZZ_p & x, const ZZ_p & a, long b);
    inline void mul (ZZ_p & x, long a, const ZZ_p & b) { mul (x, b, a); }
    void sqr (ZZ_p & x, const ZZ_p & a); // x = a^2
    ZZ_p sqr (const ZZ_p & a);
    void div (ZZ_p & x, const ZZ_p & a, const ZZ_p & b);
    void inv (ZZ_p & x, const ZZ_p & a);
    ZZ_p inv (const ZZ_p & a);
    void div (ZZ_p & x, const ZZ_p & a, long b);
    void div (ZZ_p & x, long a, const ZZ_p & b);
    void power (ZZ_p & x, const ZZ_p & a, const ZZ & e);
    inline ZZ_p power (const ZZ_p & a, const ZZ & e);
    void power (ZZ_p & x, const ZZ_p & a, long e);
    inline ZZ_p power (const ZZ_p & a, long e);
    void random (ZZ_p & x); // x = random element in ZZ_p
    inline ZZ_p random_ZZ_p ();
};

extern ZZ_pInfoT *ZZ_pInfo;	// info for current modulus, initially null



class ZZ_pContext {
private:
ZZ_pInfoT *ptr;
public:
void save();
void restore() const;

ZZ_pContext() { ptr = 0; }
ZZ_pContext(const ZZ& p);

ZZ_pContext(const ZZ_pContext&); 


ZZ_pContext& operator=(const ZZ_pContext&); 


~ZZ_pContext();


};


class ZZ_pBak
{
  private:long MustRestore;
  ZZ_pInfoT *ptr;

    ZZ_pBak (const ZZ_pBak &);	// disabled
  void operator = (const ZZ_pBak &);	// disabled

    public:void save ();
  void restore ();

    ZZ_pBak ()
  {
    MustRestore = 0;
    ptr = 0;
  }

   ~ZZ_pBak ();


};




struct ZZ_p_NoAlloc_type
{
  ZZ_p_NoAlloc_type ()
  {
  }
};
const ZZ_p_NoAlloc_type ZZ_p_NoAlloc = ZZ_p_NoAlloc_type ();


class ZZ_pTemp
{
  private:long pos;

    public:ZZ_pTemp ();
   ~ZZ_pTemp ();

    ZZ_p & val () const;
};

#define NTL_ZZ_pRegister(x)  \
   ZZ_pTemp ZZ_pTemp__ ## x; ZZ_p& x = ZZ_pTemp__ ## x . val()


class ZZ_p {

public:

ZZ _ZZ_p__rep;
  
  friend class ZZ_pInfoT;

// static data

  static void init(const ZZ&);


typedef void (*DivHandlerPtr)(const ZZ_p& a);   // error-handler for division
  static DivHandlerPtr DivHandler;


// ****** constructors and assignment

ZZ_p();

ZZ_p(const ZZ_p& a) :  _ZZ_p__rep(INIT_SIZE, ZZ_pInfo->size) { _ZZ_p__rep = a._ZZ_p__rep; }

ZZ_p(ZZ_p_NoAlloc_type) { }  // allocates no space

~ZZ_p() { } 

ZZ_p& operator=(const ZZ_p& a) { _ZZ_p__rep = a._ZZ_p__rep; return *this; }


// You can always access the _ZZ_p__representation directly...if you dare.
ZZ& LoopHole() { return _ZZ_p__rep; }

ZZ_p(ZZ_p& x, INIT_TRANS_TYPE) : _ZZ_p__rep(x._ZZ_p__rep, INIT_TRANS) { }



static const ZZ& modulus() { return ZZ_pInfo->p; }
static long ModulusSize() { return ZZ_pInfo->size; }
static long storage() { return ZZ_storage(ZZ_pInfo->size); }

//From the thread-safe version of LinBox by "hovinen" this is the function StorageSize(). Not used by our original version //shaul
//static long StorageSize() { return ZZ_pInfo->p.size () + 4; }

static const ZZ_p& zero();


ZZ_p(INIT_VAL_TYPE, const ZZ& a);
ZZ_p(INIT_VAL_TYPE, long a);

// read-only access to _ZZ_p__representation
friend const ZZ& rep(const ZZ_p& a) { return a._ZZ_p__rep; }

// ****** conversion

  inline friend void conv (ZZ_p & x, const ZZ & a) { ZZ_pInfo->conv (x, a); }

  

  friend ZZ_p to_ZZ_p (const ZZ & a)
  {
    return ZZ_p (INIT_VAL, a);
  }


  inline friend void conv (ZZ_p & x, long a) { ZZ_pInfo->conv (x, a); }

  

  friend ZZ_p to_ZZ_p (long a)
  {
    return ZZ_p (INIT_VAL, a);
  }


// ****** assignment

  ZZ_p & operator = (long a)
  {
    conv (*this, a);
    return *this;
  }


// ****** some basics


  friend void clear (ZZ_p & x)
// x = 0
  {
    clear (x._ZZ_p__rep);
  }

  friend void set (ZZ_p & x)
// x = 1
  {
    set (x._ZZ_p__rep);
  }

  friend void swap (ZZ_p & x, ZZ_p & y)
// swap x and y
  {
    swap (x._ZZ_p__rep, y._ZZ_p__rep);
  }
  
  // ****** comparison

  friend long IsZero (const ZZ_p & a)
  {
    return IsZero (a._ZZ_p__rep);
  }


  friend long IsOne (const ZZ_p & a)
  {
    return IsOne (a._ZZ_p__rep);
  }

  friend long operator == (const ZZ_p & a, const ZZ_p & b)
    { return a._ZZ_p__rep == b._ZZ_p__rep;}

  friend long operator != (const ZZ_p & a, const ZZ_p & b)
  {
    return !(a == b);
  }

  friend long operator == (const ZZ_p & a, long b);
  friend long operator == (long a, const ZZ_p & b) { return b == a; }

  friend long operator != (const ZZ_p & a, long b)
  {
    return !(a == b);
  }
  friend long operator != (long a, const ZZ_p & b)
  {
    return !(a == b);
  }

  // ****** input/output

  friend NTL_SNS ostream & operator << (NTL_SNS ostream & s, const ZZ_p & a)
  {
    return s << a._ZZ_p__rep;
  }

  friend NTL_SNS istream & operator >> (NTL_SNS istream & s, ZZ_p & x);


  // some other friends...not for general use

  friend void BlockConstruct (ZZ_p *, long);
  friend void BlockDestroy (ZZ_p *, long);
};

// ****** addition

inline void add (ZZ_p & x, const ZZ_p & a, const ZZ_p & b) {
// x = a + b
    ZZ_pInfo->add (x, a, b);
}

inline void sub (ZZ_p & x, const ZZ_p & a, const ZZ_p & b) {
// x = a - b
    ZZ_pInfo->sub (x, a, b);
}


inline void negate (ZZ_p & x, const ZZ_p & a) {
// x = -a
    ZZ_pInfo->negate (x, a);
}

// scalar versions

inline void add (ZZ_p & x, const ZZ_p & a, long b) { ZZ_pInfo->add (x, a, b); }
inline void add (ZZ_p & x, long a, const ZZ_p & b) { ZZ_pInfo->add (x, a, b); }
inline void sub (ZZ_p & x, const ZZ_p & a, long b) { ZZ_pInfo->sub (x, a, b); }
inline void sub (ZZ_p & x, long a, const ZZ_p & b) { ZZ_pInfo->sub (x, a, b); }


// ****** multiplication

inline void mul (ZZ_p & x, const ZZ_p & a, const ZZ_p & b) {
// x = a*b
    ZZ_pInfo->mul (x, a, b);
}


inline void sqr (ZZ_p & x, const ZZ_p & a) {
// x = a^2
    ZZ_pInfo->sqr (x, a);
}

inline ZZ_p sqr (const ZZ_p & a) {
    return ZZ_pInfo->sqr (a);
}


// scalar versions

inline void mul (ZZ_p & x, const ZZ_p & a, long b) { ZZ_pInfo->mul (x, a, b); }
inline void mul (ZZ_p & x, long a, const ZZ_p & b) { ZZ_pInfo->mul (x, a, b); }

// ****** division

inline void div (ZZ_p & x, const ZZ_p & a, const ZZ_p & b) { ZZ_pInfo->div (x, a, b); }
// x = a/b
// If b != 0 & b not invertible & DivHandler != 0,
// then DivHandler will be called with the offending b.
// In this case, of course, p is not really prime, and one
// can factor p by taking a gcd with rep(b).
// Otherwise, if b is not invertible, an error occurs.

inline void inv (ZZ_p & x, const ZZ_p & a) { ZZ_pInfo->inv (x, a); }

// x = 1/a
// Error handling is the same as above.

inline ZZ_p inv (const ZZ_p & a) { return ZZ_pInfo->inv (a); }

inline void div (ZZ_p & x, const ZZ_p & a, long b) { ZZ_pInfo->div (x, a, b); }
inline void div (ZZ_p & x, long a, const ZZ_p & b) { ZZ_pInfo->div (x, a, b); }

// operator notation:

inline ZZ_p operator + (const ZZ_p & a, const ZZ_p & b) {
    ZZ_p x; add (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p operator + (const ZZ_p & a, long b) {
    ZZ_p x; add (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p operator + (long a, const ZZ_p & b) {
    ZZ_p x; add (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p & operator += (ZZ_p & x, const ZZ_p & b) {
    add (x, x, b); return x;
}

inline ZZ_p & operator += (ZZ_p & x, long b) {
    add (x, x, b); return x;
}

inline ZZ_p operator - (const ZZ_p & a, const ZZ_p & b) {
    ZZ_p x; sub (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p operator - (const ZZ_p & a, long b) {
    ZZ_p x; sub (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p operator - (long a, const ZZ_p & b) {
    ZZ_p x; sub (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p & operator -= (ZZ_p & x, const ZZ_p & b) {
    sub (x, x, b); return x;
}

inline ZZ_p & operator -= (ZZ_p & x, long b) {
    sub (x, x, b); return x;
}

inline ZZ_p operator *(const ZZ_p & a, const ZZ_p & b) {
    ZZ_p x; mul (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p operator *(const ZZ_p & a, long b) {
    ZZ_p x; mul (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p operator *(long a, const ZZ_p & b) {
    ZZ_p x; mul (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p & operator *= (ZZ_p & x, const ZZ_p & b) {
    mul (x, x, b); return x;
}

inline ZZ_p & operator *= (ZZ_p & x, long b) {
    mul (x, x, b); return x;
}

inline ZZ_p operator / (const ZZ_p & a, const ZZ_p & b) {
    ZZ_p x; div (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p operator / (const ZZ_p & a, long b) {
    ZZ_p x; div (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p operator / (long a, const ZZ_p & b) {
    ZZ_p x; div (x, a, b); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p& operator/=(ZZ_p& x, const ZZ_p& b)
   { div(x, x, b); return x; } 

inline ZZ_p& operator/=(ZZ_p& x, long b)
   { div(x, x, b); return x; } 


inline ZZ_p operator-(const ZZ_p& a) 
   { ZZ_p x; negate(x, a); NTL_OPT_RETURN(ZZ_p, x); }


inline ZZ_p& operator++(ZZ_p& x) { add(x, x, 1); return x; }
inline void operator++(ZZ_p& x, int) { add(x, x, 1); }
inline ZZ_p& operator--(ZZ_p& x) { sub(x, x, 1); return x; }
inline void operator--(ZZ_p& x, int) { sub(x, x, 1); }

// ****** exponentiation

inline void power (ZZ_p & x, const ZZ_p & a, const ZZ & e) { ZZ_pInfo->power (x, a, e); }
inline ZZ_p power (const ZZ_p & a, const ZZ & e) { return ZZ_pInfo->power (a, e); }
inline void power (ZZ_p & x, const ZZ_p & a, long e) { ZZ_pInfo->power (x, a, e); }
inline ZZ_p power (const ZZ_p & a, long e) { return ZZ_pInfo->power (a, e); }

// ****** random numbers

inline void random (ZZ_p & x) {
// x = random element in ZZ_p
    ZZ_pInfo->random (x);
}

inline ZZ_p random_ZZ_p () { return ZZ_pInfo->random_ZZ_p (); }

// Inline definitions for ZZ_pInfoT methods

inline void ZZ_pInfoT::conv (ZZ_p & x, const ZZ & a) {
    rem (x._ZZ_p__rep, a, p);
}

inline void ZZ_pInfoT::add (ZZ_p & x, const ZZ_p & a, const ZZ_p & b) {
// x = a + b
    AddMod (x._ZZ_p__rep, a._ZZ_p__rep, b._ZZ_p__rep, p);
}

inline void ZZ_pInfoT::sub (ZZ_p & x, const ZZ_p & a, const ZZ_p & b) {
// x = a - b
    SubMod (x._ZZ_p__rep, a._ZZ_p__rep, b._ZZ_p__rep, p);
}

inline void ZZ_pInfoT::negate (ZZ_p & x, const ZZ_p & a) {
// x = -a
    NegateMod (x._ZZ_p__rep, a._ZZ_p__rep, p);
}

inline void ZZ_pInfoT::mul (ZZ_p & x, const ZZ_p & a, const ZZ_p & b) {
// x = a*b
    MulMod (x._ZZ_p__rep, a._ZZ_p__rep, b._ZZ_p__rep, p);
}

inline void ZZ_pInfoT::sqr (ZZ_p & x, const ZZ_p & a) {
// x = a^2
    SqrMod (x._ZZ_p__rep, a._ZZ_p__rep, p);
}

inline ZZ_p ZZ_pInfoT::sqr (const ZZ_p & a) {
    ZZ_p x; sqr (x, a); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p ZZ_pInfoT::inv (const ZZ_p & a) {
    ZZ_p x; inv (x, a); NTL_OPT_RETURN (ZZ_p, x);
}

inline void ZZ_pInfoT::power (ZZ_p & x, const ZZ_p & a, const ZZ & e) {
    PowerMod (x._ZZ_p__rep, a._ZZ_p__rep, e, p);
}

inline void ZZ_pInfoT::power (ZZ_p & x, const ZZ_p & a, long e) {
    PowerMod (x._ZZ_p__rep, a._ZZ_p__rep, e, p);
}

inline ZZ_p ZZ_pInfoT::power (const ZZ_p & a, const ZZ & e) {
    ZZ_p x; power (x, a, e); NTL_OPT_RETURN (ZZ_p, x);
}

inline ZZ_p ZZ_pInfoT::power (const ZZ_p & a, long e) {
    ZZ_p x; power (x, a, e); NTL_OPT_RETURN (ZZ_p, x);
}

inline void ZZ_pInfoT::random (ZZ_p & x) {
// x = random element in ZZ_p
    RandomBnd (x._ZZ_p__rep, p);
}

inline ZZ_p ZZ_pInfoT::random_ZZ_p () {
    ZZ_p x; random (x); NTL_OPT_RETURN (ZZ_p, x);
}

NTL_CLOSE_NNS

#endif
