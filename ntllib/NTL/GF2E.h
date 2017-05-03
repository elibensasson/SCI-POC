

#ifndef NTL_GF2E__H
#define NTL_GF2E__H

#include <NTL/GF2X.h>

NTL_OPEN_NNS

class GF2E;

class GF2EInfoT {
private:
   GF2EInfoT();                       // disabled
   GF2EInfoT(const GF2EInfoT&);   // disabled
   void operator=(const GF2EInfoT&);  // disabled
public:
   long ref_count;

   GF2EInfoT(const GF2X& NewP);
   ~GF2EInfoT() { }

   GF2XModulus p;

   long KarCross;
   long ModCross;
   long DivCross;

ZZ   _card;
   long _card_init;
   long _card_exp;   

// Operations on elements of this field
// Note: Not all of these actually need to have operations as methods of
// the field class. In particular, the add operations do not appear to depend
// on the field at all. However, it would place an undue burden on the
// programmer to have to figure out which operations are available as methods
// of the field and which are not. It also aids the concept of encapsulation
// to allow programmers to proceed without having that knowledge, so all of
// the basic operations are present here.

void add(GF2E& x, const GF2E& a, const GF2E& b);
void add(GF2E& x, const GF2E& a, GF2 b);
void add(GF2E& x, const GF2E& a, long b);
inline void add(GF2E& x, GF2 a, const GF2E& b)  { add(x, b, a); }
inline void add(GF2E& x, long a, const GF2E& b)  { add(x, b, a); }

inline void sub(GF2E& x, const GF2E& a, const GF2E& b) { add(x, a, b); }
inline void sub(GF2E& x, const GF2E& a, GF2 b) { add(x, a, b); }
inline void sub(GF2E& x, const GF2E& a, long b) { add(x, a, b); }
inline void sub(GF2E& x, GF2 a, const GF2E& b) { add(x, a, b); }
inline void sub(GF2E& x, long a, const GF2E& b) { add(x, a, b); }

void negate(GF2E& x, const GF2E& a);

void mul(GF2E& x, const GF2E& a, const GF2E& b);
void sqr(GF2E& x, const GF2E& a);
GF2E sqr(const GF2E& a);
void mul(GF2E& x, const GF2E& a, GF2 b);
void mul(GF2E& x, const GF2E& a, long b);
inline void mul(GF2E& x, GF2 a, const GF2E& b) { mul(x, b, a); }
inline void mul(GF2E& x, long a, const GF2E& b) { mul(x, b, a); }

void div(GF2E& x, const GF2E& a, const GF2E& b);
void inv(GF2E& x, const GF2E& a);
GF2E inv(const GF2E& a);
void div(GF2E& x, const GF2E& a, GF2 b);
void div(GF2E& x, const GF2E& a, long b);
void div(GF2E& x, GF2 a, const GF2E& b);
void div(GF2E& x, long a, const GF2E& b);

void power(GF2E& x, const GF2E& a, const ZZ& e);
GF2E power(const GF2E& a, const ZZ& e);
void power(GF2E& x, const GF2E& a, long e);
GF2E power(const GF2E& a, long e);

void conv(GF2E& x, const GF2X& a);
void conv(GF2E& x, long a);
void conv(GF2E& x, GF2 a);
void conv(GF2E& x, const ZZ& a);
GF2E to_GF2E(const GF2X& a);
GF2E to_GF2E(long a);
GF2E to_GF2E(GF2 a);
GF2E to_GF2E(const ZZ& a);

void trace(GF2& x, const GF2E& a);
GF2 trace(const GF2E& a);

void random(GF2E& x);
GF2E random_GF2E();
};

extern GF2EInfoT *GF2EInfo; // info for current modulus, initially null




class GF2EContext {
private:
GF2EInfoT *ptr;

public:
void save();
void restore() const;

GF2EContext() { ptr = 0; }
GF2EContext(const GF2X& p);

GF2EContext(const GF2EContext&); 


GF2EContext& operator=(const GF2EContext&); 


~GF2EContext();


};


class GF2EBak {
private:
long MustRestore;
GF2EInfoT *ptr;

GF2EBak(const GF2EBak&); // disabled
void operator=(const GF2EBak&); // disabled

public:
void save();
void restore();

GF2EBak() { MustRestore = 0; ptr = 0; }

~GF2EBak();


};



struct GF2E_NoAlloc_type { GF2E_NoAlloc_type() { } };
const GF2E_NoAlloc_type GF2E_NoAlloc = GF2E_NoAlloc_type();



class GF2E {


private:

friend class GF2EInfoT;

// static data


public:

GF2X _GF2E__rep;

// ****** constructors and assignment

GF2E() { _GF2E__rep.xrep.SetMaxLength(GF2E::WordLength()); }

GF2E(GF2E& x, INIT_TRANS_TYPE) : _GF2E__rep(x._GF2E__rep, INIT_TRANS) { }

GF2E(const GF2E& a)  
   { _GF2E__rep.xrep.SetMaxLength(GF2E::WordLength()); _GF2E__rep = a._GF2E__rep; }

GF2E(GF2E_NoAlloc_type) { }  // allocates no space

~GF2E() { } 

GF2E& operator=(const GF2E& a) { _GF2E__rep = a._GF2E__rep; return *this; }


// You can always access the _GF2E__representation directly...if you dare.
GF2X& LoopHole() { return _GF2E__rep; }

static long WordLength() { return GF2EInfo->p.WordLength(); }

static long storage() { return WV_storage(GF2E::WordLength()); }
static const GF2XModulus& modulus() { return GF2EInfo->p; }

static long KarCross() { return GF2EInfo->KarCross; }
static long ModCross() { return GF2EInfo->ModCross; }
static long DivCross() { return GF2EInfo->DivCross; }

static long degree() { return GF2EInfo->p.n; }

static const GF2E& zero();

static const ZZ& cardinality() { return GF2EInfo->_card; }

static void init(const GF2X& NewP);

// read-only access to GF2E representation
friend const GF2X& rep(const GF2E& a) { return a._GF2E__rep; }

friend void clear(GF2E& x)
// x = 0
   { clear(x._GF2E__rep); }

friend void set(GF2E& x)
// x = 1
   { set(x._GF2E__rep); }

friend void swap(GF2E& x, GF2E& y)
// swap x and y

   { swap(x._GF2E__rep, y._GF2E__rep); }

// ****** addition

friend void add(GF2E& x, const GF2E& a, const GF2E& b)
   { add(x._GF2E__rep, a._GF2E__rep, b._GF2E__rep); }

friend void add(GF2E& x, const GF2E& a, GF2 b)
   { add(x._GF2E__rep, a._GF2E__rep, b); }

friend void add(GF2E& x, const GF2E& a, long b)
   { add(x._GF2E__rep, a._GF2E__rep, b); }

friend void add(GF2E& x, GF2 a, const GF2E& b)  { add(x, b, a); }
friend void add(GF2E& x, long a, const GF2E& b)  { add(x, b, a); }

friend void sub(GF2E& x, const GF2E& a, const GF2E& b) { add(x, a, b); }
friend void sub(GF2E& x, const GF2E& a, GF2 b) { add(x, a, b); }
friend void sub(GF2E& x, const GF2E& a, long b) { add(x, a, b); }
friend void sub(GF2E& x, GF2 a, const GF2E& b) { add(x, a, b); }
friend void sub(GF2E& x, long a, const GF2E& b) { add(x, a, b); }

friend void negate(GF2E& x, const GF2E& a) { x = a; }


friend GF2E operator+(const GF2E& a, const GF2E& b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator+(const GF2E& a, GF2 b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator+(const GF2E& a, long b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator+(GF2 a, const GF2E& b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator+(long a, const GF2E& b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }


friend GF2E operator-(const GF2E& a, const GF2E& b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(const GF2E& a, GF2 b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(const GF2E& a, long b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(GF2 a, const GF2E& b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(long a, const GF2E& b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(const GF2E& a)
   { GF2E x; negate(x, a); NTL_OPT_RETURN(GF2E, x); } 


friend GF2E& operator+=(GF2E& x, const GF2E& b)
   { add(x, x, b); return x; }

friend GF2E& operator+=(GF2E& x, GF2 b)
   { add(x, x, b); return x; }

friend GF2E& operator+=(GF2E& x, long b)
   { add(x, x, b); return x; }


friend GF2E& operator-=(GF2E& x, const GF2E& b)
   { sub(x, x, b); return x; }

friend GF2E& operator-=(GF2E& x, GF2 b)
   { sub(x, x, b); return x; }

friend GF2E& operator-=(GF2E& x, long b)
   { sub(x, x, b); return x; }


friend GF2E& operator++(GF2E& x) { add(x, x, 1); return x; }

friend void operator++(GF2E& x, int) { add(x, x, 1); }

friend GF2E& operator--(GF2E& x) { sub(x, x, 1); return x; }

friend void operator--(GF2E& x, int) { sub(x, x, 1); }



// ****** multiplication

friend void mul(GF2E& x, const GF2E& a, const GF2E& b)
// x = a*b

   { MulMod(x._GF2E__rep, a._GF2E__rep, b._GF2E__rep, GF2E::modulus()); }


friend void sqr(GF2E& x, const GF2E& a)
// x = a^2

   { SqrMod(x._GF2E__rep, a._GF2E__rep, GF2E::modulus()); }

friend GF2E sqr(const GF2E& a)
   { GF2E x; sqr(x, a); NTL_OPT_RETURN(GF2E, x); }

friend void mul(GF2E& x, const GF2E& a, GF2 b)
   { mul(x._GF2E__rep, a._GF2E__rep, b); }

friend void mul(GF2E& x, const GF2E& a, long b)
   { mul(x._GF2E__rep, a._GF2E__rep, b); }

friend void mul(GF2E& x, GF2 a, const GF2E& b) { mul(x, b, a); }
friend void mul(GF2E& x, long a, const GF2E& b) { mul(x, b, a); }



friend GF2E operator*(const GF2E& a, const GF2E& b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator*(const GF2E& a, GF2 b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator*(const GF2E& a, long b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator*(GF2 a, const GF2E& b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator*(long a, const GF2E& b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }


friend GF2E& operator*=(GF2E& x, const GF2E& b)
   { mul(x, x, b); return x; }

friend GF2E& operator*=(GF2E& x, GF2 b)
   { mul(x, x, b); return x; }

friend GF2E& operator*=(GF2E& x, long b)
   { mul(x, x, b); return x; }



// ****** division



friend void div(GF2E& x, const GF2E& a, const GF2E& b) { GF2EInfo->div (x, a, b); }
friend void inv(GF2E& x, const GF2E& a) { GF2EInfo->inv (x, a); }

friend GF2E inv(const GF2E& a)
   { GF2E x; inv(x, a); NTL_OPT_RETURN(GF2E, x); }

friend void div(GF2E& x, const GF2E& a, GF2 b)
   { div(x._GF2E__rep, a._GF2E__rep, b); } 

friend void div(GF2E& x, const GF2E& a, long b)
   { div(x._GF2E__rep, a._GF2E__rep, b); } 

friend void div(GF2E& x, GF2 a, const GF2E& b) { GF2EInfo->div (x, a, b); }
friend void div(GF2E& x, long a, const GF2E& b) { GF2EInfo->div (x, a, b); }


friend GF2E operator/(const GF2E& a, const GF2E& b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator/(const GF2E& a, GF2 b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator/(const GF2E& a, long b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator/(GF2 a, const GF2E& b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator/(long a, const GF2E& b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }


friend GF2E& operator/=(GF2E& x, const GF2E& b)
   { div(x, x, b); return x; }

friend GF2E& operator/=(GF2E& x, GF2 b)
   { div(x, x, b); return x; }

friend GF2E& operator/=(GF2E& x, long b)
   { div(x, x, b); return x; }


// ****** exponentiation

friend void power(GF2E& x, const GF2E& a, const ZZ& e);

friend GF2E power(const GF2E& a, const ZZ& e)
   { GF2E x; power(x, a, e); NTL_OPT_RETURN(GF2E, x); }

friend void power(GF2E& x, const GF2E& a, long e);

friend GF2E power(const GF2E& a, long e)
   { GF2E x; power(x, a, e); NTL_OPT_RETURN(GF2E, x); }


// ****** conversion

friend void conv(GF2E& x, const GF2X& a);
// x = (a mod p)

//   { rem(x._GF2E__rep, a, GF2E::modulus()); }

friend void conv(GF2E& x, long a);
   //{ conv(x._GF2E__rep, a); }

friend void conv(GF2E& x, GF2 a);
//   { conv(x._GF2E__rep, a); }

friend void conv(GF2E& x, const ZZ& a);
//   { conv(x._GF2E__rep, a); }

friend GF2E to_GF2E(const GF2X& a);
//   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

friend GF2E to_GF2E(long a);
//   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

friend GF2E to_GF2E(GF2 a);
//   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

friend GF2E to_GF2E(const ZZ& a);
//   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

// ****** trace

friend void trace(GF2& x, const GF2E& a)
   { TraceMod(x, a._GF2E__rep, GF2E::modulus()); }
friend GF2 trace(const GF2E& a)
   { return TraceMod(a._GF2E__rep, GF2E::modulus()); }



// ****** random numbers

friend void random(GF2E& x);
// x = random element in GF2E

//   { random(x._GF2E__rep, GF2EInfo->p.n); }

friend GF2E random_GF2E()
   { GF2E x; random(x); NTL_OPT_RETURN(GF2E, x); }



// ****** comparison

friend long IsZero(const GF2E& a)
   { return IsZero(a._GF2E__rep); }

friend long IsOne(const GF2E& a)
   { return IsOne(a._GF2E__rep); }

friend long operator==(const GF2E& a, const GF2E& b)
   { return a._GF2E__rep == b._GF2E__rep; }

friend long operator==(const GF2E& a, GF2 b)
   { return a._GF2E__rep == b; }

friend long operator==(const GF2E& a, long b)
   { return a._GF2E__rep == b; }

friend long operator==(const GF2 a, const GF2E& b)
   { return a == b._GF2E__rep; }

friend long operator==(const long a, const GF2E& b)
   { return a == b._GF2E__rep; }


friend long operator!=(const GF2E& a, const GF2E& b) { return !(a == b); }
friend long operator!=(const GF2E& a, GF2 b) { return !(a == b); }
friend long operator!=(const GF2E& a, long b) { return !(a == b); }
friend long operator!=(GF2 a, const GF2E& b) { return !(a == b); }
friend long operator!=(long a, const GF2E& b) { return !(a == b); }


// ****** input/output

friend NTL_SNS ostream& operator<<(NTL_SNS ostream& s, const GF2E& a)
   { return s << a._GF2E__rep; }
   
friend NTL_SNS istream& operator>>(NTL_SNS istream& s, GF2E& x);


GF2E& operator=(long a) { conv(*this, a); return *this; }
GF2E& operator=(GF2 a) { conv(*this, a); return *this; }


// some other friends...not for general use

friend void BlockConstruct(GF2E*,long);
friend void BlockDestroy(GF2E*,long);

};


// Methods of GF2EInfoT that correspond to operations above

// ****** addition

inline void GF2EInfoT::add(GF2E& x, const GF2E& a, const GF2E& b)
   { NTL::add(x._GF2E__rep, a._GF2E__rep, b._GF2E__rep); }

inline void GF2EInfoT::add(GF2E& x, const GF2E& a, GF2 b)
   { NTL::add(x._GF2E__rep, a._GF2E__rep, b); }

inline void GF2EInfoT::add(GF2E& x, const GF2E& a, long b)
   { NTL::add(x._GF2E__rep, a._GF2E__rep, b); }

inline void GF2EInfoT::negate(GF2E& x, const GF2E& a) { x = a; }

// ****** multiplication

inline void GF2EInfoT::mul(GF2E& x, const GF2E& a, const GF2E& b)
// x = a*b

   { MulMod(x._GF2E__rep, a._GF2E__rep, b._GF2E__rep, p); }


inline void GF2EInfoT::sqr(GF2E& x, const GF2E& a)
// x = a^2

   { SqrMod(x._GF2E__rep, a._GF2E__rep, p); }

inline GF2E GF2EInfoT::sqr(const GF2E& a)
   { GF2E x; sqr(x, a); NTL_OPT_RETURN(GF2E, x); }

inline void GF2EInfoT::mul(GF2E& x, const GF2E& a, GF2 b)
   { NTL::mul(x._GF2E__rep, a._GF2E__rep, b); }

inline void GF2EInfoT::mul(GF2E& x, const GF2E& a, long b)
   { NTL::mul(x._GF2E__rep, a._GF2E__rep, b); }

// ****** division

inline GF2E GF2EInfoT::inv(const GF2E& a)
   { GF2E x; inv(x, a); NTL_OPT_RETURN(GF2E, x); }

inline void GF2EInfoT::div(GF2E& x, const GF2E& a, GF2 b)
   { NTL::div(x._GF2E__rep, a._GF2E__rep, b); } 

inline void GF2EInfoT::div(GF2E& x, const GF2E& a, long b)
   { NTL::div(x._GF2E__rep, a._GF2E__rep, b); } 

// ****** exponentiation

inline void GF2EInfoT::power(GF2E& x, const GF2E& a, const ZZ& e)
   { PowerMod(x._GF2E__rep, a._GF2E__rep, e, GF2E::modulus()); }

inline GF2E GF2EInfoT::power(const GF2E& a, const ZZ& e)
   { GF2E x; power(x, a, e); NTL_OPT_RETURN(GF2E, x); }

inline void GF2EInfoT::power(GF2E& x, const GF2E& a, long e)
   { PowerMod(x._GF2E__rep, a._GF2E__rep, e, GF2E::modulus()); }

inline GF2E GF2EInfoT::power(const GF2E& a, long e)
   { GF2E x; power(x, a, e); NTL_OPT_RETURN(GF2E, x); }

// ****** conversion

inline void GF2EInfoT::conv(GF2E& x, const GF2X& a)
// x = (a mod p)

   { NTL::rem(x._GF2E__rep, a, p); }

inline void GF2EInfoT::conv(GF2E& x, long a)
   { NTL::conv(x._GF2E__rep, a); }

inline void GF2EInfoT::conv(GF2E& x, GF2 a)
   { NTL::conv(x._GF2E__rep, a); }

inline void GF2EInfoT::conv(GF2E& x, const ZZ& a)
   { NTL::conv(x._GF2E__rep, a); }

inline GF2E GF2EInfoT::to_GF2E(const GF2X& a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

inline GF2E GF2EInfoT::to_GF2E(long a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

inline GF2E GF2EInfoT::to_GF2E(GF2 a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

inline GF2E GF2EInfoT::to_GF2E(const ZZ& a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

// ****** trace

inline void GF2EInfoT::trace(GF2& x, const GF2E& a)
   { TraceMod(x, a._GF2E__rep, p); }
   
inline GF2 GF2EInfoT::trace(const GF2E& a)
   { return TraceMod(a._GF2E__rep, p); }



// ****** random numbers

inline void GF2EInfoT::random(GF2E& x)
// x = random element in GF2E

   { NTL::random(x._GF2E__rep, p.n); }

inline GF2E GF2EInfoT::random_GF2E()
   { GF2E x; random(x); NTL_OPT_RETURN(GF2E, x); }

GF2E to_GF2E(const GF2X& a);
GF2E to_GF2E(long a);
GF2E to_GF2E(GF2 a);
GF2E to_GF2E(const ZZ& a);

NTL_CLOSE_NNS

#endif
