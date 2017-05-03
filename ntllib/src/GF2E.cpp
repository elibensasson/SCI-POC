//#define COUNT_ARITH 1   //Arnon

#include <NTL/GF2E.h>
#include <NTL/tools.h>

#include <NTL/new.h>

#include <omp.h>

NTL_START_IMPL


GF2EInfoT::GF2EInfoT(const GF2X& NewP)
{
   ref_count = 1;

   build(p, NewP);

   if (p.size == 1) {
      if (deg(p) <= NTL_BITS_PER_LONG/2)
         KarCross = 4;
      else
         KarCross = 8;
   }
   else if (p.size == 2)
      KarCross = 8;
   else if (p.size <= 5)
      KarCross = 4;
   else if (p.size == 6)
      KarCross = 3;
   else 
      KarCross = 2;


   if (p.size <= 1) {
      if (deg(p) <= NTL_BITS_PER_LONG/2)
         ModCross = 20;
      else
         ModCross = 40;
   }
   else if (p.size <= 2)
      ModCross = 75;
   else if (p.size <= 4)
      ModCross = 50;
   else
      ModCross = 25;

   if (p.size == 1) {
      if (deg(p) <= NTL_BITS_PER_LONG/2)
         DivCross = 100;
      else
         DivCross = 200;
   }
   else if (p.size == 2)
      DivCross = 400;
   else if (p.size <= 4)
      DivCross = 200;
   else if (p.size == 5)
      DivCross = 150;
   else if (p.size <= 13)
      DivCross = 100;
   else 
      DivCross = 75;

   _card_init = 0;
   _card_exp = p.n;

// This wasn't in original (trunk) version of NTL, may not be necessary //shaul
// Have problems to compile, commented in //shaul
//   NTL::power(_card, 2, p.n); 
/* Copied from original (trunk) version of NTL. May not be thread-safe //shaul
const ZZ& GF2E::cardinality()
{
   if (!GF2EInfo) Error("GF2E::cardinality: undefined modulus");

   if (!GF2EInfo->_card_init) {
      power(GF2EInfo->_card, 2, GF2EInfo->_card_exp);
      GF2EInfo->_card_init = 1;
   }

   return GF2EInfo->_card;
}
*/

}




GF2EInfoT *GF2EInfo = 0; 
typedef GF2EInfoT *GF2EInfoPtr;


static 
void CopyPointer(GF2EInfoPtr& dst, GF2EInfoPtr src)
{
   if (src == dst) return;

   if (dst) {
      dst->ref_count--;

      if (dst->ref_count < 0) 
         Error("internal error: negative GF2EContext ref_count");

      if (dst->ref_count == 0) delete dst;
   }

   if (src) {
      if (src->ref_count == NTL_MAX_LONG) 
         Error("internal error: GF2EContext ref_count overflow");

      src->ref_count++;

   }

   dst = src;
}
   



void GF2E::init(const GF2X& p)
{
   GF2EContext c(p);
   c.restore();
}

/////////////////////// Copied from .h file

void power(GF2E& x, const GF2E& a, const ZZ& e)
   { PowerMod(x._GF2E__rep, a._GF2E__rep, e, GF2E::modulus()); }

void power(GF2E& x, const GF2E& a, long e)
   { PowerMod(x._GF2E__rep, a._GF2E__rep, e, GF2E::modulus()); }

void conv(GF2E& x, const GF2X& a)
// x = (a mod p)
   { rem(x._GF2E__rep, a, GF2E::modulus()); }

void conv(GF2E& x, long a)
   { conv(x._GF2E__rep, a); }

void conv(GF2E& x, GF2 a)
   { conv(x._GF2E__rep, a); }

void conv(GF2E& x, const ZZ& a)
   { conv(x._GF2E__rep, a); }

void random(GF2E& x)
// x = random element in GF2E
   { random(x._GF2E__rep, GF2EInfo->p.n); }

GF2E to_GF2E(const GF2X& a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

GF2E to_GF2E(long a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

GF2E to_GF2E(GF2 a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

GF2E to_GF2E(const ZZ& a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

///////////////////////

GF2EContext::GF2EContext(const GF2X& p)
{
   ptr = NTL_NEW_OP GF2EInfoT(p);
}

GF2EContext::GF2EContext(const GF2EContext& a)
{
   ptr = 0;
   CopyPointer(ptr, a.ptr);
}

GF2EContext& GF2EContext::operator=(const GF2EContext& a)
{
   CopyPointer(ptr, a.ptr);
   return *this;
}


GF2EContext::~GF2EContext()
{
   CopyPointer(ptr, 0);
}

void GF2EContext::save()
{
   CopyPointer(ptr, GF2EInfo);
}

void GF2EContext::restore() const
{
	#pragma omp critical (GF2E_lock)
	{
		CopyPointer(GF2EInfo, ptr);
	}
}



GF2EBak::~GF2EBak()
{
   if (MustRestore)
      CopyPointer(GF2EInfo, ptr);

   CopyPointer(ptr, 0);
}

void GF2EBak::save()
{
   MustRestore = 1;
   CopyPointer(ptr, GF2EInfo);
}



void GF2EBak::restore()
{
	#pragma omp critical(GF2E_lock)
	{
		MustRestore = 0;
		CopyPointer(GF2EInfo, ptr);
	}
}



const GF2E& GF2E::zero()
{
   static GF2E z(GF2E_NoAlloc);
   return z;
}



istream& operator>>(istream& s, GF2E& x)
{
   GF2X y;

   s >> y;
   conv(x, y);

   return s;
}

void GF2EInfoT::div(GF2E& x, const GF2E& a, const GF2E& b)
{
	#if COUNT_ARITH==1
	//cout << "/"; //Arnon
	OpsCounter::increaseDiv();
	#endif
   GF2E t;

   inv(t, b);
   mul(x, a, t);
	#if COUNT_ARITH==1
	//cout << "/end"; //Arnon
	#endif
}

void GF2EInfoT::div(GF2E& x, GF2 a, const GF2E& b)
{
   inv(x, b);
   mul(x, x, a);
}

void GF2EInfoT::div(GF2E& x, long a, const GF2E& b)
{
   inv(x, b);
   mul(x, x, a);
}


void GF2EInfoT::inv(GF2E& x, const GF2E& a)
{
   InvMod(x._GF2E__rep, a._GF2E__rep, p);
}

NTL_END_IMPL
