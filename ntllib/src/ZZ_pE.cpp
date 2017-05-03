

#include <NTL/ZZ_pE.h>

#include <NTL/new.h>

#include <omp.h>

NTL_START_IMPL

ZZ_pEInfoT::ZZ_pEInfoT(const ZZ_pX& NewP)
{
   ref_count = 1;

   build(p, NewP);

   _card_init = 0;
   _card_base = ZZ_p::modulus();
   _card_exp = deg(NewP);

   NTL::power(_card, ZZ_p::modulus(), deg(NewP));
}
//This may need to be removed (thread-safety) //shaul
//const ZZ& ZZ_pE::cardinality()
//{
//   if (!ZZ_pEInfo) Error("ZZ_pE::cardinality: undefined modulus");
//
//   if (!ZZ_pEInfo->_card_init) {
//      power(ZZ_pEInfo->_card, ZZ_pEInfo->_card_base, ZZ_pEInfo->_card_exp);
//      ZZ_pEInfo->_card_init = 1;
//   }
//
//   return ZZ_pEInfo->_card;
//}




ZZ_pEInfoT *ZZ_pEInfo = 0; 
typedef ZZ_pEInfoT *ZZ_pEInfoPtr;


static 
void CopyPointer(ZZ_pEInfoPtr& dst, ZZ_pEInfoPtr src)
{
   if (src == dst) return;

   if (dst) {
      dst->ref_count--;

      if (dst->ref_count < 0) 
         Error("internal error: negative ZZ_pEContext ref_count");

      if (dst->ref_count == 0) delete dst;
   }

   if (src) {
      if (src->ref_count == NTL_MAX_LONG) 
         Error("internal error: ZZ_pEContext ref_count overflow");

      src->ref_count++;
   }

   dst = src;
}
   



void ZZ_pE::init(const ZZ_pX& p)
{
   ZZ_pEContext c(p);
   c.restore();
}


ZZ_pEContext::ZZ_pEContext(const ZZ_pX& p)
{
   ptr = NTL_NEW_OP ZZ_pEInfoT(p);
}

ZZ_pEContext::ZZ_pEContext(const ZZ_pEContext& a)
{
   ptr = 0;
   CopyPointer(ptr, a.ptr);
}

ZZ_pEContext& ZZ_pEContext::operator=(const ZZ_pEContext& a)
{
   CopyPointer(ptr, a.ptr);
   return *this;
}


ZZ_pEContext::~ZZ_pEContext()
{
   CopyPointer(ptr, 0);
}

void ZZ_pEContext::save()
{
   CopyPointer(ptr, ZZ_pEInfo);
}

void ZZ_pEContext::restore() const
{
	#pragma omp critical (ZZ_pE_lock)
	{
	   CopyPointer(ZZ_pEInfo, ptr);
	}
}


ZZ_pEBak::~ZZ_pEBak()
{
   if (MustRestore)
      CopyPointer(ZZ_pEInfo, ptr);

   CopyPointer(ptr, 0);
}

void ZZ_pEBak::save()
{
   MustRestore = 1;
   CopyPointer(ptr, ZZ_pEInfo);
}



void ZZ_pEBak::restore()
{
	#pragma omp critical (ZZ_pE_lock)
	{
	   MustRestore = 0;
	   CopyPointer(ZZ_pEInfo, ptr);
	}
}


const ZZ_pE& ZZ_pE::zero()
{
   static ZZ_pE z(ZZ_pE_NoAlloc);
   return z;
}


ZZ_pE::ZZ_pE()
{
   _ZZ_pE__rep.rep.SetMaxLength(ZZ_pE::degree());
}



istream& operator>>(istream& s, ZZ_pE& x)
{
   ZZ_pX y;

   s >> y;
   conv(x, y);

   return s;
}

void ZZ_pEInfoT::div(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
{
   ZZ_pE t;

   inv(t, b);
   mul(x, a, t);
}

void ZZ_pEInfoT::div(ZZ_pE& x, const ZZ_pE& a, long b)
{
   NTL_ZZ_pRegister(B);
   B = b;
   NTL::inv(B, B);
   mul(x, a, B);
}

void ZZ_pEInfoT::div(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
{
   NTL_ZZ_pRegister(B);
   B = b;
   NTL::inv(B, B);
   mul(x, a, B);
}

void ZZ_pEInfoT::div(ZZ_pE& x, long a, const ZZ_pE& b)
{
   ZZ_pE t;
   inv(t, b);
   mul(x, a, t);
}

void ZZ_pEInfoT::div(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
{
   ZZ_pE t;
   inv(t, b);
   mul(x, a, t);
}



void ZZ_pEInfoT::inv(ZZ_pE& x, const ZZ_pE& a)
{
   InvMod(x._ZZ_pE__rep, a._ZZ_pE__rep, p);
}

NTL_END_IMPL
