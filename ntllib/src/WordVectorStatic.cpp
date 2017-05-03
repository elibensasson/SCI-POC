
#include <NTL/WordVectorStatic.h>

#include <NTL/new.h>
#include <stdio.h>

NTL_START_IMPL


///Sets the vector's length. The function is called in the static case only.
void WordVectorStatic::DoSetLength(long n)   
{   
	long m;  

	if (n < 0) {  
		Error("negative length in vector::SetLength");  
	}  

	if (NTL_OVERFLOW(n, NTL_BITS_PER_LONG, 0)) 
		Error("length too big in vector::SetLength");

	if (n == 0) {  
		if (dynamicRep) dynamicRep[-1] = 0;  
		return;  
	}  

	if (!dynamicRep) {  
		m = ((n+NTL_WordVectorMinAlloc-1)/NTL_WordVectorMinAlloc) * NTL_WordVectorMinAlloc; 

		if (NTL_OVERFLOW(m, NTL_BITS_PER_LONG, 0))
			Error("length too big in vector::SetLength");

		_ntl_ulong *p = (_ntl_ulong *) 
			NTL_MALLOC(m, sizeof(_ntl_ulong), 2*sizeof(_ntl_ulong));

		if (!p) {  
			Error("out of memory in SetLength()");  
		}  
		dynamicRep = p+2;

		dynamicRep[-1] = n;
		dynamicRep[-2] = m << 1;

		return;
	}  

	long max_length = (dynamicRep[-2] >> 1);

	if (n <= max_length) {  
		dynamicRep[-1] = n;  
		return;
	}  

	long frozen = (dynamicRep[-2] & 1);

	if (frozen) Error("Cannot grow this WordVector");

	m = max(n, long(NTL_WordVectorExpansionRatio*max_length));

	m = ((m+NTL_WordVectorMinAlloc-1)/NTL_WordVectorMinAlloc)*NTL_WordVectorMinAlloc; 
	_ntl_ulong *p = dynamicRep - 2;

	if (NTL_OVERFLOW(m, NTL_BITS_PER_LONG, 0))
		Error("length too big in vector::SetLength");

	p = (_ntl_ulong *) 
		NTL_REALLOC(p, m, sizeof(_ntl_ulong), 2*sizeof(_ntl_ulong)); 
	if (!p) {  
		Error("out of memory in SetLength()");  
	}  
	dynamicRep = p+2;

	dynamicRep[-1] = n;
	dynamicRep[-2] = m << 1;
}  

///Sets the vector's length. In the static case - nothing to do but parameters check...
void WordVectorStatic::SetLength(long n)
{

	///Not in static mode - Regular WordVector behavior.
	if (!staticMode) {
		_ntl_ulong *x = dynamicRep;
		if (x && long(x[-2] >> 1) >= n && n >= 0)
			x[-1] = n;
		else
			DoSetLength(n);
		return;
	}

	//assert(!dynamicRep);
	//if (dynamicRep)
	//	free(dynamicRep-2);

	///In static mode, and staying in static mode.
	if (n <= STATIC_VECTOR_LENGTH)
		return;

	///In static mode, but switching to dynamic.
	//convertToDynamic();	
	dynamicRep = 0;
	DoSetLength(n);	
	for (int i=0; i<STATIC_VECTOR_LENGTH; i++)
		dynamicRep[i] = statRep[i];
	staticMode = 0;
}

///Sets the vector's length. In the static case - nothing to do but parameters check...
void WordVectorStatic::SetMaxLength(long n) 
{ 
	///Not in static mode - Regular WordVector behavior.
	if (!staticMode) {
		long OldLength = length(); 
		DoSetLength(n); 
		if (dynamicRep) dynamicRep[-1] = OldLength;
		return;
	}

	///In static mode, and staying in static mode.
	if (n <= STATIC_VECTOR_LENGTH)
		return;

	///In static mode, but switching to dynamic.
	//convertToDynamic();
	dynamicRep = 0;
	DoSetLength(n);//this->SetMaxLength(n);	///This time will enter the dynamic case.
	for (int i=0; i<STATIC_VECTOR_LENGTH; i++)
		dynamicRep[i] = statRep[i];
	staticMode = 0;
} 

///Resets the vector by switching to static mode and initializing elements to 0.
void WordVectorStatic::ZeroLength() { 

	for (int i=0; i<STATIC_VECTOR_LENGTH; i++)
		statRep[i] = 0;	   

	if (!staticMode) {		
		if (dynamicRep) {
			//printf("**1**");
			//if (dynamicRep[-2] & 1) Error("Cannot free this WordVector");
			free(dynamicRep-2);
			//printf("2");
		}
		staticMode = 1;
		dynamicRep = 0;
	}
	else
		dynamicRep = 0;

	//if (dynamicRep) dynamicRep[-1] = 0; 
}

//Assignment operator for static word vectors:  
WordVectorStatic& WordVectorStatic::operator=(const WordVectorStatic& a)  
{  

	if (this == &a) return *this;  
		
	//Static case - Easy copy
	if (a.staticMode) {
		this->staticMode = 1;
		dynamicRep = 0;
		for (int i=0; i<STATIC_VECTOR_LENGTH; i++)
			statRep[i] = a.statRep[i];
	}

	//a is dynamic but can fit in into the static case
	else if (a.length() <= STATIC_VECTOR_LENGTH) {
		this->staticMode = 1;
		dynamicRep = 0;
		for (int i=0; i<a.length(); i++)
			statRep[i] = a.dynamicRep[i];
	}

	//Dynamic case - allocations etc.
	else {

		//if (staticMode) 
		//	convertToDynamic();
		
		//staticMode = 0;
		long i, n;  
		_ntl_ulong *p;  
		const _ntl_ulong *ap;  

		n = a.length();  
		ap = a.elts();  

		if (!staticMode) {
			if (dynamicRep)
				free(dynamicRep-2);
			dynamicRep = 0;
		}
		staticMode = 0;
		SetLength(n); 	
		p = elts();  

		for (i = 0; i < n; i++)  
			p[i] = ap[i];  

	}

	return *this;
}  

///Class destructor. In the static case - does nothing.  
WordVectorStatic::~WordVectorStatic()  
{
	if (!staticMode) {
		if (!dynamicRep) return;  
		if (dynamicRep[-2] & 1) Error("Cannot free this WordVector");
		free(dynamicRep-2);
	}
}  

///Free the memory. In the static case - does nothing.  
void WordVectorStatic::kill()  
{  
	if (!staticMode) {
		if (!dynamicRep) return;  
		if (dynamicRep[-2] & 1) 
			Error("Cannot free this WordVector");
		free(dynamicRep-2);
		dynamicRep = 0; 
	}
	dynamicRep = 0;
}  

///Inform a range error
void WordVectorStatic::RangeError(long i) const  
{  
	if (staticMode)
		cerr << "index out of range in STATIC vector: ";  
	else
		cerr << "index out of range in DYNAMIC vector: "; 
	cerr << i;  
	if (!dynamicRep)  
		cerr << "(0)";  
	else  
		cerr << "(" << STATIC_VECTOR_LENGTH << ")";  
	Error("");  
}  

/***********************************/
/******** Static -> Dynamic ********/
/***********************************/

/** Converts the object to a dynamic one, keeping the value as the original. */
void WordVectorStatic::convertToDynamic() {
	
	assert(!dynamicRep);
	assert(staticMode);
	
	long oldLen = this->length();
	DoSetLength(oldLen); 
	for (int i=0; i<oldLen; i++)
		dynamicRep[i] = statRep[i];	

	staticMode = 0;	
}

//Copy swap between two static word vectors, using a temporary vector.
void CopySwap(WordVectorStatic& x, WordVectorStatic& y)
{
	static WordVectorStatic t;
	t = x;
	x = y;
	y = t;
}

//Swap between two static word vectors, using a temporary vector.
void WordVectorStatic::swap_impl(WordVectorStatic& x, WordVectorStatic& y)  
{  
	assert(x.staticMode == y.staticMode);	//If fails, need to handle this case too.

	//Static case:
	if (x.staticMode) {
		_ntl_ulong t[STATIC_VECTOR_LENGTH];  
		for (int i=0; i<STATIC_VECTOR_LENGTH; i++) {
			t[i] = x.statRep[i];
			x.statRep[i] = y.statRep[i];
			y.statRep[i] = t[i];
		}
	}

	//Dynamic case:
	else {
		if ((x.dynamicRep && (x.dynamicRep[-2] & 1)) ||
			(y.dynamicRep && (y.dynamicRep[-2] & 1))) {
				CopySwap(x, y);
				return;
		}

		_ntl_ulong* t; 
		t = x.dynamicRep;  
		x.dynamicRep = y.dynamicRep;  
		y.dynamicRep = t;  
	}
} 

//Append another word to the object. In the static case - Illegal!
void WordVectorStatic::append_impl(WordVectorStatic& v, _ntl_ulong a)  
{
	if (v.staticMode)
		Error("Illegal call to append_impl in WordVectorStatic");
	
	long l = v.length();
	v.SetLength(l+1);  
	v[l] = a;  
}  

//Append another word to the object. In the static case - Illegal!
void WordVectorStatic::append_impl(WordVectorStatic& v, const WordVectorStatic& w)  
{  
	if (v.staticMode)
		Error("Illegal call to append_impl in WordVectorStatic");
	
	long l = v.length();  
	long m = w.length();  
	long i;  
	v.SetLength(l+m);  
	for (i = 0; i < m; i++)  
	v[l+i] = w[i]; 
}

//User input function. Did not change the implementation since the function is not used by the PCPCD code.
istream & operator>>(istream& s, WordVectorStatic& a)   
{   
	WordVectorStatic ibuf;  
	long c;   
	long n;   
	if (!s) Error("bad vector input"); 

	c = s.peek();  
	while (IsWhiteSpace(c)) {  
		s.get();  
		c = s.peek();  
	}  
	if (c != '[') {  
		Error("bad vector input");  
	}  

	n = 0;   
	ibuf.SetLength(0);  

	s.get();  
	c = s.peek();  
	while (IsWhiteSpace(c)) {  
		s.get();  
		c = s.peek();  
	}  
	while (c != ']' && c != EOF) {   
		if (n % NTL_WordVectorInputBlock == 0) ibuf.SetMaxLength(n + NTL_WordVectorInputBlock); 
		n++;   
		ibuf.SetLength(n);   
		if (!(s >> ibuf[n-1])) Error("bad vector input");   
		c = s.peek();  
		while (IsWhiteSpace(c)) {  
			s.get();  
			c = s.peek();  
		}  
	}   
	if (c == EOF) Error("bad vector input");  
	s.get(); 

	a = ibuf; 
	return s;   
}    

//User output function. Prints the word vector to screen.
ostream& operator<<(ostream& s, const WordVectorStatic& a)   
{   
	if (a.staticMode)	s << "Static: ";
	else				s << "Dynamic: ";

	long i, n;   

	n = a.length();  

	s << '[';   

	for (i = 0; i < n; i++) {   
		s << a[i];   
		if (i < n-1) s << " ";   
	}   

	s << ']';   

	return s;   
}   

//Comparison operator. Did not change implementation.
long operator==(const WordVectorStatic& a, const WordVectorStatic& b) 
{  
	long n = a.length();  
	if (b.length() != n) return 0;  
	const _ntl_ulong* ap = a.elts(); 
	const _ntl_ulong* bp = b.elts(); 
	long i;  
	for (i = 0; i < n; i++) if (ap[i] != bp[i]) return 0;  
	return 1;  
} 

//Difference operator. Did not change implementation.
long operator!=(const WordVectorStatic& a, const WordVectorStatic& b) 
{  return !(a == b); }


//inner product between static word vectors. Did not change implementation.
long InnerProduct(const WordVectorStatic& a, const WordVectorStatic& b)
{
	long n = min(a.length(), b.length());
	const _ntl_ulong *ap = a.elts();
	const _ntl_ulong *bp = b.elts();

	_ntl_ulong acc;
	long i;

	acc = 0;
	for (i = 0; i < n; i++)
		acc ^= ap[i] & bp[i];

#if (NTL_BITS_PER_LONG == 32)
	acc ^= acc >> 16;
	acc ^= acc >> 8;
	acc ^= acc >> 4;
	acc ^= acc >> 2;
	acc ^= acc >> 1;
	acc &= 1;
#elif (NTL_BITS_PER_LONG == 64)
	acc ^= acc >> 32;
	acc ^= acc >> 16;
	acc ^= acc >> 8;
	acc ^= acc >> 4;
	acc ^= acc >> 2;
	acc ^= acc >> 1;
	acc &= 1;
#else
	_ntl_ulong t = acc;
	while (t) {
		t = t >> 8;
		acc ^= t;
	}

	acc ^= acc >> 4;
	acc ^= acc >> 2;
	acc ^= acc >> 1;
	acc &= 1;
#endif

	return long(acc);
}

//c = c + (a << n). Did not change implementation.
//void ShiftAdd(_ntl_ulong *cp, const _ntl_ulong* ap, long sa, long n)
//	// c = c + (a << n)
//{
//	if (sa == 0) return;
//
//	long i;
//
//	long wn = n / NTL_BITS_PER_LONG;
//	long bn = n - wn * NTL_BITS_PER_LONG;
//
//	if (bn == 0) {
//		for (i = sa+wn-1; i >= wn; i--)
//			cp[i] ^= ap[i-wn];
//	}
//	else {
//		_ntl_ulong t = ap[sa-1] >> (NTL_BITS_PER_LONG-bn);
//		if (t) cp[sa+wn] ^= t;
//		for (i = sa+wn-1; i >= wn+1; i--)
//			cp[i] ^= (ap[i-wn] << bn) | (ap[i-wn-1] >> (NTL_BITS_PER_LONG-bn));
//		cp[wn] ^= ap[0] << bn;
//	}
//}

/******************************************************************/
/***** Block Constructors \ Destructors. Used by NTL vectors ******/
/******************************************************************/

///Construct a block of word vectors. In the static case - convert to dynamic.
///n - number of elements in the vector. d - size of element.
long WV_BlockConstructAlloc(WordVectorStatic& x, long d, long n)
{
	return n;
	if (x.staticMode) 
		x.convertToDynamic();
	
	if (n <= 0)		/* check n value */
		Error("block construct: n must be positive");
	
	if (d <= 0)		/* check d value */
		Error("block construct: d must be positive");

	if (NTL_OVERFLOW(d, NTL_BITS_PER_LONG, 0) || 
		NTL_OVERFLOW(d, sizeof(_ntl_ulong), 2*sizeof(_ntl_ulong)))
		Error("block construct: d too large");

	long nwords, nbytes, AllocAmt, m, j; 
	_ntl_ulong *p, *q;

	//nwords = d + 2 + 4;	//Change in WV_BlockConstructSet as well
	nwords = d + 2;
	nbytes = nwords*sizeof(_ntl_ulong);

	AllocAmt = (NTL_MAX_ALLOC_BLOCK - sizeof(_ntl_ulong)) / nbytes;
	if (AllocAmt == 0) AllocAmt = 1;

	if (AllocAmt < n)
		m = AllocAmt;
	else
		m = n;

	//cout << "n=" << n << ", d=" << d << ", nwords=" << nwords << ", nbytes=" << nbytes << ", AllocAmt=" << AllocAmt << endl;

	p = (_ntl_ulong *) NTL_MALLOC(m, nbytes, sizeof(_ntl_ulong));		//malloc(m*nbytes + sizeof(_ntl_ulong))
	if (!p) Error("out of memory in block construct");

	*p = m;

	q = p+3;
	x.dynamicRep = q;

	for (j = 0; j < m; j++) {
		q[-2] = (d << 1) | 1;
		q[-1] = 0;
		//q[d] = q[d+1] = q[d+2] = q[d+3] = 0;
		q += nwords;
	}

	return m;
}

///In the static case - does nothing.
void WV_BlockConstructSet(WordVectorStatic& x, WordVectorStatic& y, long i)
{
	//if (x.staticMode && y.staticMode)
	//	return;
	return ;
	if (x.staticMode)		x.convertToDynamic();
	if (y.staticMode)		y.convertToDynamic();
	
	long d = x.dynamicRep[-2] >> 1;
	//long size = d + 2 + 4;
	long size = d + 2;

	y.dynamicRep = x.dynamicRep + i*size;
}

///Destroy a block of words. In the static case - does nothing.
long WV_BlockDestroy(WordVectorStatic& x)
{
	return 1;
	if (x.staticMode)
		//return STATIC_VECTOR_LENGTH;//(x.length() == 0) ? 1 : x.length();
		x.convertToDynamic();

	long m;
	_ntl_ulong *p;

	p = x.dynamicRep - 3;
	m = (long) *p;
	free(p);
	return m;
}


//long WV_storage(long d)
//{
//	return (d + 2)*sizeof(_ntl_ulong) + sizeof(WordVectorStatic);
//}




NTL_END_IMPL
