/*
 * Element.cpp
 *
 *  Created on: May 31, 2014
 *      Author: matan
 */


/*
 * SSE SUPPORT
 */

#ifdef WIN32
//#include <wmmintrin.h>
#include <emmintrin.h>
#endif	// #ifdef WIN32
#ifdef __GNUC__
#include <x86intrin.h>
#endif	// #ifdef __GNUC__
#include <stdint.h>
#include <stdint.h>
#include <bitset>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Element.h"
#include "Definitions.h"
#include <bitset>
#include <iostream>
FFF::Element c_mulXorAux;
namespace FFF {
#ifdef WIN32
 __m128i mask = {-1,-1,-1,-1,-1,-1,-1,-1};
#endif // #ifdef WIN32
#ifdef __GNUC__
       __m128i mask = { (long long int)0b1111111111111111111111111111111111111111111111111111111111111111, 0 };
#endif // #ifdef __GNUC__


/*
 * Operations counters
 */
//#define OP_COUNTERS

namespace telemetry{

#ifdef OP_COUNTERS
struct{
    size_t mul;
    size_t add;
    size_t inv;
    
    void clear(){
        mul = 0;
        add = 0;
        inv = 0;
    }

    void print(){
       std::cout<<"acamulated number of additions : "<<add<<std::endl;
       std::cout<<"acamulated number of multiplications : "<<mul<<std::endl;
       std::cout<<"acamulated number of inverse : "<<inv<<std::endl;
    }
}op_counters;

inline void addMul(){
    op_counters.mul++;
}
inline void addAdd(){
    op_counters.add++;
}
inline void addInv(){
    op_counters.inv++;
}
void reset(){
    op_counters.clear();
}
void print(){
    op_counters.print();
}
#else
inline void addMul(){};
inline void addAdd(){};
inline void addInv(){};
void reset(){};
void print(){
    std::cout<<"OP_COUNTERS macro is undefined, operations counters printing is disabled"<<std::endl;
}
#endif

} // namespace telemetry

const idx_t Element::irr_poly_index[][max_nonzero_coefs_in_mod]={
			{0,(long long unsigned int)-1,(long long unsigned int)-1},
			{0,2,3,7,32},
			{0,1,3,4,64}
		};

const cell_t Element::irr_poly_e[][Element::element_len]={
		{0},
		{0},
		{27}
};
const len_t Element::mod_len[]= {3, 5, 5};

//--------------------------------------------------------------//
	Element::Element(const Element& e){
		memcpy(this->c,e.c,sizeof(cell_t)*element_len);
	}
	void Element::operator=(const Element& e){
		memcpy(this->c,e.c,sizeof(cell_t)*element_len);
	}

void Element::deg(cell_t * a, int & i, len_t len){
	for(i=len-1;a[i]==0;--i){};
	if(i<0){
		i=0;
		return;
	}
	cell_t t = a[i];
	idx_t offset=0;
	i*=bits_in_cell;
	for(idx_t j = bits_in_cell/2 ; j >0 ; j>>=1)
		if((t>>(offset+j))!=0){
			offset+=j;
		}
	i+=offset;
}
void setZero(cell_t * a,len_t len){
	for(unsigned int i = 0 ; i < len ; ++i)
		a[i]=0;
}
void setUnit(cell_t * a, len_t len){
	setZero(a,len);
	a[0]=1;
}
bool isUnit(cell_t * a, len_t len){
	if(a[0]!=1)
		return false;
	for(idx_t i = 1; i  < len ; ++i)
		if(a[i]!=0)
			return false;
	return true;
}

void swap(cell_t * a , cell_t * b,len_t len){
	for(unsigned int i = 0 ; i <  len ; ++i)
	{
		a[i]^=b[i];
		b[i]^=a[i];
		a[i]^=b[i];
	}
}
void Element::assign_cells(cell_t* a, cell_t* b,len_t len){
	memcpy(a,b,sizeof(cell_t)*len);
}



void Element::generateOrdElement(Element* e){
	memset(e->c,0,sizeof(cell_t)*element_len);
	for(unsigned int i = 0 ; i < mod_len[ord>>5]-1;++i)
		TOGGLE_BIT(e->c,irr_poly_index[ord>>5][i]);
}

inline void Element::clmulXor(const cell_t* a,const  cell_t* b, cell_t* res)
{
	_mm_storeu_si128((__m128i*)res,
			_mm_xor_si128(
					*((__m128i*)a),
					_mm_clmulepi64_si128(
							_mm_loadu_si128((__m128i*)a),
							_mm_loadu_si128((__m128i*)b),
							0
							)
			)
	);
//	register __m128i v1,v2;
//	ui64_m128i(v1,(const unsigned long*)a);
//	ui64_m128i(v2,(const unsigned long*)b);
//	m128i_ui64((unsigned long*)res,_mm_clmulepi64_si128(v1,v2,0));

    //telemetry
    {
        telemetry::addMul();
        telemetry::addAdd();
    }
}
inline void Element::clmul(const cell_t* a,const  cell_t* b, cell_t* res)
{
			_mm_clmulepi64_si128(
					_mm_loadu_si128((__m128i*)a),
					_mm_loadu_si128((__m128i*)b),
					0
			);
    
            //telemetry
            {
                telemetry::addMul();
            }
}
/*
 * CPU Operations
 */
	void Element::c_add(const Element* const a,const Element* const b, Element* const c)
	{
		for(unsigned int i = 0 ; i < element_len ; ++i){
			c->c[i]=(a->c[i]^b->c[i]);
		}
        
        //telemetry
        {
            telemetry::addAdd();
        }
	}
	void Element::c_add(const Element& a,const Element& b, Element& c)
	{
		for(unsigned int i = 0 ; i < element_len ; ++i){
			c.c[i]=(a.c[i]^b.c[i]);
		}
        
        //telemetry
        {
            telemetry::addAdd();
        }
	}
	bool Element::equals(Element& a, Element& b){
		return !memcmp(a.c,b.c,sizeof(cell_t)*element_len);
	}
	void Element::c_mulXor(Element& a,Element& b,Element& c){
		cell_t clmul_res[element_len<<1];
		c_mul(&a,&b,(Element*)clmul_res);
		vecXor(clmul_res,c.c,c.c,element_len);
   
        //telemetry
        {
            telemetry::addMul();
            telemetry::addAdd();
        }
	}
void Element::c_mul(const Element* a, const Element* b, Element* c){
		register __m128i l = _mm_loadu_si128((__m128i*)irr_poly_e[ord>>5]);
		register __m128i t;
        {
        //copy a and b to a location where
        //they can be copied from without accesing
        //forbidden locations
        Element arr[3];
        arr[0] = *a;
        arr[1] = *b;

		t=_mm_clmulepi64_si128(
				_mm_loadu_si128((__m128i*)&(arr[0])),
				_mm_loadu_si128((__m128i*)&(arr[1])),
				0
		);
        }
		t=_mm_xor_si128(_mm_clmulepi64_si128(t,l,1),_mm_and_si128(t,mask));
		t=_mm_xor_si128(_mm_clmulepi64_si128(t,l,1),t);
		memcpy(c,&t,sizeof(cell_t)*element_len);
    
        //telemetry
        {
            telemetry::addMul();
        }
	}
void Element::c_mul(const Element& a,const Element& b, Element& c){
		__m128i l = _mm_loadu_si128((__m128i*)irr_poly_e[ord>>5]);
		__m128i t;
        
        {
        //copy a and b to a location where
        //they can be copied from without accesing
        //forbidden locations
        Element arr[3];
        arr[0] = a;
        arr[1] = b;

		t=_mm_clmulepi64_si128(
				_mm_loadu_si128((__m128i*)&(arr[0])),
				_mm_loadu_si128((__m128i*)&(arr[1])),
				0
		);
        }
		
        t=_mm_xor_si128(_mm_clmulepi64_si128(t,l,1),_mm_and_si128(t,mask));
		t=_mm_xor_si128(_mm_clmulepi64_si128(t,l,1),t);
        memcpy(&c,&t,sizeof(cell_t)*element_len);
    
        //telemetry
        {
            telemetry::addMul();
        }
	}
	void Element::c_sqr(Element* a, Element* c){
		c_mul(a,a,c);
	}
	void Element::c_sqr(Element& a, Element& c){
		c_mul(a,a,c);
	}
	void Element::left_shift(cell_t * v, unsigned int k , unsigned int len){
		 if(k>=(len<<log_bits_in_cell)){
				 setZero(v,len);
				 return;
			 }
			 unsigned int i;// = uint_e_regs[0];
			 unsigned int shift_cells;// = uint_e_regs[1];
			 shift_cells = k>>log_bits_in_cell;
			 unsigned int par;// = uint_e_regs[2];
			 par= k & andMask(log_bits_in_cell);
			 unsigned int rest;// = uint_e_regs[3];
			 rest = bits_in_cell-par;
			 for(i = len - 1 ; i >= shift_cells+1 ; --i){
				 v[i]=v[i-shift_cells]<<par | (v[i-shift_cells-1]>>(rest-1))>>1;
			 }
			 v[shift_cells]=v[0]<<par;
			 for(i=shift_cells;i>0;--i)
				 v[shift_cells-1]=0;
	}
	void Element::egcd(cell_t a[element_len+1], cell_t b[element_len+1], cell_t g[element_len+1]){
		 cell_t u[element_len+1];
		 Element::assign_cells(u,a,element_len+1);//e_regs[0]
		 cell_t v[element_len+1];
		 Element::assign_cells(v,b,element_len+1);//e_regs[1]
		 cell_t g1[element_len+1];//e_regs[2]
		 cell_t g2[element_len+1];//e_regs[3]
		 setZero(g2,element_len+1);
		 setUnit(g1,element_len+1);//g1
		 while(!isUnit(u,element_len+1)){
			 int i; //uint_e_regs[20]
			 int j; //uint_e_regs[21]
			 deg(u,i,element_len+1);
			 deg(v,j,element_len+1);
			 cell_t t[element_len+1];//e_regs[6]
			 if(i<j){
				 swap(u,v,element_len+1);
				 swap(g1,g2,element_len+1);
				 i=j-i;//j
			 } else {
				 i-=j;//j
			 }
			 Element::assign_cells(t,v,element_len+1);//v*z^j
			 left_shift(t,i,element_len+1);
			 Element::vecXor(u,t,element_len+1);//u=u+v*z^j

			 Element::assign_cells(t,g2,element_len+1);//g2*z^j
			 left_shift(t,i,element_len+1);
			 Element::vecXor(g1,t,element_len+1);//g1=g1+g2*z^j

		 }
		Element::assign_cells(g,g1,element_len);
	}

	void Element::c_inv(Element& a, Element& res){
		cell_t c[element_len+1];
		cell_t e[element_len+1];
		assign_cells(e,a.c,element_len);//=e_regs[8]
		e[element_len]=0;
		Element e_ord;
		Element::generateOrdElement(&e_ord);
		assign_cells(c,e_ord.c,element_len);
		c[element_len]=1;
		egcd(e,c,res.c);
    
        //telemetry
        {
            telemetry::addInv();
        }
	}
	void Element::c_inv(Element* a, Element* res){
		cell_t c[element_len+1];
		cell_t e[element_len+1];
		assign_cells(e,a->c,element_len);//=e_regs[8]
		e[element_len]=0;
		Element e_ord;
		Element::generateOrdElement(&e_ord);
		assign_cells(c,e_ord.c,element_len);
		c[element_len]=1;
		egcd(e,c,res->c);
    
        //telemetry
        {
            telemetry::addInv();
        }
	}
	void Element::c_setZero(Element* a){
		memset(a->c,0,sizeof(Element));
	}
	void Element::c_setZero(Element& a){
		memset(a.c,0,sizeof(Element));
	}
	void Element::c_setUnit(Element* a){
		c_setZero(a);
		a->c[0]=1;
	}
	void Element::c_setUnit(Element& a){
		c_setZero(a);
		a.c[0]=1;
	}
	bool Element::c_isUnit(Element* a){
		if(a->c[0]!=1)
			return false;
		for(unsigned int i = 1 ; i < element_len ; ++i)
			if(a->c[i]!=0)
				return false;
		return true;
	}
	bool Element::c_isUnit(Element& a){
		if(a.c[0]!=1)
			return false;
		for(unsigned int i = 1 ; i < element_len ; ++i)
			if(a.c[i]!=0)
				return false;
		return true;
	}

void Element::c_exp(Element a, unsigned long long exp, Element& res){
	Element::c_setUnit(res);
	for(; exp != 0 ; exp>>=1){
		if(1&exp)
			Element::c_mul(res,a,res);
		Element::c_sqr(a,a);
	}
}
void Element::vecXor(cell_t * a, cell_t * b, len_t len){
		for(unsigned int i = 0 ; i <  len ; ++i)
			a[i]^=b[i];
	}
void Element::vecXor(cell_t * a, cell_t * b, cell_t* c, len_t len){
		for(unsigned int i = 0 ; i <  len ; ++i)
			c[i]=(a[i]^b[i]);
	}
void Element::printCells(const cell_t * x, len_t l){
	for(unsigned int i = l-1 ; i < (~0) ; --i){
		std::cout <<std::bitset<bits_in_cell>(x[i]);
		if(i!=0)
			std::cout << ",";
	}

}
void Element::printElement(const Element& x){
	printCells(x.c,element_len);
}
void Element::assign(Element (*a),const Element (*b)){
	memcpy(a,b,sizeof(cell_t)*element_len);
}
void Element::assign(Element& a,const Element &b){
	memcpy(&a,&b,sizeof(cell_t)*element_len);
}

void Element::c_vecMul(Element* a, Element* b, Element* c, len_t l){
	for(unsigned int i = 0  ; i < l ; ++i){
		Element::c_mul(a[i],b[i],c[i]);
	}
}

bool Element::c_isZero(Element& a){
	for(unsigned int i = 0 ; i  < element_len ; ++i)
		if(a.c[i]!=0)
			return false;
	return true;
}
void Element::naiveMul(const cell_t* a, const cell_t* b, cell_t* c){
	cell_t clmul[element_len<<1];
	for(unsigned int i = 0 ; i < (element_len<<1) ; ++i)
		clmul[i]=0;
	for(unsigned int i = 0 ; i < ord ; ++i)
		for(unsigned int j = 0 ; j < ord ; ++j){
			if(EXTRACT_BIT(a,i) &&
					EXTRACT_BIT(b,j) )
				clmul[(i+j)>>log_bits_in_cell]^=(((cell_t)1)<<((i+j) & andMask(log_bits_in_cell)));
		}
	for(unsigned int i = (ord<<1)-1 ; i>=ord ; --i){
		if(EXTRACT_BIT(clmul,i)){
			for(unsigned int j = 0 ; j <  mod_len[ord>>5];++j)
				TOGGLE_BIT(clmul,i-ord + irr_poly_index[ord>>5][j]);
		}
	}
	unsigned int i;
	for(i = 0 ; i < element_len ; ++i)
		c[i]=clmul[i];
}

void Element::naiveMul(const Element& a, const Element& b, Element& c){
	naiveMul(a.c,b.c,c.c);
}
} /* namespace FFF */
