#include "AlgebraCommon.hpp"

#include <vector>
#include <set>


/* The purpose of this class is to concentrate various convenient methods related to the algebraic classes.
* In particular, many useful conversions between the older algebraic objects used mainly
in the PCP code, to newer algebraic objects used a lot in the ACSP and BREX->ACSP code (and perhaps elsewhere)
*/




using Algebra::details::Basis;
using Algebra::FieldElement;

namespace Algebra{

	/**********************CONVERSIONS ****************************************/

	vector<FieldElement> oldBasisToElementVector(const Algebra::details::Basis& b){
		vector<FieldElement> v;
		for (size_t i = 0; i < b.getSizeOfBasis(); i++)
			v.push_back(FieldElement(b.getBasisElementByIndex(i)));
		return v;
	}

	/*vector<FieldElement> FElemVectorToFieldElementVector(vector<FElem> v){
		vector<FieldElement> vec;
		for (auto it = v.begin(); it != v.end(); it++){
			vec.push_back(FieldElement(NTL::GF2E(*it)));
		}
		return vec;
	}*/
	//converts FieldElement Array to FieldElement Vector - TODO:check if can do this faster using Vector constructor - but then need to be careful not erasing c
	vector<FieldElement> fieldPointArrayToPointVector(const FieldElement* c, size_t size){
		vector<FieldElement> v;
		for (size_t i = 0; i < size; i++)
			v.push_back(c[i]);
		return v;
	}


	


	/** returns the vector/coeffs of polynomial over GF2 corresponding to the field element e
	subtlety : if an element e of GF_{2^m} has only the first l coeffs non-zero NTL will store
	e as a vector of length l. On the other hand, this method always returns a vector of length m. e.g., in the case of
	such e, the vector returned will end with m-l zeros. */
	NTL::vec_GF2 fieldPointToGF2Vector(const FieldElement& e){

		using NTL::to_vec_GF2;
		NTL::vec_GF2 a = to_vec_GF2(::NTL::GF2E(e)._GF2E__rep);
		NTL::vec_GF2 b(NTL::INIT_SIZE, ::NTL::GF2E(e).degree());
		for (long i = 0; i < a.length(); i++)
			b[i] = a[i];

		return b;
	}

	vector<FieldElement> fieldPointArrayToElementVector(const FieldElement* c, size_t size){
		vector<FieldElement> v;
		for (size_t i = 0; i < size; i++)
			v.push_back(FieldElement(c[i]));
		return v;
	}


	//generates FieldELement x
	FieldElement xElement(){
		NTL::GF2X a(0, 0);
		NTL::SetCoeff(a, 1);
		return FieldElement(NTL::to_GF2E(a));
	}

	//naive powering-don't use for large powers
	FieldElement pow(const FieldElement& a, int b){
		FieldElement res = one();
		for (auto i = 0; i < b; i++)
			res *= a;
		return res;
	}



	/** returns the vector/coeffs of polynomial over GF2 corresponding to the field element e
	subtlety : if an element e of GF_{2^m} has only the first l coeffs non-zero NTL will store
	e as a vector of length l. On the other hand, this method always returns a vector of length m. e.g., in the case of
	such e, the vector returned will end with m-l zeros. */
	NTL::vec_GF2 fieldPointToFullVector(const FieldElement& e){
		NTL::GF2X f = ::NTL::GF2E(e)._GF2E__rep;
		NTL::vec_GF2 a = to_vec_GF2(f);
		NTL::vec_GF2 b(NTL::INIT_SIZE, getNTLFieldDegree());
		for (long i = 0; i < a.rep.length(); i++)
			b.rep[i] = a.rep[i];

		return b;
	}

	//this method not working right now!
	FieldIndex fieldElementToInteger(const FieldElement& e) {
		FieldIndex a = 0;
		NTL::GF2X& x = (NTL::GF2E(e)).LoopHole();//getting NTL poly represenation of e
		long replen = x.xrep.length();
		//_COMMON_ASSERT(replen <= 1, "fieldElementToInteger method currently assumes only 32 non-zero bits in e")
		//if (replen > 0)
			a = x.xrep[1];
		return a;
	}

}//namespace Algebra
