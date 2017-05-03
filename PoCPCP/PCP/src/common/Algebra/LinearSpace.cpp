#include "LinearSpace.hpp"
#include "vanishingPolynomialWrapper.hpp"
#include <omp.h>

namespace Algebra{

#ifdef _MSC_VER // Visual C++
	/* Workaround (shaul): This #ifdef is due to (yet another) violation of the OpenMP standard in the VC++ compiler.
	The threadprivate clause in VC++ is implemented directly as the application of the thread __declspec attribute */
	extern "C" __declspec(thread) unsigned int ntlGlobalThreadNum;  /**< Used to pass a unique thread number to NTL for memory management purposes
																	Used by GF2XRegister in NTL/GF2X1.cpp
																	Defined as extern "C" to avoid name mangling when used with C++.  */
#else
	extern "C" unsigned int ntlGlobalThreadNum;                     /**< Used to pass a unique thread number to NTL for memory management purposes
																	Used by GF2XRegister in NTL/GF2X1.cpp
																	Defined as extern "C" to avoid name mangling when used with C++.  */
#pragma omp threadprivate(ntlGlobalThreadNum)
#endif

namespace {
    elementsSet_t convertFromLegacyBasis(const details::Basis& BasisH){
        elementsSet_t newBasis;

        for(size_t i=0; i< BasisH.getSizeOfBasis() ; i++){
            newBasis.insert(FieldElement(BasisH.getBasisElementByIndex(i)));
        }

        return newBasis;
    }
}

LinearSpace::LinearSpace(details::basisShiftPair BasisH) : BasisH_(BasisH), subspacePoly_(convertFromLegacyBasis(BasisH.basis)){};

bool LinearSpace::exist(const std::unique_ptr<const FieldElementPredicate>& pred)const { 
    details::LinearSubset H(BasisH_.basis);

    const size_t sizeOfH = H.getSizeOfField();
    bool isFound = false;
#ifndef _DEBUG
#pragma omp parallel
#endif
    {
        const size_t numThreads = omp_get_num_threads();
        const size_t currThreadId = omp_get_thread_num();
        ntlGlobalThreadNum = currThreadId;
        for (size_t i = 0; ((i+currThreadId)< sizeOfH) && (!isFound); i+=numThreads){
            const auto currLocation = i+currThreadId;
            FieldElement e(H.getFieldElementByIndex(currLocation) + BasisH_.affineShift);
            if (pred->test(e) == true){
                isFound=true;
            }
            
#ifdef PRINT_PROGRESS
	    if (currLocation%100 == 0){
                std::cout<<currLocation<<"/"<<sizeOfH<<",";
                std::fflush(stdout);
            }
#endif //PRINT_PROGRESS
        }
    }
    return isFound;
}

size_t LinearSpace::size()const {
    details::LinearSubset H(BasisH_.basis);
    return H.getSizeOfField();
}

std::unique_ptr<DivisorPolynomial> LinearSpace::vanishingPoly()const{
    return std::unique_ptr<DivisorPolynomial>(new vanishingPolynomialWrapper(BasisH_.basis, BasisH_.affineShift));
}


bool LinearSpace::contains(const FieldElement& e)const{
    // "e" is in the linear space
    // if and only if it is a root of the
    // minimal polynomial that vanishes over the linear space
    return zero() == subspacePoly_.eval(e);
}

} // namespace Algebra
