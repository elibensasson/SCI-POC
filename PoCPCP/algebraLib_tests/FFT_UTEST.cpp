#include <algebraLib/FFT.hpp>
#include <gtest/gtest.h>

namespace{

using Algebra::FieldElement;
using Algebra::zero;
using Algebra::elementsSet_t;
using Algebra::generateRandom;
using Algebra::getStandartBasis;
using Algebra::IFFT;
using std::vector;

FieldElement getElementFromOrderedBasis(const vector<FieldElement>& orderedBasis, size_t elementIndex){
    FieldElement res = zero();
    for(const auto& b: orderedBasis){
        if (elementIndex%2 == 1) res+=b;
        elementIndex>>=1;
    }

    return res;
}

FieldElement HornerEvaluation(const vector<FieldElement>& poly, const FieldElement x){
    FieldElement res = zero();

    for(long i = poly.size()-1; i>=0 ; i--){
        res = (res*x) + poly[i];
    }

    return res;
}

TEST(Algebra,IFFT){
    const unsigned short basisSize = 10;
    
    //construct the basis (just the standart basis)
    elementsSet_t basisUnordered = getStandartBasis(basisSize);
    vector<FieldElement> orderedBasis;
    for(const auto& b : basisUnordered) orderedBasis.push_back(b);
    const FieldElement shift = generateRandom();

    //construct the evaluation
    vector<FieldElement> vals;
    for(size_t i=0; (i>>basisSize)==0; i++){
        vals.push_back(generateRandom());
    }

    //compute the IFFT
    vector<FieldElement> polyCoeffs = IFFT(vals,orderedBasis,shift);

    //Test result
    for(size_t i=0; (i>>basisSize)==0; i++){
        const FieldElement x = shift + getElementFromOrderedBasis(orderedBasis, i);
        const FieldElement y = HornerEvaluation(polyCoeffs,x);
        EXPECT_EQ(vals[i],y);
    }
}

TEST(Algebra,FFT){
    const unsigned short basisSize = 10;
    const size_t polyDeg = rand() % (1<<basisSize);
    
    //construct the basis (just the standart basis)
    elementsSet_t basisUnordered = getStandartBasis(basisSize);
    vector<FieldElement> orderedBasis;
    for(const auto& b : basisUnordered) orderedBasis.push_back(b);
    const FieldElement shift = generateRandom();

    //construct the evaluation
    vector<FieldElement> polyCoeffs;
    for(size_t i=0; i < polyDeg; i++){
        polyCoeffs.push_back(generateRandom());
    }

    //compute the FFT
    vector<FieldElement> vals = FFT(polyCoeffs,orderedBasis,shift);

    //Test result
    for(size_t i=0; (i>>basisSize)==0; i++){
        const FieldElement x = shift + getElementFromOrderedBasis(orderedBasis, i);
        const FieldElement y = HornerEvaluation(polyCoeffs,x);
        EXPECT_EQ(vals[i],y);
    }
}

TEST(Algebra,FFT_LDE){
    const unsigned short basisSize_src = 10;
    const unsigned short basisSize_dst = basisSize_src + (rand() % 3);
     
    //construct the basis (just the standart basis)
    vector<FieldElement> orderedBasis_src;
    {
    elementsSet_t basisUnordered = getStandartBasis(basisSize_src);
    for(const auto& b : basisUnordered) orderedBasis_src.push_back(b);
    }
    const FieldElement shift_src = generateRandom();
    
    vector<FieldElement> orderedBasis_dst;
    {
    elementsSet_t basisUnordered = getStandartBasis(basisSize_dst);
    for(const auto& b : basisUnordered) orderedBasis_dst.push_back(b);
    }
    const FieldElement shift_dst = generateRandom();

    //construct the evaluation
    vector<FieldElement> vals;
    for(size_t i=0; (i>>basisSize_src)==0; i++){
        vals.push_back(generateRandom());
    }

    //comput the LDE
    vector<FieldElement> LDE_res = LDE(vals, orderedBasis_src, shift_src, orderedBasis_dst, shift_dst);

    //compute the reference
    vector<FieldElement> polyCoeffs = IFFT(vals,orderedBasis_src,shift_src);
    vector<FieldElement> ref = FFT(polyCoeffs,orderedBasis_dst,shift_dst);

    EXPECT_EQ(ref.size(),LDE_res.size());
    for(size_t i=0; i < ref.size(); i++){
        EXPECT_EQ(ref[i], LDE_res[i]);
    }

}

}
