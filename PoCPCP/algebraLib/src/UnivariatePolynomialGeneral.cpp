#include "algebraLib/UnivariatePolynomialGeneral.hpp"
#include "algebraLib/FFT.hpp"
#include "algebraLib/ErrorHandling.hpp"

#include <Polynomials.h>

#include <memory>
#include <math.h>

using std::unique_ptr;
using std::vector;
using std::max;

namespace Algebra{
			
UnivariatePolynomialGeneral::UnivariatePolynomialGeneral(const UnivariatePolynomialInterface& src):polynomial_(src.getCoefficients()){ }

unique_ptr<PolynomialInterface> UnivariatePolynomialGeneral::clone()const{
    return unique_ptr<PolynomialInterface>(new UnivariatePolynomialGeneral(*this));
}

UnivariatePolynomialGeneral::UnivariatePolynomialGeneral(const FieldElement& constant):polynomial_(1) {
	setCoefficient(0,constant);
}

UnivariatePolynomialGeneral::UnivariatePolynomialGeneral(const vector<FieldElement>& coeffs):polynomial_(coeffs) {};
UnivariatePolynomialGeneral::UnivariatePolynomialGeneral(vector<FieldElement>&& coeffs):polynomial_(std::move(coeffs)) {};

UnivariatePolynomialGeneral::UnivariatePolynomialGeneral(const elementsSet_t& roots):polynomial_(1) {
	typedef elementsSet_t set;
	setCoefficient(0,one());

	/**
	 * For each root \f$r\f$ generate
	 * the affine polynomial \f$q(x) = x - r\f$
	 * and multiply current polynomial 
	 * by it
	 */
	for(const auto& element : roots){
		
        UnivariatePolynomialGeneral affinePoly(element);
        affinePoly.setCoefficient(1,one());
        multiply(affinePoly);
	}
}

UnivariatePolynomialGeneral::UnivariatePolynomialGeneral(const vector<FieldElement>& evaluation, const vector<FieldElement>& orderedBasis, const FieldElement& spaceShift): polynomial_(IFFT(evaluation,orderedBasis,spaceShift)){};

UnivariatePolynomialGeneral::UnivariatePolynomialGeneral(const evaluation_t& evaluationTabel, const elementsSet_t& spaceBasis){
    
    const short basisSize = spaceBasis.size();
    const size_t spaceSize = (1<<basisSize);
    vector<FieldElement> orderedBasis;
    for(const auto& b : spaceBasis) orderedBasis.push_back(b);

    //construct the evaluation
    vector<FieldElement> vals(spaceSize);
    for(size_t i=0; (i>>basisSize)==0; i++){
        const auto iter = evaluationTabel.find(getSpaceElementByIndex(orderedBasis, zero(), i));
        ALGEBRALIB_ASSERT(iter != evaluationTabel.end(),"Bad evaluation table, imposible interpolation");
        vals[i] = iter->second;
    }

    //compute the IFFT
    polynomial_ = IFFT(vals,orderedBasis,zero());
}

UnivariatePolynomialGeneral::UnivariatePolynomialGeneral(const evaluation_t& evaluationTabel){

    NTL::vec_GF2E x_values, y_values;

    //translate the evaluationTable parameter
    //  details::UniPolynomialEvaluation d_table;
    for(const auto& point : evaluationTabel){
        NTL::append(x_values,NTL::GF2E(point.first));
        NTL::append(y_values,NTL::GF2E(point.second));
    }

    const auto ntlPoly = NTL::interpolate(x_values,y_values);

    //move the ntl coefficient to the model
    polynomial_.resize(1+deg(ntlPoly));
    for(size_t i=0; i< polynomial_.size(); i++){
        polynomial_[i] = FieldElement(coeff(ntlPoly,i));
    }
}
                        

UnivariatePolynomialGeneral::UnivariatePolynomialGeneral(const NTL::GF2EX& p):polynomial_(1+deg(p)){
    for(size_t i=0; i< polynomial_.size(); i++){
        polynomial_[i] = FieldElement(coeff(p,i));
    }
}

UnivariatePolynomialGeneral::operator const NTL::GF2EX()const{
    
    NTL::GF2EX res;

    for(size_t i=0; i< polynomial_.size(); i++){
        SetCoeff(res,i,NTL::GF2E(polynomial_[i]));
    }

    return res;
}

FieldElement UnivariatePolynomialGeneral::eval(const FieldElement& x)const {
	//use horner evaluation
    FieldElement res = zero();
    for(int i=polynomial_.size()-1; i>=0; i--){
        res = res*x + polynomial_[i];
    }
    
    return res;
}

std::vector<FieldElement> UnivariatePolynomialGeneral::eval(const std::vector<FieldElement>& orderedBasis, const FieldElement& shift)const{
    return FFT(polynomial_,orderedBasis,shift);
}

const std::vector<FieldElement> UnivariatePolynomialGeneral::getCoefficients()const{
    return polynomial_;
}

UnivariatePolynomialInterface* UnivariatePolynomialGeneral::eval(const UnivariatePolynomialInterface& p)const{
	UnivariatePolynomialGeneral x(p);
	UnivariatePolynomialGeneral* result = new UnivariatePolynomialGeneral();

	//uses Horner Evaluation
	for (PolynomialDegree::integral_t i= PolynomialDegree::integral_t(getDegree()); i >= 0; i--){
		(*result) = getCoefficient(i) + (x*(*result));
	}

	return result;
}

void UnivariatePolynomialGeneral::setCoefficient(const unsigned int index,const FieldElement& coefficient) {
    if (index >= polynomial_.size()){
        polynomial_.resize(index+1,zero());
    }
    
    polynomial_[index] = coefficient;
}

FieldElement UnivariatePolynomialGeneral::getCoefficient(const unsigned int index)const {
    if(index >= polynomial_.size()){
        return zero();
    }
    return polynomial_[index];
}

PolynomialDegree UnivariatePolynomialGeneral::getDegree()const {
    
    //if there are no coefficients this is the ZERO polynomial
    if(polynomial_.size() == 0){
        return PolynomialDegree::getZeroPolyDegree();
    }
    
    //find the highest non zero coefficient index
    for(int i=polynomial_.size()-1; i>=0; i--){
        if (polynomial_[i] != zero()){
            return PolynomialDegree(i);
        }
    }

    //no non-zero coefficients found
    return PolynomialDegree::getZeroPolyDegree();
}
	
unique_ptr<UnivariatePolynomialInterface> UnivariatePolynomialGeneral::divideByMe(const UnivariatePolynomialInterface& dividend)const {
	//generate an instance of an NTL polynomials
	//not very efficient, but saves programming time,
	//later division can be done directly
	NTL::GF2EX tmp;
	NTL::GF2EX q;
	//ARIEL:doing division only if it's not the zero polynomial (in which case degree is -infinity)
	if (dividend.getDegree().isInteger()){
		for (size_t i = 0; i <= PolynomialDegree::integral_t(dividend.getDegree()); i++){
			SetCoeff(tmp, i, (NTL::GF2E)dividend.getCoefficient(i));
		}

		//calculate the division
		div(q, tmp, NTL::GF2EX(*this));
	}
	else{
		ALGEBRALIB_FATAL("Tried to divide by zero polynomial");
    }
	//return result
	return unique_ptr<UnivariatePolynomialInterface>(new UnivariatePolynomialGeneral((NTL::GF2EX)q));	
}

void UnivariatePolynomialGeneral::multiply(const UnivariatePolynomialGeneral& other){
    //use FFT for multiplication:
    //first evaluate both polynomials on a big enough space,
    //multiplicate the evaluations,
    //and interpolate it back using IFFT
	
    const PolynomialDegree degreeBound = degreeOfProduct(getDegree(),other.getDegree());

    //if expected degree is -infinity, this is the zero polynomial
    if(!degreeBound.isInteger()){
        polynomial_.resize(0);
        return;
    }
    
    //else find a big enough evaluation space
    const size_t space_dim = ceil(log2(1+ PolynomialDegree::integral_t(degreeBound)));
    const auto basis = getStandartBasis(space_dim);
    vector<FieldElement> orderedBasis(basis.begin(),basis.end());

    //evaluate both polys
    vector<FieldElement> eval1 = FFT(polynomial_,orderedBasis,zero());
    const vector<FieldElement> eval2 = FFT(other.polynomial_,orderedBasis,zero());


    //compute the product evaluation
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
    for(size_t i=0; i< eval2.size(); i++){
        eval1[i] *= eval2[i];
    }

    //interpolate to get the coefficients
    polynomial_ = IFFT(eval1,orderedBasis,zero());
}

void UnivariatePolynomialGeneral::add(const UnivariatePolynomialGeneral& other){
    const PolynomialDegree degOfOther = other.getDegree();

    //if the other is the zero poly, nothing should be done
    if(!degOfOther.isInteger())return;

    //otherwise, the other polynomial has effective coefficients.
    //extend corrent poly len if needed
    const size_t otherNumCoeffs = 1+PolynomialDegree::integral_t(degOfOther);
    if(polynomial_.size() < otherNumCoeffs){
        polynomial_.resize(otherNumCoeffs,zero());
    }

    //add the coefficients pointwise
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
    for(long long i=0; i< otherNumCoeffs ; i++){
        polynomial_[i] += other.polynomial_[i];
    }
}

UnivariatePolynomialGeneral operator+(const UnivariatePolynomialGeneral& a, const UnivariatePolynomialGeneral& b) { 
    UnivariatePolynomialGeneral res(a);
    res.add(b);
    return res;
}

UnivariatePolynomialGeneral& operator+=(UnivariatePolynomialGeneral& x, const UnivariatePolynomialGeneral& b) {
    x.add(b);
    return x;
}

UnivariatePolynomialGeneral operator-(const UnivariatePolynomialGeneral& a, const UnivariatePolynomialGeneral& b)
{ return a+b; }

UnivariatePolynomialGeneral operator*(const UnivariatePolynomialGeneral& a, const UnivariatePolynomialGeneral& b) {
    UnivariatePolynomialGeneral res(a);
    res.multiply(b);
    return res;
}

UnivariatePolynomialGeneral& operator*=(UnivariatePolynomialGeneral& x, const UnivariatePolynomialGeneral& b) { 
    x.multiply(b);
    return x;
}

bool operator==(const UnivariatePolynomialGeneral& a, const UnivariatePolynomialGeneral& b) { 
    //if the degrees differ, they can't be the same
    const PolynomialDegree aDeg = a.getDegree();
    const PolynomialDegree bDeg = b.getDegree();
    if(aDeg != bDeg)return false;

    const size_t numCoeffs = max(PolynomialDegree::integral_t(0),1+PolynomialDegree::integral_t(aDeg));
    for(size_t i=0; i< numCoeffs; i++){
        if(a.getCoefficient(i) != b.getCoefficient(i))return false;
    }

    //if difference not found
    return true;
}

bool operator!=(const UnivariatePolynomialGeneral& a, const UnivariatePolynomialGeneral& b)
{ return !(a == b); }

UnivariatePolynomialGeneral power(const UnivariatePolynomialGeneral& a, long e)
{ return UnivariatePolynomialGeneral(power(NTL::GF2EX(a),e)); }

std::ostream& operator<<(std::ostream& s, const UnivariatePolynomialGeneral& a)
{ return s << NTL::GF2EX(a); }

} //namespace Algebra
