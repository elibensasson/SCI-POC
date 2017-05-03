#include "LDERunner.hpp"

#include "FFTCommon.hpp"
#include <omp.h>
#if defined _WIN64 || defined _WIN32
#include <windows.h>
#else	// #if defined _WIN64 || defined _WIN32
#endif	// #if defined _WIN64 || defined _WIN32

#include <algebraLib/FFT.hpp>
#include <NTL/GF2E.h>

namespace Algebra{

using NTL::vec_GF2E;
using NTL::GF2EX;
using Infrastructure::POW2;

//naive LDE for small dimensions, where not worth converting format and running FFT
void naiveLDE(const LightUniPolyEval &interpolate_eval,
	const Algebra::details::Basis					&interpolate_basis,
	const Algebra::FieldElement				&interpolate_shift,
	Algebra::FieldElement* eval_eval,
	const Algebra::details::Basis                   &eval_basis,
	const Algebra::FieldElement              &eval_shift,
	const FieldIndex start_index, //for Light eval need to know the index of first element of subspace span(ntl_interpolate_basis) +ntl_interpolate_shift in evaluation
	const FieldIndex increment

	) {

	GF2EX poly;
	/** Naive interpolation used for small inputs. */

	vec_GF2E x_values, y_values;
	FieldIndex i;
	FieldIndex size = POW2(interpolate_basis.getSizeOfBasis());
	FieldElement currentPoint;
	x_values.SetLength(::Infrastructure::safeConvert(size));
	y_values.SetLength(::Infrastructure::safeConvert(size));
	FieldIndex currentIndex = start_index;
	for (i = 0; i < size; i++) {
		currentPoint = interpolate_eval.getPoint(currentIndex);
		x_values[::Infrastructure::safeConvert(i)] = ::NTL::GF2E(currentPoint);
		y_values[::Infrastructure::safeConvert(i)] = ::NTL::GF2E(interpolate_eval.queryAtPoint(currentIndex));
		//cout << "Point " << i << ": " << x_values[::Infrastructure::safeConvert(i)] << endl;
		//cout << "Value " << i << ": " << y_values[::Infrastructure::safeConvert(i)] << endl;
		currentIndex += increment;
	}
	NTL::interpolate(poly, x_values, y_values);
	/** Naive evaluation*/
	size = POW2(eval_basis.getSizeOfBasis());
	for (i = 0; i < size; i++) {
		currentPoint = getSpaceElementByIndex(eval_basis.asVector(),eval_shift,i);
		eval_eval[i] = FieldElement(eval(poly, ::NTL::GF2E(currentPoint)));
	}
}

//naive interpolate for small dimensions, where not worth converting format and running FFT
void naiveInterpolate(const LightUniPolyEval &eval,
	const Algebra::details::Basis					&basis,
	const Algebra::FieldElement				&shift,
	details::UnivariatePolynomial &P,
	const FieldIndex start_index, //for Light eval need to know the index of first element of subspace span(ntl_interpolate_basis) +ntl_interpolate_shift in evaluation
	const FieldIndex increment

	) {

	GF2EX poly;
	/** Naive interpolation used for small inputs. */

	vec_GF2E x_values, y_values;
	FieldIndex i;
	FieldIndex size = POW2(basis.getSizeOfBasis());
	FieldElement currentPoint;
	x_values.SetLength(::Infrastructure::safeConvert(size));
	y_values.SetLength(::Infrastructure::safeConvert(size));
	FieldIndex currentIndex = start_index;
	for (i = 0; i < size; i++) {
		currentPoint = eval.getPoint(currentIndex);
		x_values[::Infrastructure::safeConvert(i)] = ::NTL::GF2E(currentPoint);
		y_values[::Infrastructure::safeConvert(i)] = ::NTL::GF2E(eval.queryAtPoint(currentIndex));
	//	cout << "Point " << i << ": " << x_values[::Infrastructure::safeConvert(i)] << endl;
	//	cout << "Value " << i << ": " << y_values[::Infrastructure::safeConvert(i)] << endl;
		currentIndex += increment;
	}
	NTL::interpolate(poly, x_values, y_values);
	P.setPoly(poly);
	
}

using Algebra::LightUniPolyEval;


//same as above - but uses LightUniPolyEval rather than UniPolynomialEvaluation 
void LDERunner(FieldElement const * const  ntl_interpolate_eval,
	const vector<FieldElement>					&ntl_interpolate_basis,
	const Algebra::FieldElement				&ntl_interpolate_shift,
	Algebra::FieldElement* ntl_eval_eval,
	const vector<FieldElement>                   &ntl_eval_basis,
	const Algebra::FieldElement              &ntl_eval_shift

	) {

    const size_t num_of_elements_src = POW2(ntl_interpolate_basis.size());
    const size_t num_of_bytes_src = num_of_elements_src * sizeof(FieldElement);
    memcpy(ntl_eval_eval,ntl_interpolate_eval,num_of_bytes_src);
    LDE_inplace(ntl_eval_eval,ntl_interpolate_basis,ntl_interpolate_shift,ntl_eval_basis,ntl_eval_shift);
}

//just IFFT. Gets light evaluation and stores result as poly.
void IFFTRunner(
        const LightUniPolyEval &ntl_eval,
        const Algebra::details::Basis& ntl_basis,
        const Algebra::FieldElement& ntl_shift,
        details::UnivariatePolynomial& ntl_poly) 

{
	if (ntl_basis.getSizeOfBasis() < MATAN_IFFT_LOWER_BOUND){
		naiveInterpolate(ntl_eval, ntl_basis, ntl_shift, ntl_poly, 0, 1);
		return;
	}

    std::vector<FieldElement> basis(ntl_basis.getSizeOfBasis());
    for (int i=0; i<ntl_basis.getSizeOfBasis(); i++){
        basis[i] = ntl_basis.getBasisElementByIndex(i);
    }

    std::vector<FieldElement> coeffs = Algebra::IFFT(ntl_eval.getTable(),basis,ntl_shift); 
	
    for (long long i = 0; i < coeffs.size(); ++i) {
		ntl_poly.setCoeff(i,coeffs[i]);
	}
}

} //namespace Algebra
