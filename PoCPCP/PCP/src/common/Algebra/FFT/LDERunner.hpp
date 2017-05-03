#pragma once

#ifndef __LDERUNNER_HPP__
#define __LDERUNNER_HPP__

#include "common/Algebra/LightUniPolyEval.hpp"

namespace Algebra{
/*generates a low degree extension using Matan's code - first converts to Matan's format,
* applies FFT and inverse FFT, and converts back to NTL. 
*/

//same as above - but uses LightUniPolyEval rather than UniPolynomialEvaluation. This one also uses navie LDE for small instances
void LDERunner(FieldElement const * const ntl_interpolate_eval,
	const std::vector<Algebra::FieldElement>					&ntl_interpolate_basis,
	const Algebra::FieldElement				&ntl_interpolate_shift,
	Algebra::FieldElement* ntl_eval_eval,
	const std::vector<Algebra::FieldElement>                   &ntl_eval_basis,
	const Algebra::FieldElement              &ntl_eval_shift

	);

//just IFFT. Gets light evaluation and stores result in ntl_poly
void IFFTRunner(const LightUniPolyEval &ntl_eval,
	const Algebra::details::Basis					&ntl_basis,
	const Algebra::FieldElement				&ntl_shift,
	details::UnivariatePolynomial& ntl_poly
	);
} //namespace Algebra
#endif	// #ifndef __LDERUNNER_HPP__
