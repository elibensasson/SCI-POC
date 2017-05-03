/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "algebraLib/BitExtract.hpp"
#include <NTL/mat_GF2.h>

#define DIM 64
//define DIM  (NTL::GF2E::degree()) 
//#define xFE (Algebra::FieldElement(NTL::to_GF2E(NTL::GF2X(1, 1))))
#define zFE (mapIntegerToFieldElement(0, DIM, 0))

namespace Algebra {

	const FieldElement invExtrConsts[ /* 16 */ ] = {
		mapIntegerToFieldElement(0, DIM, 8+2),
		mapIntegerToFieldElement(0, DIM, 8+4),
		mapIntegerToFieldElement(0, DIM, 16+8),
		mapIntegerToFieldElement(0, DIM, 0x5CB972E5CB972E52),
		zFE, zFE, zFE, zFE, zFE,
		zFE, zFE, zFE, zFE, zFE,
		zFE, mapIntegerToFieldElement(0, DIM, 0x6DB6DB6DB6DB0005)
	};
	//const FieldElement invExtrConsts[] = { xFE*xFE*xFE + xFE, xFE*xFE*xFE + xFE*xFE,
	//	xFE*xFE*xFE*xFE + xFE*xFE*xFE, Algebra::mapIntegerToFieldElement(0, 64, 6681497853469601362), Algebra::zero() };

NTL::mat_GF2 matForBitExtr(const int bitNum) {
	NTL::mat_GF2 M;
	NTL::vec_GF2 v;
	const NTL::GF2E invExtractBit = (NTL::GF2E)(invExtrConsts[bitNum]);
	const NTL::GF2E x = NTL::to_GF2E(NTL::GF2X(1, 1));
	M.SetDims(DIM - 1, DIM - 1);
	NTL::GF2E y = x; //not one;
	for (int i = 0; i < DIM - 1; ++i) {
		NTL::GF2E tmp = invExtractBit*(y*y + y);
		v = to_vec_GF2(tmp.LoopHole());
		v.SetLength(DIM);
		for (int j = 0; j < DIM - 1; ++j)
			if (j < bitNum)
				M[i][j] = v[j];
			else
				M[i][j] = v[j + 1];
		//M[i] = v;
		y *= x;
	}
	return M;
}

//TODO: half-trace method ii.2.4
FieldElement extractBit(const FieldElement& elem, const int bitNum) {
	static const NTL::GF2E x = NTL::to_GF2E(NTL::GF2X(1, 1));
	const NTL::mat_GF2 invM = inv(matForBitExtr(bitNum));
	NTL::vec_GF2 tvec = to_vec_GF2(((NTL::GF2E)elem).LoopHole());
	tvec.SetLength(DIM);
	NTL::vec_GF2 vcoeffs;
	vcoeffs.SetLength(DIM - 1);
	for (int j = 0; j < DIM - 1; ++j)
		if (j < bitNum)
			vcoeffs[j] = tvec[j];
		else
			vcoeffs[j] = tvec[j+1];
	NTL::vec_GF2 v = vcoeffs * invM;
	//std::cout << v << std::endl;
	NTL::GF2E y = NTL::GF2E::zero();
	NTL::GF2E tmp = x; //1coeff dontcare
	for (int j = 0; j < DIM - 1; ++j) {
		if (NTL::GF2::zero() != v[j])
			y += tmp;
		tmp *= x;
	}
	return FieldElement(y);
}

} // namespace Algebra
