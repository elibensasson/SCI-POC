/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef ALGEBRALIB_TRACE_HPP_
#define ALGEBRALIB_TRACE_HPP_

#include <algebraLib/FieldElement.hpp>

namespace Algebra {

extern const FieldElement invExtrConsts[];

FieldElement extractBit(const FieldElement& elem, const int bitNum);

} // namespace Algebra

#endif // ALGEBRALIB_TRACE_HPP_
