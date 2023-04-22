/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _CHALLENGER_SPLINE_H
#define _CHALLENGER_SPLINE_H

#include "LinxCore/Raster.h"

#include <complex>

namespace Challenger {

class ComplexSpline2D {
public:
  template <typename TSeq>
  Spline2D(const TSeq& x, const TSeq& y, const TSeq& u, const TSeq& v);

  template <typename TMap>
  Linx::Raster<std::complex<double>, 2> operator()(const TMap& f) const;
};

} // namespace Challenger

#endif