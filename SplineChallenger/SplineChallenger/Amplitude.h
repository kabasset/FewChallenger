/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLINECHALLENGER_AMPLITUDE_H
#define _SPLINECHALLENGER_AMPLITUDE_H

#include "LinxCore/Raster.h"
#include "SplineChallenger/Interpolant.h"

#include <array>

namespace SplineChallenger {

static constexpr Linx::Index Ny = 50;
static constexpr Linx::Index Ne = 33;

struct Pe {
  const double p;
  const double e;
  const double y;

  explicit Pe(double p_ = 0, double e_ = 0) : p {p_}, e {e_}, y {std::log(p - 2. * e - 2.1)} {}
};

bool operator==(const Pe& lhs, const Pe& rhs) {
  return lhs.p == rhs.p && lhs.e == rhs.e;
}

class AmplitudeCarrier {

public:
  Linx::Index m_lmax;
  Linx::Index m_nmax;
  Linx::AlignedRaster<std::complex<double>, 2> m_z;
  // Interpolant m_interpolant;

  AmplitudeCarrier(Linx::Index lmax, Linx::Index nmax) : m_lmax(lmax), m_nmax(nmax), m_z({Ny, Ne}, nullptr, 0), {}
};

} // namespace SplineChallenger

#endif // _FEW_AMPLITUDE_H
