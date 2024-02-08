/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLINECHALLENGER_AMPLITUDE_H
#define _SPLINECHALLENGER_AMPLITUDE_H

#include "Linx/Data/Raster.h"
#include "MemoryChallenger/Amplitude.h" // Ellipse, Mode
#include "Splider/Lagrange.h"

#include <array>

namespace SplineChallenger {

using Spline = Splider::Lagrange; // Simple local bicubic interpolation

static constexpr Linx::Index Ny = 50;
static constexpr Linx::Index Ne = 33;

Linx::Vector<std::vector<double>, 2> loadGrid() {

  std::vector<double> ys(Ny);
  for (Linx::Index i = 0; i < Ny; i++) {
    const double p = 10 * (i + 1);
    const double e = 0;
    ys[i] = std::log(p - 2. * e - 2.1);
  }

  std::vector<double> es(Ne);
  for (Linx::Index i = 0; i < Ne; i++) {
    es[i] = i;
  }

  return {ys, es};
}

Linx::Raster<std::complex<double>> loadModeData(const MemoryChallenger::Mode& mode) {
  Linx::Raster<std::complex<double>> raw({Ny, Ne});
  for (int e = 0; e < Ne; ++e) {
    for (int y = 0; y < Ny; ++y) {
      raw[{y, e}] = {y + e, 0};
    }
  }
  // FIXME return raw.flip<0>()
  Linx::Raster<std::complex<double>> out(raw.shape());
  for (int e = 0; e < Ne; ++e) {
    for (int y = 0; y < Ny; ++y) {
      out[{y, e}] = raw[{Ny - 1 - y, e}];
    }
  }
  return out;
}

} // namespace SplineChallenger

#endif // _FEW_AMPLITUDE_H
