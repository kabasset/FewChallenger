/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLINECHALLENGER_AMPLITUDE_H
#define _SPLINECHALLENGER_AMPLITUDE_H

#include "LinxCore/Raster.h"
#include "MemoryChallenger/Amplitude.h" // Ellipse, Mode
#include "SplineChallenger/Interpolant.h"

#include <array>

namespace SplineChallenger {

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
  Linx::Raster<std::complex<double>> out({Ny, Ne});
  for (int e = 0; e < Ne; ++e) {
    for (int y = 0; y < Ny; ++y) {
      out[{y, e}] = {y + e, 0};
    }
  }
  return out;
}

} // namespace SplineChallenger

#endif // _FEW_AMPLITUDE_H
