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

class AmplitudeCarrier {

public:
  Linx::AlignedRaster<std::complex<double>, 2> m_z;
  Interpolant m_interpolant;

  tempate<typename TGrid, typename TTraj>
  AmplitudeCarrier(const TGrid& grid, const TTraj& trajectory) : m_interpolant(grid, trajectory) {}
};

} // namespace SplineChallenger

#endif // _FEW_AMPLITUDE_H
