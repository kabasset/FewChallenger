/// @copyright 2020, Niels Warburton, Michael L. Katz, Alvin J.K. Chua, Scott A. Hughes
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _CHALLENGER_AMPLITUDE_H
#define _CHALLENGER_AMPLITUDE_H

#include "Challenger/Interpolant.h"
#include "LinxCore/Raster.h"

#include <array>

namespace Challenger {

// The 11 below means the lmax = 10
struct Amplitudes {
  static constexpr Linx::Index lmax = 10;
  std::array<Linx::Raster<Interpolant, 3>, lmax + 1> re;
  std::array<Linx::Raster<Interpolant, 3>, lmax + 1> im;
  // FIXME 4D raster
};

class AmplitudeCarrier {

public:
  Amplitudes m_amplitudes; // FIXME Linx::Raster<Interpolant, 4> m_amplitudes;
  Linx::Index m_nmax;

  AmplitudeCarrier(Linx::Index lmax, Linx::Index nmax) :
      m_amplitudes {load_and_interpolate_amplitude_data(lmax, nmax, m_amplitudes)}, m_nmax {nmax} {}

  template <typename Tpe, typename Tlmn>
  Linx::Raster<std::complex<double>> interpolate(const Tpe& pes, const Tlmn& lmns) {
    // FIXME is lmns a mere list of indices?
    Linx::Raster<std::complex<double>> amplitude({lmns.size(), pes.size()});
    auto it = amplitude.begin();
    for (const auto& pe : pes) { // FIXME structured binding
      const auto p = pe.p;
      const auto e = pe.e;
      const auto y = std::log(p - 2. * e - 2.1);
      for (const auto& lmn : lmns) { // FIXME structured binding
        const auto l = lmn.l;
        const auto m = lmn.m;
        const auto n = lmn.n;
        *it = m_amplitudes.re[{l, m, n + m_nmax}](y, e), m_amplitudes.im[{l, m, n + m_nmax}](y, e);
        // FIXME *it = m_amplitudes[{l, m, n +m_nmax}](y, e);
        ++it;
      }
    }
  }
};

} // namespace Challenger

#endif // _FEW_AMPLITUDE_H
