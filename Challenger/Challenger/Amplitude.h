/// @copyright 2020, Niels Warburton, Michael L. Katz, Alvin J.K. Chua, Scott A. Hughes
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _CHALLENGER_AMPLITUDE_H
#define _CHALLENGER_AMPLITUDE_H

#include "Challenger/Interpolant.h"
#include "LinxCore/Raster.h"

#include <array>

namespace Challenger {

struct PE {
  double p;
  double e;

  inline double y() const {
    return std::log(p - 2. * e - 2.1);
  }
};

struct LMN {
  Linx::Index l; // 2..l_max
  Linx::Index m; // 0..l
  Linx::Index n; // -n_max..n_max
};

Linx::Raster<ComplexInterpolant2D, 3> load_and_interpolate_amplitude_data(lmax, nmax, m_interpolants); // FIXME

class AmplitudeCarrier {

public:
  Linx::Index m_lMax;
  Linx::Index m_nMax;
  Linx::Raster<ComplexInterpolant2D, 3> m_interpolants; // FIXME Single ComplexSpline2D?

  AmplitudeCarrier(Linx::Index lmax, Linx::Index nmax) :
      m_interpolants(load_and_interpolate_amplitude_data(lmax, nmax, m_interpolants)), m_lMax(lmax), m_nMax {nmax} {}

  template <typename Tpe, typename Tlmn>
  Linx::Raster<std::complex<double>> interpolate(const Tpe& pes, const Tlmn& lmns) {
    // FIXME is lmns a mere list of indices?
    Linx::Raster<std::complex<double>> amplitude({lmns.size(), pes.size()});
    auto it = amplitude.begin();
    for (const auto& pe : pes) { // FIXME structured binding
      const auto p = pe.p;
      const auto e = pe.e;
      const auto y = pe.y();
      Linx::Position<3> p;
      for (Linx::Index l = 2; l <= m_lMax; ++l) {
        for (Linx::Index m = 0; m <= l; ++m) {
          for (Linx::Index n = -m_nMax; n <= m_nMax; ++n) {
            p = {l - 2, m, n + m_nMax}; // FIXME compactify
            amplitude[p] = m_interpolants[p](y, e);
          }
        }
      }
    }
  }
};

} // namespace Challenger

#endif // _FEW_AMPLITUDE_H
