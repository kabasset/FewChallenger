/// @copyright 2020, Niels Warburton, Michael L. Katz, Alvin J.K. Chua, Scott A. Hughes
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _CHALLENGER_AMPLITUDE_H
#define _CHALLENGER_AMPLITUDE_H

#include "Challenger/Interpolant.h"
#include "LinxCore/Raster.h"

#include <array>

namespace Challenger {

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

/**
 * l: 2..l_max
 * m: 0..l
 * n: -n_max..n_max
 */
using Lmn = Linx::Position<3>;

template <typename TSeq, typename TMap>
ComplexInterpolant2D createInterpolant(const Linx::Position<3>& lmn, const TSeq& ys, const TSeq& es, TMap&& z) {

  (void)lmn; // Silent unused warning (would be used to load from file)

  Linx::Raster<double, 3> zipped({2, Ny, Ne});
  for (Linx::Index i = 0; i < Ny * Ne; ++i) { // Not range() for faire comparison with FEW
    zipped[2 * i] = i;
    zipped[2 * i + 1] = -i;
  }

  // FIXME can we avoid this, by refactoring the interpolants?
  for (int e = 0; e < Ne; ++e) {
    for (int y = 0; y < Ny; ++y) {
      std::forward<TMap>(z)[{y, e, 0}] = zipped[{0, Ny - 1 - y, e}];
      std::forward<TMap>(z)[{y, e, 1}] = zipped[{1, Ny - 1 - y, e}];
    }
  }

  return ComplexInterpolant2D(ys, es, z);
}

template <typename TMap, typename TInterpolants>
void fillInterpolantRaster(Linx::Index lmax, Linx::Index nmax, TMap& z, TInterpolants& interpolants) {

  Linx::Raster<double, 1> ys({Ny});
  for (Linx::Index i = 0; i < Ny; i++) { // Not range() for faire comparison with FEW
    ys[i] = i;
  }

  Linx::Raster<double, 1> es({Ne});
  for (Linx::Index i = 0; i < Ne; i++) { // Not range() for faire comparison with FEW
    es[i] = i;
  }

  Linx::Position<3> pos;
  Linx::Index i = 0;
  for (Linx::Index l = 2; l <= lmax; ++l) { // Idem
    for (Linx::Index m = 0; m <= l; ++m) {
      for (Linx::Index n = -nmax; n <= nmax; ++n) {
        pos = {n + nmax, m, l - 2};
        auto section = std::forward<TMap>(z).section(i);
        interpolants[pos] = createInterpolant({l, m, n}, ys, es, section);
        ++i;
      }
    }
  }
}

constexpr Linx::Index modeCount(Linx::Index lmax, Linx::Index nmax) {
  Linx::Index out = 0;
  for (Linx::Index l = 2; l <= lmax; ++l) { // Idem
    for (Linx::Index m = 0; m <= l; ++m) {
      for (Linx::Index n = -nmax; n <= nmax; ++n) {
        ++out;
      }
    }
  }
  return out;
}

class AmplitudeCarrier {

public:
  Linx::Index m_lmax;
  Linx::Index m_nmax;
  Linx::AlignedRaster<double, 4> m_z;
  Linx::AlignedRaster<ComplexInterpolant2D, 3> m_interpolants; // FIXME Single ComplexSpline2D?

  AmplitudeCarrier(Linx::Index lmax, Linx::Index nmax) :
      m_lmax(lmax), m_nmax(nmax), m_z({Ny, Ne, 2, modeCount(m_lmax, m_nmax)}, nullptr, 0),
      m_interpolants({m_nmax * 2 + 1, m_lmax + 1, m_lmax - 1}, nullptr, 0) {
    fillInterpolantRaster(m_lmax, m_nmax, m_z, m_interpolants);
  }

  template <typename Tpe, typename Tlmn>
  Linx::Raster<std::complex<double>> interpolate(const Tpe& pes, const Tlmn& lmns) {
    Linx::Raster<std::complex<double>> amplitude({lmns.size(), pes.size()});
    auto it = amplitude.begin();
    for (const auto& pe : pes) {
      for (const auto& lmn : lmns) {
        *it = m_interpolants[{lmn[2] + m_nmax, lmn[1], lmn[0] - 2}](pe.y, pe.e);
        ++it;
      }
    }
    return amplitude;
  }
};

} // namespace Challenger

#endif // _FEW_AMPLITUDE_H
