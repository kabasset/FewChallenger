/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _MEMORYCHALLENGER_AMPLITUDE_H
#define _MEMORYCHALLENGER_AMPLITUDE_H

#include "Linx/Data/Raster.h"
#include "MemoryChallenger/Interpolant.h"

#include <array>

namespace MemoryChallenger {

static constexpr Linx::Index Ny = 50;
static constexpr Linx::Index Ne = 33;

class Ellipse {

public:
  Ellipse(double p = 0, double e = 0) : m_p {p}, m_e {e}, m_y {std::log(p - 2. * e - 2.1)} {}

  double p() const {
    return m_p;
  }

  double e() const {
    return m_e;
  }

  double y() const {
    return m_y;
  }

  friend bool operator==(const Ellipse& lhs, const Ellipse& rhs) {
    return lhs.m_p == rhs.m_p && lhs.m_e == rhs.m_e;
  }

private:
  double m_p;
  double m_e;
  double m_y;
};

/**
 * l: 2..l_max
 * m: 0..l
 * n: -n_max..n_max
 */
using Mode = Linx::Position<3>;

template <typename TSeq, typename TMap>
ComplexInterpolant2D createInterpolant(const Linx::Position<3>& lmn, const TSeq& ys, const TSeq& es, TMap&& z) {

  (void)lmn; // Silent unused warning (would be used to load from file)

  Linx::Raster<double, 3> zipped({2, Ny, Ne});
  for (int e = 0; e < Ne; ++e) {
    for (int y = 0; y < Ny; ++y) {
      zipped[{0, y, e}] = y + e;
      zipped[{1, y, e}] = 0;
    }
  }

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
  for (Linx::Index i = 0; i < Ny; i++) {
    const double p = 10 * (i + 1);
    const double e = 0;
    ys[i] = std::log(p - 2. * e - 2.1);
  }

  Linx::Raster<double, 1> es({Ne});
  for (Linx::Index i = 0; i < Ne; i++) {
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
  Linx::AlignedRaster<ComplexInterpolant2D, 3> m_interpolants;

  AmplitudeCarrier(Linx::Index lmax, Linx::Index nmax) :
      m_lmax(lmax), m_nmax(nmax), m_z({Ny, Ne, 2, modeCount(m_lmax, m_nmax)}, nullptr, 0),
      m_interpolants({m_nmax * 2 + 1, m_lmax + 1, m_lmax - 1}, nullptr, 0) {
    fillInterpolantRaster(m_lmax, m_nmax, m_z, m_interpolants);
  }

  template <typename TE, typename TM>
  Linx::Raster<std::complex<double>> interpolate(const TE& ellipses, const TM& modes) {
    Linx::Raster<std::complex<double>> amplitude({modes.size(), ellipses.size()});
    auto it = amplitude.begin();
    for (const auto& py : ellipses) {
      for (const auto& lmn : modes) {
        *it = m_interpolants[{lmn[2] + m_nmax, lmn[1], lmn[0] - 2}](py.y(), py.e());
        ++it;
      }
    }
    return amplitude;
  }
};

} // namespace MemoryChallenger

#endif // _FEW_AMPLITUDE_H
