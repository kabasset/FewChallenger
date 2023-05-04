/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLINECHALLENGER_SPLINE_H
#define _SPLINECHALLENGER_SPLINE_H

#include "LinxCore/Raster.h"
#include "LinxCore/Tiling.h"

#include <vector>

namespace Spline {

class SplineDomain {
  friend class SplineKnots;
  friend class Spline;

public:
  template <typename T>
  SplineDomain(const T& t) : m_t(t.begin(), t.end()), m_h(m_t.size() - 1), m_g(m_t.size() - 1) {
    // TODO assert size > 2
    for (std::size_t i = 0; i < size; ++i) {
      const double hi = m_t[i + 1] - m_t[i];
      // FIXME assert hi > 0
      m_h[i] = hi;
      m_g[i] = 1. / hi;
    }
  }

  std::size_t index(double x) const {
    if (x < m_t[0]) {
      throw std::exception("x is too small!");
    }
    if (x > m_t[m_t.size() - 1]) {
      throw std::exception("x is too large!");
    }
    std::size_t i = 0;
    while (x >= m_t[i]) {
      ++i;
    }
    return i;
  }

private:
  std::vector<double> m_t; ///< The knot positions
  std::vector<double> m_h; ///< The knot spacings
  std::vector<double> m_g; ///< The inverse of the knot spacings
};

template <typename T>
class SplineKnots {
  friend class Spline;

public:
  SplineKnots(const SplineDomain& t, const TIn& f) : m_domain(t), m_f(f.begin(), f.end()), m_z(m_f.size()) {
    const auto size = m_f.size();
    std::vector<T> d(size); // First derivatives
    d[0] = (m_f[1] - m_f[0]) * m_domain.m_g[0]; // Because next loop starts at 1
    const auto* fIt = m_f.data() + 1;
    const auto* gIt = m_domain.m_g.data() + 1;
    const auto* dIt = d.data() + 1;
    for (auto zIt = &m_z[0] + 1; zIt != &m_z[size - 1]; ++zIt, ++fIt, ++gIt, ++dIt) {
      // d[i] = (m_f[i + 1] - m_f[i]) / m_h[i];
      *dIt = (*(fIt + 1) - *fIt) * *gIt;
      // z[i] = (m_f[i + 1] - m_f[i]) / m_h[i] - (m_f[i] - m_f[i - 1]) / m_h[i - 1];
      *zIt = *dIt - *(dIt - 1);
    }
    // z[0] and z[size-1] are left at 0 for natural splines
  }

private:
  const SplineDomain& m_domain; ///< The knots domain
  std::vector<T> m_f; ///< The knot values
  std::vector<T> m_z; ///< The knot second derivatives
}

template <typename T>
class Spline {
public:
  Spline(const SplineKnots<T>& knots) : m_knots(knots) {}

  T operator()(double x) {
    const auto i = m_knots.m_domain.index(x);
    const double cf0 = right / hi;
    const double cf1 = left / hi;
    const double cz0 = right * right * right / (6. * hi) - hi * right / 6.;
    const double cz1 = left * left * left / (6. * hi) - hi * left / 6.;
    return m_knots.m_f[i] * cf0 + m_knots.m_f[i + 1] * cf1 + m_knots.m_z[i] * cz0 + m_knots.m_z[i + 1] * cz1;
  }

private:
  const SplineKnots<T>& m_knots;
};

} // namespace Spline

#endif
