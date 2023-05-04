/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLINECHALLENGER_SPLINE_H
#define _SPLINECHALLENGER_SPLINE_H

#include "LinxCore/Raster.h"
#include "LinxCore/Tiling.h"

#include <stdexcept>
#include <vector>

namespace Focs {

class SplineDomain {
  template <typename>
  friend class Spline;

public:
  template <typename TIt>
  SplineDomain(const TIt& begin, const TIt& end) : m_u(begin, end), m_h(m_u.size() - 1), m_g(m_u.size() - 1) {
    // TODO assert size > 2
    for (std::size_t i = 0; i < m_h.size(); ++i) {
      const double h = m_u[i + 1] - m_u[i];
      // FIXME assert h > 0
      m_h[i] = h;
      m_g[i] = 1. / h;
    }
  }

  SplineDomain(std::initializer_list<double> u) : SplineDomain(u.begin(), u.end()) {}

  std::size_t index(double x) const {
    if (x < m_u[0]) {
      throw std::runtime_error("x is too small!");
    }
    if (x > m_u[m_u.size() - 1]) {
      throw std::runtime_error("x is too large!");
    }
    std::size_t i = 0;
    while (x >= m_u[i]) {
      ++i;
    }
    return i;
  }

private:
  std::vector<double> m_u; ///< The knot positions
  std::vector<double> m_h; ///< The knot spacings
  std::vector<double> m_g; ///< The inverse of the knot spacings
};

template <typename T>
class Spline {

public:
  template <typename TValues>
  Spline(const SplineDomain& domain, const TValues& v) : m_domain(domain), m_v(v.begin(), v.end()), m_z(m_v.size()) {
    const auto size = m_v.size();
    // FIXME check size == doman.m_u.size()
    std::vector<T> d(size); // First derivatives
    d[0] = (m_v[1] - m_v[0]) * m_domain.m_g[0]; // Because next loop starts at 1
    auto* dIt = d.data() + 1;
    const auto* vIt = m_v.data() + 1;
    const auto* gIt = m_domain.m_g.data() + 1;
    for (auto zIt = &m_z[0] + 1; zIt != &m_z[size - 1]; ++zIt, ++vIt, ++gIt, ++dIt) {
      // d[i] = (m_v[i + 1] - m_v[i]) / m_h[i];
      *dIt = (*(vIt + 1) - *vIt) * *gIt;
      // z[i] = (m_v[i + 1] - m_v[i]) / m_h[i] - (m_v[i] - m_v[i - 1]) / m_h[i - 1];
      *zIt = *dIt - *(dIt - 1);
    }
    // z[0] and z[size-1] are left at 0 for natural splines
  }

  T operator()(double x) const {
    // FIXME separate dependency on x and on v
    const auto i = m_domain.index(x);
    const double left = x - m_domain.m_u[i];
    const double right = m_domain.m_u[i + 1] - x;
    const double h = left + right;
    const double cv0 = right / h;
    const double cv1 = left / h;
    const double cz0 = right * right * right / (6. * h) - h * right / 6.;
    const double cz1 = left * left * left / (6. * h) - h * left / 6.;
    return m_v[i] * cv0 + m_v[i + 1] * cv1 + m_z[i] * cz0 + m_z[i + 1] * cz1;
  }

private:
  const SplineDomain& m_domain; ///< The knots domain
  std::vector<T> m_v; ///< The knot values
  std::vector<T> m_z; ///< The knot second derivatives
};

} // namespace Focs

#endif
