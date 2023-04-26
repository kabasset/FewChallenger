/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _CHALLENGER_SPLINE_H
#define _CHALLENGER_SPLINE_H

#include "LinxCore/Raster.h"

#include <complex>

namespace Challenger {

/**
 * @brief Spline interpolation and integration.
 * @details
 * This class takes at construction the position of the input knots (u) and output samples (x).
 * Most spline coefficients are precomputed.
 * The call operator then takes as input the value of the knots (v) and return that of the samples (y).
 * Second derivatives of the knots (z) are available byproducts.
 */
class OrientedSpline {

public:
  /**
   * @brief Constructor.
   * @param u The input knots positions (strictly increasing)
   * @param x The output samples positions (strictly increasing and in u's interval)
   * @details
   * Precompute spline coefficients.
   */
  OrientedSpline(std::vector<double> u, std::vector<double> x) :
      m_size(u.size() - 1), m_u(std::move(u)), m_h(m_size), m_g(m_size), m_x(std::move(x)), m_n(m_size),
      m_k(m_x.size()) {
    // FIXME assert m_size > 3
    std::size_t j = 0;
    for (std::size_t i = 0; i < m_size; ++i) {
      const double hi = m_u[i + 1] - m_u[i];
      // FIXME assert hi > 0
      m_h[i] = hi;
      m_g[i] = 1. / hi;
      while (j < m_k.size() && m_x[j] < m_u[i + 1]) {
        const double left = m_x[j] - m_u[i];
        const double right = m_u[i + 1] - m_x[j];
        m_k[j].z1 = left * left * left / (6. * hi) - hi * left / 6.;
        m_k[j].z0 = right * right * right / (6. * hi) - hi * right / 6.;
        m_k[j].v1 = left / hi;
        m_k[j].v0 = right / hi;
        ++m_n[i];
        ++j;
      }
    } // Could be optimized, e.g. with iterators, but executed only once and more readable this way

    m_h0 = m_h[0];
    m_h1 = m_h[1];
    m_h02 = m_h0 * m_h0;
    m_h12 = m_h1 * m_h1;
    m_h01 = m_h0 * m_h1;
    m_hm0 = m_h[m_size - 1];
    m_hm1 = m_h[m_size - 2];
    m_hm02 = m_hm0 * m_hm0;
    m_hm12 = m_hm1 * m_hm1;
    m_hm01 = m_hm0 * m_hm1;
    m_g0 = m_g[0];
    m_zm0Factor = 1. / (2. * m_hm0 * m_hm0 + 3 * m_hm0 * m_hm1 + m_hm1 * m_hm1);
  }

  /**
   * @brief Get the knot x's.
   */
  const std::vector<double>& knotPositions() const {
    return m_u;
  }

  /**
   * @brief Get the interpolation x's.
   */
  const std::vector<double>& samplePositions() const {
    return m_x;
  }

  /**
   * @brief Compute the second derivatives at knots.
   */
  template <typename T>
  std::vector<T> knotSecondDerivatives(const T* v) const {
    std::vector<T> s(m_size);
    std::vector<T> z(m_size + 1);
    s[0] = (v[1] - v[0]) * m_g0; // Because next loop starts at 1
    auto* vIt = v + 1;
    auto* gIt = m_g.data() + 1;
    auto* sIt = s.data() + 1;
    auto* zIt = z.data() + 1;
    for (auto i = m_size - 1; i--; ++vIt, ++gIt, ++sIt, ++zIt) {
      // s[i] = (v[i + 1] - v[i]) / m_h[i];
      *sIt = (*(vIt + 1) - *vIt) * *gIt;
      // z[i] = (v[i + 1] - v[i]) / m_h[i] - (v[i] - v[i - 1]) / m_h[i - 1];
      *zIt = *sIt - *(sIt - 1);
    }

    // Not-a-knot at 0
    const auto s0 = s[0];
    const auto s1 = s[1];
    const auto z1 = z[1];
    const auto z2 = z[2];
    z[0] = (6. * (s1 - s0) - 2. * (m_h0 + m_h1) * z1 - m_h1 * z2) * m_g0;

    // Not-a-knot at m_size - 1
    // FIXME Should it be at m_size?
    const auto sm0 = s[m_size - 1];
    const auto sm1 = s[m_size - 2];
    const auto zm1 = z[m_size - 1];
    z[m_size - 1] = (6. * m_hm0 * (sm0 - sm1) - (m_hm02 - m_hm12) * zm1) * m_zm0Factor;

    return z;
  }

  /**
   * @brief Interpolate.
   */
  template <typename T>
  std::vector<T> interpolate(const T* v, const T* z) const {
    std::vector<T> s(m_x.size());
    std::size_t j = 0;
    for (std::size_t i = 0; i < m_size; ++i) {
      for (long n = 0; n < m_n[i]; ++n) {
        s[j] = z[i + 1] * m_k[j].z1 + z[i] * m_k[j].z0 + v[i + 1] * m_k[j].v1 + v[i] * m_k[j].v0;
        ++j;
      }
    }
    return s;
  }

  /**
   * @brief Interpolate.
   */
  template <typename T>
  std::vector<T> operator()(const T* v) const {
    const auto z = knotSecondDerivatives(v);
    return interpolate(v, z.data());
  }

  /**
   * @brief Interpolate and integrate.
   * @param v The values at knots
   * @param z The second derivatives at knots
   * @param w The weights
   */
  template <typename T>
  T integrate(const T* v, const T* z, const double* w) const {
    T sum {};
    auto* vIt = v;
    auto* zIt = z;
    auto* wIt = w;
    auto kIt = m_k.begin();
    for (auto i = m_size; i--; ++vIt, ++zIt) {
      for (auto n = m_n[i]; n--; ++wIt, ++kIt) {
        // sum += w[j] * (z[i + 1] * m_k[i].z1 + z[i] * m_k[i].z0 + v[i + 1] * m_k[i].v1 + v[i] * m_k[i].v0);
        sum += *wIt * (*(zIt + 1) * kIt->z1 + *zIt * kIt->z0 + *(vIt + 1) * kIt->v1 + *vIt * kIt->v0);
      }
    }
    return sum;
  }

private:
  /**
   * @brief The coefficients of a polynom.
   */
  struct Coefficients {
    double z1; ///< Coefficient of z[i + 1]
    double z0; ///< Coefficient of z[i]
    double v1; ///< Coefficient of y[i + 1]
    double v0; ///< Coefficient of y[i]
  };

  /**
   * @brief The number of splines.
   */
  std::size_t m_size;

  /**
   * @brief The sample x's.
   */
  std::vector<double> m_u;

  /**
   * @brief The sample steps.
   */
  std::vector<double> m_h;
  double m_h0; ///< m_h[0]
  double m_h1; ///< m_h[1]
  double m_h02; ///< m_h0 * m_h0
  double m_h12; ///< m_h1 * m_h1
  double m_h01; ///< m_h0 * m_h1
  double m_hm0; ///< m_h[m_size - 1]
  double m_hm1; ///< m_h[m_size - 2]
  double m_hm02; ///< m_hm0 * m_hm0
  double m_hm12; ///< m_hm1 * m_hm1
  double m_hm01; ///< m_hm0 * m_hm1
  std::vector<double> m_g; ///< Inverses of m_h[i]'s
  double m_g0; ///< m_g[0]
  double m_zm0Factor; ///< Normalization factor for z[m_size - 1]

  /**
   * @brief The integration x's.
   */
  std::vector<double> m_x;

  /**
   * @brief The number of integration x's per spline.
   */
  std::vector<long> m_n;

  /**
   * @brief The pre-computed constants.
   */
  std::vector<Coefficients> m_k;
};

} // namespace Challenger

#endif