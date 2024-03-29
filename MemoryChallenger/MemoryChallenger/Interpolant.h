/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _MEMORYCHALLENGER_INTERPOLANT_H
#define _MEMORYCHALLENGER_INTERPOLANT_H

#include "Linx/Data/Raster.h"

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <vector>

namespace MemoryChallenger {

class Interpolant1D {

public:
  template <typename TSeq>
  Interpolant1D(const TSeq& x, TSeq& y) :
      m_xacc {gsl_interp_accel_alloc()}, m_spline {gsl_spline_alloc(gsl_interp_cspline, x.size())} {
    gsl_spline_init(m_spline, x.data(), y.data(), x.size());
  }

  ~Interpolant1D() {
    gsl_spline_free(m_spline);
    gsl_interp_accel_free(m_xacc);
  }

  double operator()(double x) const {
    return gsl_spline_eval(m_spline, x, m_xacc);
  }

  bool operator==(const Interpolant1D& rhs) const {
    if (this == &rhs) {
      return true;
    }
    return this->m_spline == rhs.m_spline;
  }

private:
  gsl_interp_accel* m_xacc;
  gsl_spline* m_spline;
};

class Interpolant2D {

public:
  template <typename TSeq, typename TMap>
  Interpolant2D(const TSeq& x, const TSeq& y, const TMap& z) :
      m_xacc {gsl_interp_accel_alloc()}, m_yacc {gsl_interp_accel_alloc()},
      m_spline {gsl_spline2d_alloc(gsl_interp2d_bicubic, x.size(), y.size())} {
    gsl_spline2d_init(m_spline, x.data(), y.data(), z.data(), x.size(), y.size());
  }

  Interpolant2D(Interpolant2D&& rhs) : m_xacc(rhs.m_xacc), m_yacc(rhs.m_yacc), m_spline(rhs.m_spline) {
    rhs.m_xacc = nullptr;
    rhs.m_yacc = nullptr;
    rhs.m_spline = nullptr;
  }

  Interpolant2D& operator=(Interpolant2D&& rhs) {
    m_xacc = rhs.m_xacc;
    m_yacc = rhs.m_yacc;
    m_spline = rhs.m_spline;
    rhs.m_xacc = nullptr;
    rhs.m_yacc = nullptr;
    rhs.m_spline = nullptr;
    return *this;
  }

  ~Interpolant2D() {
    gsl_spline2d_free(m_spline);
    gsl_interp_accel_free(m_xacc);
    gsl_interp_accel_free(m_yacc);
  }

  double operator()(double x, double y) const {
    return gsl_spline2d_eval(m_spline, x, y, m_xacc, m_yacc);
  }

  bool operator==(const Interpolant2D& rhs) const {
    if (this == &rhs) {
      return true;
    }
    return this->m_spline == rhs.m_spline;
  }

private:
  gsl_interp_accel* m_xacc;
  gsl_interp_accel* m_yacc;
  gsl_spline2d* m_spline;
};

class ComplexInterpolant2D {

public:
  template <typename TSeq, typename TMap>
  ComplexInterpolant2D(const TSeq& x, const TSeq& y, const TMap& z) :
      m_real(x, y, z.section(0)), m_imag(x, y, z.section(1)) {}

  std::complex<double> operator()(double x, double y) const {
    return std::complex<double>(m_real(x, y), m_imag(x, y));
  }

  bool operator==(const ComplexInterpolant2D& rhs) const {
    if (this == &rhs) {
      return true;
    }
    return this->m_real == rhs.m_real && this->m_imag == rhs.m_imag;
  }

private:
  Interpolant2D m_real;
  Interpolant2D m_imag;
};

} // namespace MemoryChallenger

#endif // _CHALLENGER_INTERPOLANT_H
