/// @copyright 2020, Niels Warburton, Michael L. Katz, Alvin J.K. Chua, Scott A. Hughes
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _FEW_INTERPOLANT_H
#define _FEW_INTERPOLANT_H

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <vector>

namespace Few {

typedef std::vector<double> Vector;

class Interpolant {
public:
  // 1D interpolation
  Interpolant(Vector x, Vector f);
  double eval(double x);

  // 2D interpolation
  Interpolant(Vector x, Vector y, Vector f);
  double eval(double x, double y);

  // Destructor
  ~Interpolant();

private:
  int interp_type; // Set to 1 for 1D interpolation and 2 for 2D interpolation

  gsl_spline* spline;
  gsl_spline2d* spline2d;
  gsl_interp_accel* xacc;
  gsl_interp_accel* yacc;
};

} // namespace Few

#endif // _FEW_INTERPOLANT_H
