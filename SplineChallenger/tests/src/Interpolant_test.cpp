/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SplineChallenger/Interpolant.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Interpolant_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(oriented_real_test) {
  const std::vector<double> u {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  SplineChallenger::OrientedSpline spline(u, x);
  const std::vector<double> v {10, 20, 30, 40};
  const auto y = spline(v.data());
  BOOST_TEST(y.size() == x.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    BOOST_TEST(y[i] > v[i]);
    BOOST_TEST(y[i] < v[i + 1]);
  }
}

BOOST_AUTO_TEST_CASE(oriented_complex_test) {
  const std::vector<double> u {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  SplineChallenger::OrientedSpline spline(u, x);
  const std::vector<std::complex<double>> v {{10, -1}, {20, -2}, {30, -3}, {40, -4}};
  const auto y = spline(v.data());
  BOOST_TEST(y.size() == x.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    BOOST_TEST(y[i].real() > v[i].real());
    BOOST_TEST(y[i].real() < v[i + 1].real());
    BOOST_TEST(y[i].imag() < v[i].imag());
    BOOST_TEST(y[i].imag() > v[i + 1].imag());
  }
}

BOOST_AUTO_TEST_CASE(separable_real_test) {
  const std::vector<double> u0 {1, 2, 3, 4};
  const std::vector<double> u1 {1, 10, 100};
  const std::vector<double> x0 {1.1, 2.5, 3.9};
  const std::vector<double> x1 {2., 50.};
  SplineChallenger::SeparableSpline<2> spline({u0, u1}, {x0, x1});
  const Linx::Raster<double> v({u0.size(), u1.size()}, {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400});
  const auto y = spline(v);
  BOOST_TEST(y.size() == x0.size() * x1.size());
  for (const auto& p : y.domain()) {
    BOOST_TEST(y[p] > v[p]);
    BOOST_TEST(y[p] < v[p + 1]);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
