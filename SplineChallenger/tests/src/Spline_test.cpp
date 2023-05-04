/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SplineChallenger/Spline.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Spline_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(oriented_real_test) {
  const Focs::SplineDomain u {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  const std::vector<double> v {10, 20, 30, 40};
  Focs::Spline<double> spline(u, v);
  std::vector<double> y;
  for (const auto& e : x) {
    y.push_back(spline(e));
  }
  BOOST_TEST(y.size() == x.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    BOOST_TEST(y[i] > v[i]);
    BOOST_TEST(y[i] < v[i + 1]);
  }
}

BOOST_AUTO_TEST_CASE(oriented_complex_test) {
  const Focs::SplineDomain u {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  const std::vector<std::complex<double>> v {{10, -1}, {20, -2}, {30, -3}, {40, -4}};
  Focs::Spline<std::complex<double>> spline(u, v);
  std::vector<std::complex<double>> y;
  for (const auto& e : x) {
    y.push_back(spline(e));
  }
  BOOST_TEST(y.size() == x.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    BOOST_TEST(y[i].real() > v[i].real());
    BOOST_TEST(y[i].real() < v[i + 1].real());
    BOOST_TEST(y[i].imag() < v[i].imag());
    BOOST_TEST(y[i].imag() > v[i + 1].imag());
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
