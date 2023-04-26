/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Challenger/Spline.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Spline_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(oriented_real_test) {
  const std::vector<double> u {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  Challenger::OrientedSpline spline(u, x);
  const std::vector<double> v {10, 20, 30, 40};
  const auto out = spline(v.data());
  BOOST_TEST(out.size() == x.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    BOOST_TEST(out[i] > v[i]);
    BOOST_TEST(out[i] < v[i + 1]);
  }
}

BOOST_AUTO_TEST_CASE(oriented_complex_test) {
  const std::vector<double> u {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  Challenger::OrientedSpline spline(u, x);
  const std::vector<std::complex<double>> v {{10, -1}, {20, -2}, {30, -3}, {40, -4}};
  const auto out = spline(v.data());
  BOOST_TEST(out.size() == x.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    BOOST_TEST(out[i].real() > v[i].real());
    BOOST_TEST(out[i].real() < v[i + 1].real());
    BOOST_TEST(out[i].imag() < v[i].imag());
    BOOST_TEST(out[i].imag() > v[i + 1].imag());
  }
}

BOOST_AUTO_TEST_CASE(separable_real_test) {
  const std::vector<double> u0 {1, 2, 3, 4};
  const std::vector<double> u1 {1, 10, 100};
  const std::vector<double> x0 {1.1, 2.5, 3.9};
  const std::vector<double> x1 {2., 50.};
  Challenger::SeparableSpline<2> spline({u0, u1}, {x0, x1});
  const Linx::Raster<double> v({u0.size(), u1.size()}, {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400});
  const auto out = spline(v);
  BOOST_TEST(out.size() == x0.size() * x1.size());
  // for (const auto& p : out.domain()) {
  //   BOOST_TEST(out[p] > v[p - {0, 1}]);
  //   BOOST_TEST(out[p] > v[p - {1, 0}]);
  //   BOOST_TEST(out[p] > v[p - {1, 1}]);
  //   BOOST_TEST(out[p] < v[p + {0, 1}]);
  //   BOOST_TEST(out[p] < v[p + {1, 0}]);
  //   BOOST_TEST(out[p] < v[p + {1, 1}]);
  // }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
