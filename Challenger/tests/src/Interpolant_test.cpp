/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Challenger/Interpolant.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Interpolant_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(real_test) {
  const std::vector<double> x {1, 2, 3, 4};
  const std::vector<double> y {1, 10, 100, 1000};
  const Linx::Raster<double> z(
      {x.size(), y.size()},
      {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000});
  Challenger::Interpolant2D interpolant(x, y, z);
  const auto out = interpolant(2.5, 25);
  BOOST_TEST(out > 20);
  BOOST_TEST(out < 300);
}

BOOST_AUTO_TEST_CASE(complex_test) {
  const std::vector<double> x {1, 2, 3, 4};
  const std::vector<double> y {1, 10, 100, 1000};
  const Linx::Raster<double, 3> z(
      {x.size(), y.size(), 2},
      {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000,
       0, 0, 0, 0, 0,  0,  0,  0,  0,   0,   0,   0,   0,    0,    0,    0});
  Challenger::ComplexInterpolant2D interpolant(x, y, z);
  const auto out = interpolant(2.5, 25);
  BOOST_TEST(out.real() > 20);
  BOOST_TEST(out.real() < 300);
  BOOST_TEST(out.imag() == 0);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
