/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MemoryChallenger/Interpolant.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Interpolant_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(real_test) {
  const std::vector<double> u0 {1, 2, 3, 4};
  const std::vector<double> u1 {1, 10, 100, 1000};
  const std::vector<Linx::Vector<double, 2>> x {{1.1, 2.}, {2.5, 20.}, {3.9, 50.}};
  const Linx::Raster<double, 2> v(
      {u0.size(), u1.size()},
      {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000});
  MemoryChallenger::Interpolant2D interpolant(u0, u1, v);
  for (std::size_t i = 0; i < x.size(); ++i) {
    const auto y = interpolant(x[i][0], x[i][1]);
    std::cout << x[i] << " - " << y << "\n";
    Linx::Position<2> p;
    for (std::size_t j = 1; j < u0.size(); ++j) {
      if (u0[j] < x[i][0]) {
        p[0] = j;
      }
    }
    for (std::size_t j = 1; j < u1.size(); ++j) {
      if (u1[j] < x[i][1]) {
        p[1] = j;
      }
    }
    BOOST_TEST(y > v[p]);
    BOOST_TEST(y < v[p + 1]);
  }
}

BOOST_AUTO_TEST_CASE(complex_test) {
  const std::vector<double> x {1, 2, 3, 4};
  const std::vector<double> y {1, 10, 100, 1000};
  const Linx::Raster<double, 3> z(
      {x.size(), y.size(), 2},
      {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000,
       0, 0, 0, 0, 0,  0,  0,  0,  0,   0,   0,   0,   0,    0,    0,    0});
  MemoryChallenger::ComplexInterpolant2D interpolant(x, y, z);
  const auto out = interpolant(2.5, 25);
  BOOST_TEST(out.real() > 20);
  BOOST_TEST(out.real() < 300);
  BOOST_TEST(out.imag() == 0);
}

BOOST_AUTO_TEST_CASE(raster_test) {
  const std::vector<double> x {1, 2, 3, 4};
  const std::vector<double> y {1, 10, 100, 1000};

  Linx::AlignedRaster<MemoryChallenger::ComplexInterpolant2D, 3> interpolants({2, 2, 2}, nullptr, 0);
  for (auto& e : interpolants) {
    const Linx::Raster<double, 3> z(
        {x.size(), y.size(), 2},
        {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000,
         0, 0, 0, 0, 0,  0,  0,  0,  0,   0,   0,   0,   0,    0,    0,    0});
    e = MemoryChallenger::ComplexInterpolant2D(x, y, z);
  }

  for (const auto& e : interpolants) {
    const auto out = e(2.5, 25);
    BOOST_TEST(out.real() > 20);
    BOOST_TEST(out.real() < 300);
    BOOST_TEST(out.imag() == 0);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
