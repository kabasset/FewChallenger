/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Few/Amplitude.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Interpolant_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(sum_test) {
  static constexpr int lmax = 10;
  static constexpr int nmax = 10; // FIXME
  Few::AmplitudeCarrier carrier(lmax, nmax, "");
  std::vector<std::complex<double>> out; // FIXME size

  // FIXME
  std::vector<double> p_arr;
  std::vector<double> e_arr;
  std::vector<int> l_arr;
  std::vector<int> m_arr;
  std::vector<int> n_arr;
  int num;
  int num_modes;

  carrier.Interp2DAmplitude(
      out.data(),
      p_arr.data(),
      e_arr.data(),
      l_arr.data(),
      m_arr.data(),
      n_arr.data(),
      num,
      num_modes);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
