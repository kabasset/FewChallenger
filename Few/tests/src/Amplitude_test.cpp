/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Few/Amplitude.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Amplitude_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(sum_test) {
  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  Few::AmplitudeCarrier carrier(lmax, nmax, "");

  // FIXME many values
  std::vector<double> ps {1.};
  std::vector<double> es {1.};
  std::vector<int> ls {5};
  std::vector<int> ms {5};
  std::vector<int> ns {5};
  const auto n = ps.size();
  const auto nmodes = ls.size();
  std::vector<std::complex<double>> out(n * nmodes);

  carrier.Interp2DAmplitude(out.data(), ps.data(), es.data(), ls.data(), ms.data(), ns.data(), n, nmodes);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()