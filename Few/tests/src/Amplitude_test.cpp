/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Few/Amplitude.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Amplitude_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(sum_test) {
  std::cout << "Creating carrier..." << std::endl;
  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  Few::AmplitudeCarrier carrier(lmax, nmax, "");
  std::cout << "  done." << std::endl;

  // FIXME many values
  std::cout << "Interpolating..." << std::endl;
  std::vector<double> ps {4.5};
  std::vector<double> es {14.5};
  std::vector<int> ls {5};
  std::vector<int> ms {5};
  std::vector<int> ns {5};
  const auto n = ps.size();
  const auto nmodes = ls.size();
  std::vector<std::complex<double>> out(n * nmodes);

  carrier.Interp2DAmplitude(out.data(), ps.data(), es.data(), ls.data(), ms.data(), ns.data(), n, nmodes);

  std::cout << "  ";
  for (const auto& e : out) {
    std::cout << e << " ";
  }
  std::cout << std::endl;
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
