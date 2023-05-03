/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MemoryChallenger/Amplitude.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Amplitude_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(singleton_test) {
  std::cout << "Creating carrier..." << std::endl;
  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  MemoryChallenger::AmplitudeCarrier carrier(lmax, nmax);
  std::cout << "  done." << std::endl;

  std::cout << "Interpolating..." << std::endl;
  std::vector<MemoryChallenger::Ellipse> ellipses {MemoryChallenger::Ellipse {10., 0.}};
  std::vector<MemoryChallenger::Mode> modes {MemoryChallenger::Mode {5, 5, 5}};

  const auto out = carrier.interpolate(ellipses, modes);

  std::cout << "  " << out << std::endl;
}

BOOST_AUTO_TEST_CASE(full_lmn_test) {
  std::cout << "Creating carrier..." << std::endl;
  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  MemoryChallenger::AmplitudeCarrier carrier(lmax, nmax);
  std::cout << "  done." << std::endl;

  std::cout << "Interpolating..." << std::endl;
  std::vector<MemoryChallenger::Ellipse> ellipses {MemoryChallenger::Ellipse {10., 0.}};
  std::vector<MemoryChallenger::Mode> modes;
  for (Linx::Index l = 2; l <= lmax; ++l) { // Idem
    for (Linx::Index m = 0; m <= l; ++m) {
      for (Linx::Index n = -nmax; n <= nmax; ++n) {
        modes.push_back({l, m, n});
      }
    }
  }
  const auto out = carrier.interpolate(ellipses, modes);

  std::cout << "  " << out << std::endl;
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
