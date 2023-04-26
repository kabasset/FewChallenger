/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Challenger/Amplitude.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Amplitude_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(singleton_test) {
  std::cout << "Creating carrier..." << std::endl;
  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  Challenger::AmplitudeCarrier carrier(lmax, nmax);
  std::cout << "  done." << std::endl;

  // FIXME many values
  std::cout << "Interpolating..." << std::endl;
  std::vector<Challenger::Pe> pes {Challenger::Pe {4.5, 14.5}};
  std::vector<Challenger::Lmn> lmns {Challenger::Lmn {5, 5, 5}};

  const auto out = carrier.interpolate(pes, lmns);

  std::cout << "  " << out << std::endl;
}

BOOST_AUTO_TEST_CASE(full_lmn_test) {
  std::cout << "Creating carrier..." << std::endl;
  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  Challenger::AmplitudeCarrier carrier(lmax, nmax);
  std::cout << "  done." << std::endl;

  // FIXME many values
  std::cout << "Interpolating..." << std::endl;
  std::vector<Challenger::Pe> pes {Challenger::Pe {4.5, 14.5}};
  std::vector<Challenger::Lmn> lmns;
  for (Linx::Index l = 2; l <= lmax; ++l) { // Idem
    for (Linx::Index m = 0; m <= l; ++m) {
      for (Linx::Index n = -nmax; n <= nmax; ++n) {
        lmns.push_back({l, m, n});
      }
    }
  }
  const auto out = carrier.interpolate(pes, lmns);

  std::cout << "  " << out << std::endl;
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
