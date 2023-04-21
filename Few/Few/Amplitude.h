/// @copyright 2020, Niels Warburton, Michael L. Katz, Alvin J.K. Chua, Scott A. Hughes
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _FEW_AMPLITUDE_H
#define _FEW_AMPLITUDE_H

#include "Few/Interpolant.h"

#include <complex>

namespace Few {

// The 11 below means the lmax = 10
struct waveform_amps {
  Interpolant*** re[11];
  Interpolant*** im[11];
};

class AmplitudeCarrier {
public:
  struct waveform_amps* amps;
  int lmax, nmax;

  AmplitudeCarrier(int lmax_, int nmax_, std::string few_dir);
  void Interp2DAmplitude(
      std::complex<double>* amplitude_out,
      double* p_arr,
      double* e_arr,
      int* l_arr,
      int* m_arr,
      int* n_arr,
      int num,
      int num_modes);

  void dealloc();
};

} // namespace Few

#endif // _FEW_AMPLITUDE_H
