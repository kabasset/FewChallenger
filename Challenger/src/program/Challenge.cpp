/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Challenger/Challenger.h"
#include "ElementsKernel/ProgramHeaders.h"
#include "LinxRun/Chronometer.h"
#include "LinxRun/ProgramOptions.h"

#include <map>
#include <string>

using Duration = std::chrono::milliseconds;

std::vector<std::complex<double>> runFew(const std::vector<double>& ps, const std::vector<double>& es) {

  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  Few::AmplitudeCarrier carrier(lmax, nmax, "");

  std::vector<int> ls;
  std::vector<int> ms;
  std::vector<int> ns;
  for (int l = 2; l <= lmax; ++l) {
    for (int m = 0; m <= l; ++m) {
      for (int n = -nmax; n <= nmax; ++n) {
        ls.push_back(l);
        ms.push_back(m);
        ns.push_back(n);
      }
    }
  }
  const auto n = ps.size();
  const auto nmodes = ls.size();
  std::vector<std::complex<double>> out(n * nmodes);

  carrier.Interp2DAmplitude(out.data(), ps.data(), es.data(), ls.data(), ms.data(), ns.data(), n, nmodes);

  return out;
}

std::vector<std::complex<double>> runMemoryChallenger(const std::vector<MemoryChallenger::Ellipse>& ellipses) {

  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  MemoryChallenger::AmplitudeCarrier carrier(lmax, nmax);

  std::vector<MemoryChallenger::Mode> modes;
  for (Linx::Index l = 2; l <= lmax; ++l) {
    for (Linx::Index m = 0; m <= l; ++m) {
      for (Linx::Index n = -nmax; n <= nmax; ++n) {
        modes.push_back({l, m, n});
      }
    }
  }

  return carrier.interpolate(ellipses, modes).container();
}

std::vector<std::complex<double>> runSplineChallenger(const std::vector<MemoryChallenger::Ellipse>& ellipses) {

  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  std::vector<Linx::Vector<double, 2>> x(ellipses.size());
  for (std::size_t i = 0; i < ellipses.size(); ++i) {
    x[i] = {ellipses[i].y(), ellipses[i].e()};
  }
  const auto u = SplineChallenger::loadGrid();
  const Splider::SplineIntervals u0(u[0]);
  const Splider::SplineIntervals u1(u[1]);
  Splider::BiSplineResampler<std::complex<double>> resample(u0, u1, x);

  std::vector<std::complex<double>> out;
  for (Linx::Index l = 2; l <= lmax; ++l) {
    for (Linx::Index m = 0; m <= l; ++m) {
      for (Linx::Index n = -nmax; n <= nmax; ++n) {
        const auto modal = resample(SplineChallenger::loadModeData({l, m, n}));
        out.insert(out.end(), modal.begin(), modal.end());
      }
    }
  }

  return out;
}

class Challenge : public Elements::Program {

public:
  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override {
    Linx::ProgramOptions options;
    options.named("case", "Test case: f (Few), m (MemoryChallenger), s (SplineChallenger)", 'f');
    options.named("size", "Number of trajectory points", 1);
    return options.asPair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {

    Logging logger = Logging::getLogger("Challenge");
    const int size = args["size"].as<int>();

    std::vector<double> ps(size);
    std::vector<double> es(size);
    std::vector<MemoryChallenger::Ellipse> ellipses(size);
    for (std::size_t i = 0; i < size; ++i) {
      ps[i] = 490. / size * (i + 1) + 10;
      es[i] = 32. / size * i;
      ellipses[i] = {ps[i], es[i]};
      logger.debug() << ps[i] << ", " << es[i] << ", " << ellipses[i].y();
    }

    logger.info() << "Running benchmark...";
    Linx::Chronometer<Duration> chrono;
    std::vector<std::complex<double>> out;
    chrono.start();
    switch (args["case"].as<char>()) {
      case 'f':
        out = runFew(ps, es);
        break;
      case 'm':
        out = runMemoryChallenger(ellipses);
        break;
      case 's':
        out = runSplineChallenger(ellipses);
        break;
      default:
        logger.error("Unknown test case.");
        return ExitCode::NOT_OK;
    }
    const auto duration = chrono.stop();

    logger.info() << "  " << out[0] << " ... " << out[out.size() - 1];
    logger.info() << "  Done in " << duration.count() << "ms";

    return ExitCode::OK;
  }
};

MAIN_FOR(Challenge)
