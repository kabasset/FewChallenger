/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Challenger/Challenger.h"
#include "ElementsKernel/ProgramHeaders.h"
#include "Linx/Run/Chronometer.h"
#include "Linx/Run/ProgramOptions.h"

#include <map>
#include <string>

using Duration = std::chrono::milliseconds;

std::size_t out_size(int lmax, int nmax, std::size_t sample_count) {
  return 3843 * sample_count; // FIXME compute from lmax and nmax
}

std::vector<std::complex<double>> runFew(const std::vector<double>& ps, const std::vector<double>& es) {

  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  Few::AmplitudeCarrier carrier(lmax, nmax, "");

  const auto size = out_size(lmax, nmax, ps.size());
  std::vector<int> ls;
  ls.reserve(size);
  std::vector<int> ms;
  ms.reserve(size);
  std::vector<int> ns;
  ns.reserve(size);

  for (int l = 2; l <= lmax; ++l) {
    for (int m = 0; m <= l; ++m) {
      for (int n = -nmax; n <= nmax; ++n) {
        ls.push_back(l);
        ms.push_back(m);
        ns.push_back(n);
      }
    }
  }
  std::vector<std::complex<double>> out(size);

  carrier.Interp2DAmplitude(out.data(), ps.data(), es.data(), ls.data(), ms.data(), ns.data(), ps.size(), ls.size());

  return out;
}

std::vector<std::complex<double>> runMemoryChallenger(const std::vector<MemoryChallenger::Ellipse>& ellipses) {

  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  MemoryChallenger::AmplitudeCarrier carrier(lmax, nmax);

  const auto size = out_size(lmax, nmax, ellipses.size());
  std::vector<MemoryChallenger::Mode> modes;
  modes.reserve(size);
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
  std::size_t size = out_size(lmax, nmax, ellipses.size());

  std::vector<Linx::Vector<double, 2>> x(ellipses.size());
  for (std::size_t i = 0; i < ellipses.size(); ++i) {
    x[i] = {ellipses[i].y(), ellipses[i].e()};
  }

  const auto u = SplineChallenger::loadGrid();
  const auto build = SplineChallenger::Spline::Multi::builder(u[0], u[1]);
  auto cospline = build.cospline<std::complex<double>>(x);

  std::vector<std::complex<double>> out; // FIXME out(size) (see below)
  out.reserve(size);
  std::vector<std::complex<double>> modal; // FIXME write in out directly (see below)
  modal.reserve(ellipses.size());
  for (Linx::Index l = 2; l <= lmax; ++l) {
    for (Linx::Index m = 0; m <= l; ++m) {
      for (Linx::Index n = -nmax; n <= nmax; ++n) {
        modal = cospline(SplineChallenger::loadModeData({l, m, n})); // FIXME cospline.transform(load(), &out[i])
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
    return options.as_pair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {

    Logging logger = Logging::getLogger("Challenge");
    const int size = args["size"].as<int>();

    logger.info() << "Generating trajectory...";
    std::vector<double> ps(size);
    std::vector<double> es(size);
    std::vector<MemoryChallenger::Ellipse> ellipses(size);
    for (std::size_t i = 0; i < size; ++i) {
      ps[i] = 490. / size * (i + 1) + 10;
      es[i] = 32. / size * i;
      ellipses[i] = {ps[i], es[i]};
      logger.debug() << "  " << ps[i] << ", " << es[i] << ", " << ellipses[i].y();
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

    logger.info() << "  Number of evaluations: " << out.size();
    logger.info() << "  " << out[0] << " ... " << out[out.size() - 1];
    logger.info() << "  Done in: " << duration.count() << "ms";
    for (const auto& e : out) {
      logger.debug() << "  " << e;
    }

    return ExitCode::OK;
  }
};

MAIN_FOR(Challenge)
