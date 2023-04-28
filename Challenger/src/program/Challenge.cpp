/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Challenger/Challenger.h"
#include "ElementsKernel/ProgramHeaders.h"
#include "LinxRun/Chronometer.h"
#include "LinxRun/ProgramOptions.h"

#include <map>
#include <string>

using Duration = std::chrono::milliseconds;

std::complex<double> runFew(int size) {

  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  Few::AmplitudeCarrier carrier(lmax, nmax, "");

  std::vector<double> ps(size, 4.5);
  std::vector<double> es(size, 14.5);
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

  return std::accumulate(out.begin(), out.end(), std::complex<double> {});
}

std::complex<double> runMemoryChallenger(int size) {

  static constexpr int lmax = 10;
  static constexpr int nmax = 30;
  MemoryChallenger::AmplitudeCarrier carrier(lmax, nmax);

  std::vector<MemoryChallenger::Pe> pes(size, MemoryChallenger::Pe {4.5, 14.5});
  std::vector<MemoryChallenger::Lmn> lmns;
  for (Linx::Index l = 2; l <= lmax; ++l) {
    for (Linx::Index m = 0; m <= l; ++m) {
      for (Linx::Index n = -nmax; n <= nmax; ++n) {
        lmns.push_back({l, m, n});
      }
    }
  }

  const auto out = carrier.interpolate(pes, lmns);

  return std::accumulate(out.begin(), out.end(), std::complex<double> {});
}

class Challenge : public Elements::Program {

public:
  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override {
    Linx::ProgramOptions options;
    options.named("case", "Test case: f (Few), c (Challenger)", 'f');
    options.named("size", "Number of trajectory points", 1);
    return options.asPair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {

    Logging logger = Logging::getLogger("Challenge");
    const int size = args["size"].as<int>();

    logger.info() << "Running benchmark...";
    Linx::Chronometer<Duration> chrono;
    chrono.start();
    switch (args["case"].as<char>()) {
      case 'f':
        runFew(size);
        break;
      case 'm':
        runMemoryChallenger(size);
        break;
      default:
        logger.error("Unknown test case.");
        return ExitCode::NOT_OK;
    }
    const auto duration = chrono.stop();

    logger.info() << "  Done in " << duration.count() << "ms";

    return ExitCode::OK;
  }
};

MAIN_FOR(Challenge)
