// -*- coding: utf-8 -*-

/**
 * Dynamic Range Evolution and Diversification simulation with dynamic
 * speciation
 *
 */

#include <iostream>
#include <cstdlib>

#include "exceptions.h"
#include "model-config.h"
#include "simulation.h"

using namespace DreadDs;

int main(int argc, const char *argv[]) {
  try {
    Config conf(argc, argv);
    Simulation sim(conf);
    int s = sim.run(conf.n_iterations);
    if (conf.verbosity > 0)
      std::cout << "Final step count: " << s << std::endl;
  }
  catch (const UsageException &ae) {
    std::cerr <<
      "Usage: " << argv[0] << "--option=value [--option=value …]" << std::endl <<
      ae.what()  << std::endl;
    std::exit(0);
  }
  catch (const ConfigError &ae) {
    std::cerr << ae.what()  << std::endl <<
      "See " << argv[0] << " --help for usage." << std::endl;
    std::exit(1);
  }
  catch (const ApplicationException &ae) {
    std::cerr << ae.what()  << std::endl;
    std::exit(2);
  }
}
