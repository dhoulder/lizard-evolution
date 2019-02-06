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
#include "model.h"

using namespace DreadDs;

int main(int argc, const char *argv[]) {
  try {
    Config conf(argc, argv);
    Model model(conf);
    for (int i =0; i < conf.n_iterations; ++i) {
      int n = model.do_step();
      model.save();
      if (n == 0)
	break;
    }
    if (conf.verbosity > 0)
      std::cout << "Final step count: " << model.step << std::endl;
  }
  catch (const UsageException &ae) {
    std::cerr <<
      "Usage: " << argv[0] << "--option=value [--option=value …]\n"  <<
      ae.what()  << std::endl;
    std::exit(0);
  }
  catch (const ConfigError &ae) {
    std::cerr << ae.what() <<
      "\nSee " << argv[0] << " --help for usage." << std::endl;
    std::exit(1);
  }
  catch (const ApplicationException &ae) {
    std::cerr << ae.what()  << std::endl;
    std::exit(2);
  }
}
