// -*- coding: utf-8 -*-

#include <fstream>

#include "boost/multi_array.hpp"

#include "model-limits.h"
#include "model-config.h"
#include "environment.h"

using namespace DreadDs;


static void process_header(Environment *env,
			   std::ifstream &ifs) {
  int rows = 5;
  int cols = 6; // FIXME STUB


  // first file determines
  env->values.resize(boost::extents[rows][cols]);
}

Environment::Environment(const EnvParamsVec &env_inputs) {

  bool first = true;
  for (const auto &ep: env_inputs) {
    // FIXME WIP load cells from file

    // first file defines bb. subsequent files must match
    // first file creates env and v[0], subsequent files fill in v[1], v[2], ...


    std::ifstream ifs(ep.grid_filename);
    if(ifs.fail())
      throw ConfigError("Can't open environment input file " + ep.grid_filename);

    if (first) {
      first = false;
      process_header(this, ifs);


    } else {
      // FIXME WIP
    }
  }

  // FIXME WIP fail if empty (first still true)

}
