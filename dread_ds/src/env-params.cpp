// -*- coding: utf-8 -*-

#include <string>
#include <boost/filesystem.hpp>
#include <stdexcept>
#include <math.h>
#include <iostream>

#include "env-params.h"
#include "exceptions.h"
#include "environment.h"

namespace bfs = boost::filesystem;

using namespace std;
using namespace DreadDs;

void TsEnvParams::scan_ts_dir() {
  /**
   * Build a map of directories in which we look for input environment
   * files. This approach allows directories to be named in a fairly
   * flexible manner (e.g. 10, 010, 10-something) and facilitates
   * interpolating between time steps if required.
   */
  try {
    for (auto &&de : bfs::directory_iterator(ts_dir)) {
      auto &&p = de.path();
      int n;
      try {
	n = stoi(p.filename().native());
      }
      catch (invalid_argument &e) {
	// Non-numeric directory name. Ignore it.
	continue;
      }
      if (!bfs::is_regular_file(p / grid_filename)) {
	throw ConfigError("No file called " + grid_filename + " in " + p.native());
      }
      auto &&v = ts_dir_table.emplace(n, p);
      if (!v.second)
	throw ConfigError("Time-series directory named " +
			  v.first->second.native() +
			  " conflicts with  directory named " +
			  p.native());
    }
  }
  catch (bfs::filesystem_error &e) {
    throw ConfigError("Error scanning " + ts_dir + ": " + e.what());
  }
}

string TsEnvParams::get_filename(int step_offset) const {
  int n = ts_start + step_offset*ts_step;
  try {
    return (ts_dir_table.at(n) / grid_filename).native();
  }
  catch (const out_of_range &e) {
    throw ConfigError("No time-series directory for time " +
		      to_string(n) + " in " + ts_dir);
  }
}


float TsEnvParams::update_environment(Environment *env,
				      int step_offset,
				      int layer) const {
  if (env->current_step_offset != step_offset) {
    // this file not already loaded
    EnvReader er(get_filename(step_offset));
    env->check_coordinates(er);
    env->load(er, layer);
  }
  return 0.0f;
}


string ExEnvParams::get_filename(int step_offset) const {
  return grid_filename;
}

float ExEnvParams::update_environment(Environment *env,
				      int step_offset,
				      int layer) const {
  float delta = (step_offset * ramp) +
    (sine_amplitude * sin(2 * M_PI *
			  ((double)sine_offset +
			   (double)step_offset/sine_period)));
  if (env->conf.verbosity > 1)
    std::cout << "Environment variable " << layer <<
      " delta=" << delta << " at step " << step_offset+1 << std::endl;
  return delta;
}
