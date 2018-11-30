// -*- coding: utf-8 -*-

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

#include "model-config.h"
#include "exceptions.h"

namespace po = boost::program_options;
namespace bfs = boost::filesystem;

using namespace std;
using namespace DreadDs;


static vector<float> comma_split(string csv) {
  vector<float> v;
  boost::tokenizer<boost::escaped_list_separator<char>> tok(csv);
  for (const auto &t : tok)
    v.push_back(stof(t));
  return v;
}


Config::Config(int ac, const char *av[]) {
  /**
   * Configuration values are obtained using Boost.Program_options See
   * https://www.boost.org/doc/libs/1_68_0/doc/html/program_options.html
   * These take the form of "$key = $value " strings. For sets of
   * configuration options that can be specified multiple times, the
   * individual "key = value" pairs are repeated as many times as
   * necessary. See
   * https://www.boost.org/doc/libs/1_68_0/doc/html/program_options/overview.html#id-1.3.31.5.10.2
   * for config file format

   */

  try {
    vector <double> env_ramp;
    vector <string> env_grids, env_mode, env_ts_dirs,
      species_niche, species_niche_breadth;
    vector <int> env_ts_start, env_ts_step;
    vector <float>
      env_sine_period, env_sine_offset, env_sine_amplitude,
      species_north, species_east, species_south, species_west,
      species_max_dispersal_radius;

    po::options_description cmdline_options("Options available only on the command line");
    cmdline_options.add_options()
            ("help,h", "Show this help message.")
            ("config,c", po::value<string>(),
	     "Pathname of configuration file containing options described below. "
	     "Command-line arguments take precedence. Arguments of the form "
	     "section.name=value can be shortened to name=value under a [section] line. "
	     "Use the # character to introduce comments.");

    po::options_description config_options("Configuration options");
    config_options.add_options()
      ("verbosity,v", po::value<int>(&verbosity)->default_value(1),
       "Verbosity level, 0 (quiet) to 3.")

      ("genetic-dims", po::value<int>(&genetic_dims)->default_value(max_genetic_dims),
       "Number of genetic dimensions")

      ("gene-flow-threshold", po::value<float>(&gene_flow_clip)->required(),
       "Gene flow probabilities that are closer to 0 or 1 by this amount or "
       "less are treated as 0 or 1 respectively")

      ("gene-flow-max-distance", po::value<float>(&gene_flow_zero_distance)->required(),
       "Gene flow probability approaches zero for genetic distances greater than this")

      ("gene-drift-rate", po::value<float>(&gene_drift_sd)->required(),
       "Standard deviation of Brownian motion gene drift per time-step")

      ("niche-evolution-rate", po::value<float>(&niche_evolution_rate)->required(),
       "Speed of adaptation to niche as a fraction of the "
       "distance to niche centre per time-step")

      ("dispersal-min", po::value<float>(&dispersal_min)->required(),
       "Minimum dispersal amount (0…1)")

      ("output-dir,o",  po::value<string>(&output_dir)->required(),
       "Pathname of directory for output files.")

      ("output-file-prefix,p",  po::value<string>(&output_file_prefix),
       "Output filenames will all start with this string.")

      ("iterations,n", po::value<int>(&n_iterations)->required(),
       "Number of time steps to simulate.");

    po::options_description env_options("Input environment (one or more sets)");
    env_options.add_options()
      ("env.mode", po::value<vector<string>>(&env_mode),
       "ts (time-series) or ex (extrapolate). For ts, a set of input files "
       "containing the environment at different time steps is supplied. "
       "In ex mode, a base environment file is supplied, and the environment "
       "for each time step is extrapolated using the ramp and sine options below.")

      ("env.grid", po::value<vector<string>>(&env_grids)->required(),
       "Pathname of grid file containing input environment values. "
       "When used in conjunction with env.mode=ts, this path "
       "is relative to each numbered sub-directory of ts-dir. "
       "Can be specified multiple times for multiple environment dimensions.")

      ("env.ts-dir", po::value<vector<string>>(&env_ts_dirs),
       "(ts mode) Pathname of a directory containing numerically named "
       "sub-directories containing input environment data ")

      ("env.start-dir", po::value<vector<int>>(&env_ts_start),
       "(ts mode) Specifies the number of the numerically named "
       "sub-directory of ts-dir to use for the first time step.")

      ("env.step-by", po::value<vector<int>>(&env_ts_step),
       "(ts mode) Specifies the increment to use for each successive "
       "time step. Can be negative.")

      ("env.ramp", po::value<vector<double>>(&env_ramp),
       "(ex mode) Environmental change per time step.")

      ("env.sine-period", po::value<vector<float>>(&env_sine_period),
       "(ex mode) Length of sinusoidal environmental change in time steps.")

      ("env.sine-offset", po::value<vector<float>>(&env_sine_offset),
       "(ex mode) Offset of sinusoidal environmental change in time-steps")

      ("env.sine-amplitude", po::value<vector<float>>(&env_sine_amplitude),
       "(ex mode) Peak amplitude of sinusoidal environmental change.");


    po::options_description species_options("Initial species (one or more sets)");
    species_options.add_options()
      ("species.niche", po::value<vector<string>>(&species_niche)->required(),
       "Species niche as comma-separated values, one per environment dimension")

      ("species.niche-breadth",
       po::value<vector<string>>(&species_niche_breadth)->required(),
       "Species niche breadth")

      ("species.max-dispersal-radius",
       po::value<vector<float>>(&species_max_dispersal_radius)->required(),
       "Maximum dispersal radius in grid cells")

      ("species.north", po::value<vector<float>>(&species_north)->required(),
       "Northern boundary of initial species bounding rectangle in grid coordinates")

      ("species.east", po::value<vector<float>>(&species_east)->required(),
       "Eastern boundary of initial species bounding rectangle in grid coordinates")

      ("species.south", po::value<vector<float>>(&species_south)->required(),
       "Southern boundary of initial species bounding rectangle in grid coordinates")

      ("species.west", po::value<vector<float>>(&species_west)->required(),
       "Western boundary of initial species bounding rectangle in grid coordinates");

    config_options.add(env_options).add(species_options);
    cmdline_options.add(config_options);

    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).
              options(cmdline_options).run(), vm);

    if (vm.count("help")) {
      std::ostringstream os;
      os << cmdline_options;
      throw UsageException(os.str());
    }

    if (vm.count("config")) {
      string config_file = vm["config"].as<string>();
      if (!config_file.empty()) {
	ifstream ifs(config_file);
	if(ifs.fail())
	  throw ConfigError("Can't open configuration file " +
			    config_file);
	po::store(po::parse_config_file(ifs, config_options), vm);
      }
    }
    po::notify(vm);

    env_dims = env_grids.size();
    if (env_dims > max_env_dims)
      throw ConfigError("Too many environment grids");
    if (env_dims < 1)
      throw ConfigError("No environment grids");

    if (dispersal_min < 0.0 || dispersal_min > 1.0)
      throw ConfigError("dispersal-min must be between 0 and 1");

    if (env_mode.size() != env_dims)
      throw ConfigError("Each env.grid needs a corresponding env.mode");

    int ts_i =0, ex_i = 0;
    std::string mode;
    for (int i = 0; i < env_dims; ++i) {
      EnvParams ep;
      ep.grid_filename = env_grids[i];
      try {
	mode = env_mode[i];
      }
      catch (const out_of_range& oor) {
	throw ConfigError("Each env.grid needs a corresponding env.mode");
      }
      if (mode == "ts") {
	try {
	  ep.ts_dir = env_ts_dirs.at(ts_i);
	  ep.ts_start = env_ts_start.at(ts_i);
	  ep.ts_step = env_ts_step.at(ts_i);
	}
	catch (const out_of_range& oor) {
	  throw ConfigError("For env.mode=ts, each env.grid needs a corresponding "
			    "env.ts-dir, env.start-dir, env.step-by.");
	}
	ep.scan_ts_dir();
	++ts_i;

      } else if (mode == "ex")  {
	// base + offset mode
	try {
	  ep.ramp = env_ramp.at(ex_i);
	  ep.sine_period = env_sine_period.at(ex_i);
	  ep.sine_offset = env_sine_offset.at(ex_i);
	  ep.sine_amplitude = env_sine_amplitude.at(ex_i);
	}
	catch (const out_of_range& oor) {
	  throw ConfigError("For env.mode=ex, each env.grid needs a corresponding "
			    "env.ramp, env.sine-period, env.sine-offset, "
			    "env.sine-amplitude");
	}
	++ex_i;
      }
      else throw ConfigError("env.mode must be ex or ts.");

      env_params.push_back(ep);
    }

    try {

      for (int i = 0; i < species_niche.size(); ++i) {
	SpeciesParameters sp;
	sp.north = species_north.at(i);
	sp.south = species_south.at(i);
	sp.east = species_east.at(i);
	sp.west = species_west.at(i);
	sp.max_dispersal_radius = species_max_dispersal_radius.at(i);
	if (sp.max_dispersal_radius < 0.0)
	  throw ConfigError("species.max-dispersal-radius must not be negative");

	// Need species niche for each environment dimension, so
	// they're supplied as CSV in the config file
	try {
	  vector<float> nc = move(comma_split(species_niche.at(i)));
	  vector<float> nt = move(comma_split(species_niche_breadth.at(i)));
	  if (nc.size() != nt.size())
	    throw ConfigError("niche and niche-breadth have different lengths");
	  for (int i=0; i < nc.size(); ++i) {
	    sp.niche.push_back(NicheSpec {nc[i],
		  nt[i] * 0.5f }); // *0.5 to convert breadth to tolerance
	  }
	}
	catch (const invalid_argument &e) {
	  throw ConfigError("Bad niche value. Expected a number");
	}

	if (sp.niche.size() < env_dims)
	  throw ConfigError("Not enough niche values. "
			    "Need one for each environment dimension.");

	// FIXME set genetics? (currently all 0.0)
	initial_species.push_back(sp);
      }
    }
    catch (const out_of_range &e) {
      throw ConfigError("Each species.niche needs a corresponding "
			"species.niche-breadth, species.max-dispersal-radius, "
			"species.north, species.east, species.south, species.west");
    }

    if (initial_species.size() < 1)
      throw ConfigError("No initial species");

    if (genetic_dims >max_genetic_dims)
      throw ConfigError("Too many genetic dimensions");

    if (gene_flow_clip < 0.0 or gene_flow_clip > 1.0)
      throw ConfigError("Gene-flow threshold must be between 0 and 1");

    if (gene_flow_zero_distance <= 0.0f)
      throw ConfigError("Gene-flow max distance must be greater than 0");

    if (gene_drift_sd < 0.0f)
      throw ConfigError("Gene-flow drift rate cannot be negative");

    if (niche_evolution_rate < 0.0 or niche_evolution_rate > 1.0)
      throw ConfigError("Niche evolution rate must be between 0 and 1");
  }

  catch(po::error &e)
    {
      throw ConfigError("Configuration error: " + string(e.what()));
    }

}


void EnvParams::scan_ts_dir() {
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

string EnvParams::get_filename(int step_offset) const {
  if (ts_dir.empty())
    return grid_filename;
  else {
    int n = ts_start + step_offset*ts_step;
    try {
      return (ts_dir_table.at(n) / grid_filename).native();
    }
    catch (const out_of_range &e) {
      throw ConfigError("No time-series directory for time " +
			to_string(n) + " in " + ts_dir);
    }
  }
}
