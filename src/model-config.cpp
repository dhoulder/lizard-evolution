// -*- coding: utf-8 -*-
/**
 * Load and check configuration options
 */

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <chrono>
#include <mutex>
#include <cmath>
#include <algorithm>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

#include "model-config.h"
#include "exceptions.h"
#include "env-params.h"
#include "environment.h"
#include "random.h"


using namespace std;
using namespace DreadDs;
namespace po = boost::program_options;

typedef std::chrono::seconds Sec;
typedef std::chrono::high_resolution_clock HrClock;
template<class Duration>
using TimePoint = std::chrono::time_point<HrClock, Duration>;

typedef boost::tokenizer<boost::escaped_list_separator<char>> CsvTokenizer;

static vector<float> comma_split(string csv) {
  vector<float> v;
  CsvTokenizer tok(csv);
  for (const auto &t : tok)
    v.push_back(stof(t));
  return v;
}

static mutex rng_mutex;
static bool rng_seeded = false;
static rng_eng_t rng;

Config::Config(int ac, const char *av[]) {
  /**
   * Configuration values are obtained using Boost.Program_options See
   * https://www.boost.org/doc/libs/1_68_0/doc/html/program_options.html
   * These take the form of "$key = $value " strings.
   *
   * For sets of configuration options that can be specified multiple
   * times (e.g. input environment, initial species), the individual
   * "key = value" pairs are repeated as many times as necessary and
   * we rely on ordinal correspondence between those sequences to
   * form the sets.
   *
   * See ../examples/example.conf and
   * https://www.boost.org/doc/libs/1_68_0/doc/html/program_options/overview.html#id-1.3.31.5.10.2
   * for the config file format
   */
  try {
    vector <double> env_ramp;
    vector <string> env_grids, env_mode, env_ts_dirs,
      species_niche, species_niche_breadth, initial_recatangle;
    vector <int> env_ts_start, env_ts_step;
    vector <float>
      env_sine_period, env_sine_offset, env_sine_amplitude,
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

      ("csv-precision", po::value<int>(&csv_precision)->default_value(3),
       "Number of significant digits in numeric CSV field output")

      ("genetic-dims", po::value<int>(&genetic_dims)->default_value(max_genetic_dims),
       "Number of genetic dimensions")

      ("gene-flow-threshold", po::value<float>(&gene_flow_clip)->
       default_value(default_gene_flow_clip),
       "Gene flow probabilities that are closer to 0 or 1 by this amount or "
       "less are treated as 0 or 1 respectively")

      ("gene-flow-max-distance", po::value<float>(&gene_flow_zero_distance)->required(),
       "Gene flow probability approaches zero for genetic distances greater than this")

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
       "Number of time steps to simulate.")

      ("check-speciation",
       po::value<int>(&check_speciation)->default_value(1),
       "Number of timesteps between speciation checks. "
       "Use 0 to disable speciation, 1 to check after every step")
      ;

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
      ("species.niche-centre",
       po::value<vector<string>>(&species_niche)->required(),
       "Species niche centre as comma-separated values, one per environment "
       "dimension, or the string random-cell to derive niche-centre from a random "
       "cell in the landscape. Also see species.initial-rectangle.")

      ("species.niche-breadth",
       po::value<vector<string>>(&species_niche_breadth)->required(),
       "Species niche breadth")

      ("species.max-dispersal-radius",
       po::value<vector<float>>(&species_max_dispersal_radius)->required(),
       "Maximum dispersal radius in grid cell widths")

      ("species.initial-rectangle",
       po::value<vector<string>>(&initial_recatangle)->required(),
       "Initial bounding rectangle for species. Can be specified as either "
       "4 comma-separated grid coordinates (xmin, xmax, ymin, ymax), or (if "
       "niche-centre = random-cell) \"random, <min>, <max>\" to generate a "
       "random starting rectangle with side lengths between <min> and <max> "
       "grid cell widths");

    config_options.add(env_options).add(species_options);
    cmdline_options.add(config_options);

    // Parse command-line options first. This makes them take
    // precedence over options in the the config file.
    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).
              options(cmdline_options).run(), vm);

    if (vm.count("help")) {
      ostringstream os;
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

    if (csv_precision < 0 || csv_precision > 8)
      throw ConfigError("csv-precision must be between 0 and 8");

    env_dims = env_grids.size();
    if (env_dims > max_env_dims)
      throw ConfigError("Too many environment grids");
    if (env_dims < 1)
      throw ConfigError("No environment grids");
    if (env_mode.size() != env_dims)
      throw ConfigError("Each env.grid needs a corresponding env.mode");

    int ts_i =0, ex_i = 0;
    string mode;
    for (int i = 0; i < env_dims; ++i) {
      mode = env_mode[i];

      if (mode == "ts") {
	// Input environment specified as time series, one file per time step.
	auto ep = make_shared<TsEnvParams>(env_grids[i]);
	try {
	  ep->ts_dir = env_ts_dirs.at(ts_i);
	  ep->ts_start = env_ts_start.at(ts_i);
	  ep->ts_step = env_ts_step.at(ts_i);
	}
	catch (const out_of_range& oor) {
	  throw ConfigError("For env.mode=ts, each env.grid needs a corresponding "
			    "env.ts-dir, env.start-dir, env.step-by.");
	}
	ep->scan_ts_dir();
	env_params.push_back(ep);

	++ts_i;

      } else if (mode == "ex")  {
	// Input environment specified as initial environment and
	// extrapolated to each time step using a linear ramp and sine
	// function.
	auto ep = make_shared<ExEnvParams>(env_grids[i]);
	try {
	  ep->ramp = env_ramp.at(ex_i);
	  ep->sine_period = env_sine_period.at(ex_i);
	  ep->sine_offset = env_sine_offset.at(ex_i);
	  ep->sine_amplitude = env_sine_amplitude.at(ex_i);
	}
	catch (const out_of_range& oor) {
	  throw ConfigError("For env.mode=ex, each env.grid needs a corresponding "
			    "env.ramp, env.sine-period, env.sine-offset, "
			    "env.sine-amplitude");
	}
	if (ep->sine_period <= 0.0f)
	  throw ConfigError("sine_period must be greater than 0");
	env_params.push_back(ep);

	++ex_i;
      }
      else throw ConfigError("env.mode must be ex or ts.");
    }

    try {
      for (int i = 0; i < species_niche.size(); ++i) {
	SpeciesParameters sp;
        CsvTokenizer tok(initial_recatangle.at(i));
        vector<string> irs(tok.begin(), tok.end());
        try {
          if (irs.at(0) == "random") {
            sp.bounds_mode = sp.RANDOM_SIZE;
            sp.random_rect_min = stof(irs.at(1));
            sp.random_rect_max = stof(irs.at(2));
            if (sp.random_rect_min > sp.random_rect_max)
              throw ConfigError(
                  "min_length exceeds max_length in "
                  "\"random, <min_length>, <max_length>\"");
          } else {
            sp.bounds_mode = sp.NSEW;
            sp.west = stof(irs.at(0));
            sp.east = stof(irs.at(1));
            sp.south = stof(irs.at(2));
            sp.north = stof(irs.at(3));
            if (sp.west >= sp.east || sp.south >= sp.north)
              throw ConfigError(
                  "Bad initial-rectangle coordinates .Must be "
                  "<xmin>, <xmax>, <ymin>, <ymax>");
          }
        }
        catch (const logic_error &e) {
          throw ConfigError(
              "Bad initial-rectangle. Must be \"<xmin>, <xmax>, <ymin>, <ymax>\" "
              "or \"random, <min_length>, <max_length>\"");
        }

	sp.max_dispersal_radius = species_max_dispersal_radius.at(i);
	if (sp.max_dispersal_radius < 0.0)
	  throw ConfigError("species.max-dispersal-radius must not be negative");

	// Need species niche for each environment dimension, so
	// they're supplied as $value,$value,… in the argument value
	try {
	  vector<float> nc;
	  vector<float> nb = move(comma_split(species_niche_breadth.at(i)));

	  string sn = species_niche.at(i);
	  if (sn == "random-cell") {
	    sp.niche_mode = sp.RANDOM_CELL;
	    nc = move(vector<float>(nb.size(), NAN)); // placeholders
	  } else {
	    sp.niche_mode = sp.NICHE_CENTRE;
	    nc = move(comma_split(sn));
	    if (nc.size() != nb.size())
	      throw ConfigError("niche-centre and niche-breadth have different lengths");
	  }

	  for (int i=0; i < nb.size(); ++i) {
	    sp.niche.push_back(
	      NicheSpec(nc[i],
			nb[i] * 0.5f )); // *0.5 to convert breadth to tolerance
	  }
	}
	catch (const invalid_argument &e) {
	  throw ConfigError("Bad niche value. Expected a number");
	}

	if (sp.niche.size() < env_dims)
	  throw ConfigError("Not enough niche values. Need one for each environment "
			    "dimension. (Separate values with a comma)");

	// FIXME set genetics? (currently all 0.0)
	initial_species.push_back(sp);
      }
    }
    catch (const out_of_range &e) {
      throw ConfigError("Each species.niche-centre needs a corresponding "
			"species.niche-breadth, species.max-dispersal-radius, "
			"species.initial-rectangle");
    }

    if (dispersal_min < 0.0 || dispersal_min > 1.0)
      throw ConfigError("dispersal-min must be between 0 and 1");

    if (initial_species.size() < 1)
      throw ConfigError("No initial species");

    if (genetic_dims >max_genetic_dims)
      throw ConfigError("Too many genetic dimensions");

    if (gene_flow_clip < 0.0 or gene_flow_clip > 1.0)
      throw ConfigError("Gene-flow threshold must be between 0 and 1");

    if (gene_flow_zero_distance <= 0.0f)
      throw ConfigError("Gene-flow max distance must be greater than 0");

    if (niche_evolution_rate < 0.0 or niche_evolution_rate > 1.0)
      throw ConfigError("Niche evolution rate must be between 0 and 1");
  }

  catch(po::error &e)
    {
      throw ConfigError("Configuration error: " + string(e.what()));
    }

  {
    lock_guard<mutex> lock(rng_mutex);
    if (!rng_seeded) {
      rng_seeded = true;
      // seed the RNG with 0…1E9
      auto now = HrClock::now();
      TimePoint<Sec> now_sec = std::chrono::time_point_cast<Sec>(now);
      unsigned int seed = std::chrono::duration_cast<std::chrono::nanoseconds>(
	now-now_sec).count();
      rng.seed(seed);
    }
  }
}

void Config::set_params_from_env(SpeciesParameters &sp,
				const Environment &env) const {
  /**
   * Set species parameters that depend on the input environment.
   */

  auto &&env_shape = env.values.shape();
  long n, s, e, w;   // rows and columns, not coordinates

  float max_y = env.to_ns(0);
  float min_y = env.to_ns(env_shape[0] - 1);
  float min_x = env.to_ew(0);
  float max_x = env.to_ew(env_shape[1] - 1);

  if (sp.bounds_mode == sp.NSEW) {
    sp.north = max(min_y, min(max_y, sp.north));
    sp.south = max(min_y, min(max_y, sp.south));
    sp.east = max(min_x, min(max_x, sp.east));
    sp.west = max(min_x, min(max_x, sp.west));
    n = env.row(sp.north);
    s = env.row(sp.south);
    e = env.col(sp.east);
    w = env.col(sp.west);
    if (n >= s || w >= e)
      throw ConfigError("species.initial-rectangle is not "
			"inside input environment grid");
  } else if (sp.niche_mode != sp.NICHE_CENTRE) {
    n = w = 0;
    s = env_shape[0] - 1;
    e =  env_shape[1] - 1;
  } else
    throw ConfigError("species.initial-rectangle = random "
		      "can only be used with niche-centre = random-cell");

  if (sp.niche_mode == sp.NICHE_CENTRE)
    return;

  // Need to set niche centre from random cell in specified initial
  // rectangle or entire landscape.

  // Find the non-NAN locations within nsew
  std::deque<Location> locations;
  Location loc;
  for (loc.y = n; loc.y <= s; ++loc.y)
    for (loc.x = w; loc.x <= e; ++loc.x) {
      if (env.no_data(loc))
	continue;
      locations.push_back(loc);
    }
  if (locations.size() < 1)
    throw ConfigError("No usable cells within initial rectangle for "
		      "species.niche-centre = random-cell");

  // choose cell and set niche centre to match
  uniform_int_distr_t distr(0, locations.size()-1);
  int random_loc;
  {
    lock_guard<mutex> lock(rng_mutex);
    random_loc = distr(rng);
  }
  const auto &cell_loc = locations[random_loc];
  if (verbosity > 0)
    cout << "Chose random starting cell x=" <<
      env.to_ew(cell_loc.x)  << " (col=" << cell_loc.x <<
      "), y=" << env.to_ns(cell_loc.y) << " (row=" << cell_loc.y << ")" << endl;
  EnvCell ec;
  env.get(cell_loc, ec);
  float const *v = ec;
  for (auto &&ns: sp.niche) {
    ns.centre = *(v++);
    if (verbosity > 0)
      cout << "  niche centre " << v-ec << " at cell = " << ns.centre << endl;
  }

  if (sp.bounds_mode != sp.NSEW) {
    // Choose random area around cell within specified limits.
    uniform_real_distr_t distr(
      // convert side lengths in cell widths to NSEW coordinate
      // offset from chosen cell
      env.geo_transform[1] * sp.random_rect_min * 0.5f,
      env.geo_transform[1] * sp.random_rect_max * 0.5f);
    float lx = env.to_ew(cell_loc.x);
    float ly = env.to_ns(cell_loc.y);
    {
      lock_guard<mutex> lock(rng_mutex);
      // add random NSEW padding around chosen cell, clamped to
      // centres of boundary cells.
      sp.north = min((float)env.to_ns(n),
		     ly +  distr(rng));
      sp.south = max((float)env.to_ns(s),
		     ly - distr(rng));
      sp.east = min((float)env.to_ew(e),
		    lx +  distr(rng));
      sp.west = max((float)env.to_ew(w),
		    lx - distr(rng));
    }
  }
}


std::vector<SpeciesParameters> Config::get_initial_species(
  const Environment &env) const {
  // Some parts of initial_species depend on the environment, so we
  // make a copy of initial_species and finish configuring it.
  std::vector<SpeciesParameters> spv = initial_species;
  for (auto &sp: spv)
    set_params_from_env(sp, env);
  return spv;
}
