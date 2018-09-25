/**
   Implementation of DREaD_ds model. See ./README.md
*/
#include <memory>
#include <vector>
#include <map>
#include <math.h>

#include <Eigen/Dense>

// FIXME for debugging
#include <iostream>
#include <iomanip>

#include "dread_ds.h"

using Eigen::MatrixXd;

typedef int Timestep;

struct Niche {
  /**
     Describes a niche on an environmental variable for a species.
  */
  // mean and sd of niche position of all demes of this species
  float position_mean;
  float position_sd;
  float breadth_mean;
  float breadth_sd;
  // max and min values of position - (breadth  /2)
  float max;
  float min;
};

struct Genetics {
  /**
     Holds the genetic position (on an abstract genetic trait) and
     variance of all the demes of a species.
   */
  float position;
  float variance;
};

struct Range {
  int cell_count; // number of demes (occupied cells).
  float population; // total population across all occupied cells
};

struct Characteristics {
  // Describes a species
  // FIXME std::vector may be overkill (constant lengths)
  std::vector <Niche> niche;   // Niches derived from all demes of this species.
  std::vector <Genetics> genetics;
  Range range;
};

struct Deme {
  /**
     Holds characteristics of a genetically homogeneous population in
     a cell.
  */
  // FIXME std::vector may be overkill (constant lengths)
  std::vector <float> niche_position;
  std::vector <float> niche_sd;
  float amount; // population per cell
  std::vector <float> genetic_position; // genetic position in n-dimensional space. See struct Genetics
};

// Cells occupied by demes of a species, indexed by flattened cell
// index of environment matrix.  Several demes can occupy a cell,
// hence std::vector
typedef std::map <int, std::vector <Deme>> DemeMap;

struct Spacies;
struct Species {
  /**
     Describes a species and its phylogeny
   */
  Species *parent;
  std::shared_ptr <Species> left_child = NULL;
  std::shared_ptr <Species> right_child = NULL;

  Timestep extinction; // Time step of extinction, or -1 if extant
  Timestep split; // Time of speciation. parent->split is species
		  // origin time. -1 if not speciated
  Characteristics initial; // At species origin (i.e. split from parent)
  Characteristics latest; // Updated at each time step. Frozen after speciation.
  DemeMap demes; // Cells occupied by this species.
};

struct DispersalWeight {
  // Describes dispersal propensity for (x,y) offset from origin cell
  // at (0,0) due to distance cost.
  int x;
  int y;
  float weight; // 0 to 1.0
};

typedef std::vector<DispersalWeight> DispersalKernel;

class Config {
public:
  // Parameters for a simulation run.
  int debug = 0;
  float max_dispersal_radius;
  float dispersal_floor; // dispersal weight lower than this is treated as 0
  double env_ramp; // linear environment change per time step
  float sine_period; // wavelength of sinusoidal environment change in time steps
  float sine_offset; // shift sinusoidal change by this fraction of the
		     // wavelength. Use 0.25 for cos()
  float sine_amplitude; // maximum swing of sinusoidal environment change

  float env_change(int time_step) {
    return (time_step * env_ramp) +
      (sine_amplitude * sin(2 * M_PI *
			    ((double)sine_offset + (double)time_step/sine_period)));
  }
};

class DreadDs::Impl {
public:

  Impl(int cols, int rows): env(MatrixXd(cols, rows)) {
    // FIXME WIP
  }

  ~Impl() {
    //FIXME WIP
  }

  void setup_dispersal();
  void disperse(Species &species);


  void update_environment(int step) {
    // This offset gets applied to every cell in env
    env_offset = conf.env_change(step);

    if (conf.debug > 3)
      std::cout << step << " env " << env_offset << std::endl;
  }

  Config conf;
  MatrixXd env;
  float env_offset = 0.0f;
  DispersalKernel dk;
  std::vector <Species> roots; // Initial species
  std::vector <Species> tips; // extant leaf species
};

static float dispersal_distance(long x, long y) {
  // Replace with some other distance metric if required.
  return sqrt(x*x + y*y);
}

void DreadDs::Impl::setup_dispersal() {
  // Calculate dispersal kernel
  // TODO only really have to store one quadrant (perhaps only one octant)
  int r = (int) (ceil(conf.max_dispersal_radius) + 0.5);
  for (int i = -r; i <= r; ++i)
    for (int j = -r; j <= r; ++j) {
      float v = (1.0f - ((dispersal_distance(i, j)) / conf.max_dispersal_radius));
      if (v >= conf.dispersal_floor) {
	dk.push_back(DispersalWeight {i, j, v});
      }
    }
  if (conf.debug > 3)
    for (std::vector<DispersalWeight>::iterator it = dk.begin() ; it != dk.end(); ++it)
      std::cout <<  it->x << ", " << it->y <<  " " << it->weight << std::endl;
};

void DreadDs::Impl::disperse(Species &species) {
  // FIXME WIP

  // FIXME setup target DemeMap


  for (DemeMap::iterator d_it = species.demes.begin();
       d_it != species.demes.end();
       ++d_it) {

    // FIXME disperse to target


  }

}


DreadDs::DreadDs(int cols, int rows): current_step(0), impl(new Impl(cols, rows)) {

  // TODO load config
  impl->conf = { // FIXME STUB
    4, // debug level
    2, // dispersal radius
    0.2, // dispersal floor
    0, // env_ramp
    4, // sine_period
    0, // sine_offset,
    0 // sine_amplitude
  };

  impl->setup_dispersal();

  // TODO initialise env

  // TODO initialise roots, tips
}

DreadDs::~DreadDs() = default;


int DreadDs::step() {
  impl->update_environment(current_step);

  for (std::vector<Species>::iterator s_it = impl->tips.begin();
       s_it != impl->tips.end();
       ++s_it) {
    impl->disperse(*s_it);


    // TODO collapse demes in cell within genetic tolerance.
    // TODO? can do this row-lagged in the dispersal loop, providing we
    // stay far enough behind the dispersal area
    // TODO?? do this in another thread?

    // TODO niche evolution

    // TODO check for extinction

    // TODO genetic drift


    // TODO speciate.
    // Make sure iterator doesn't see this.

  }

  ++current_step;
}


int DreadDs::run(int n_steps) {
  for (int final = current_step + n_steps;
       current_step < final; ) {
    // FIXME STUB. needs some way to stop when done
    step();
  }
  return current_step;
}
