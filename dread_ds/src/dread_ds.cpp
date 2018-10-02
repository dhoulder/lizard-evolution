// -*- coding: utf-8 -*-

/**
   Implementation of DREaD_ds model. See ./README.md
*/
#include <memory>
#include <vector>
#include <map>
#include <math.h>


#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

// FIXME for debugging
#include <iostream>
#include <iomanip>

#include "dread_ds.h"

static const int max_env_dims = 2; // Max environment layers (e.g. 2
				   // for temperature, precipitation)
static const int max_genetic_dims= 3; // Max number of abstract genetic axes

struct EnvCell {
  // has to be a class|struct to keep boost::multi_array allocator happy
  float v[max_env_dims];
};


typedef boost::multi_array<EnvCell, 2> EnvMatrix;
typedef EnvMatrix::index EnvIndex;

typedef int Timestep;

struct Niche {
  /**
     Describes a niche on an environmental variable for a species.
  */
  // mean and sd of niche position of all demes of this species
  float position_mean = 0.0f;
  float position_sd = 0.0f;
  float breadth_mean = 0.0f;
  float breadth_sd = 0.0f;
  // max and min values of position - (breadth  /2)
  float max = 0.0f;
  float min = 0.0f;
};

struct Genetics {
  /**
     Holds the genetic position (on an abstract genetic trait) and
     variance of all the demes of a species.
   */
  float position = 0.0f;
  float variance = 0.0f;
};

struct Range {
  int cell_count; // number of demes (occupied cells).
  float population; // total population across all occupied cells
};

struct Characteristics {
  // Describes a species

  Niche niche[max_env_dims];   // Niches derived from all demes of this species.
  Genetics genetics[max_genetic_dims];
  Range range;
};

class Deme {
  /**
     Describes a genetically homogeneous population in a cell.
  */
public:
  float niche_centre[max_env_dims];
  float niche_tolerance[max_env_dims];
  float amount; // population per cell
  float genetic_position[max_genetic_dims]; // genetic position in n-dimensional space. See struct Genetics

  Deme(const Deme &from, float new_amount): Deme(from) {
    // FIXME only need to do env_dims, not max_env_dims. Maybe use float xxxx[nnnn] = {0.0f} just to be sure? or pass in nnnn
    amount = new_amount;
  }
};

struct Location {
  int x;
  int y;

  // Used as key in a map, so needs an ordering
  friend bool operator< (const Location &a, const Location &b) {
    return (a.x < b.x) || (a.x == b.x && a.y < b.y);
  }
};

// Cells occupied by demes of a species, Several demes can occupy a
// cell, hence std::vector
typedef std::map <Location, std::vector <Deme>> DemeMap;

struct DispersalWeight {
  // Describes dispersal propensity for (x,y) offset from origin cell
  // at (0,0) due to distance cost.
  int x;
  int y;
  float weight; // 0 to 1.0
};

typedef std::vector <DispersalWeight> DispersalKernel;

class Spacies;
class Species {
public:
  /**
     Describes a species and its phylogeny
   */
  Species *parent = NULL;
  std::shared_ptr <Species> left_child = NULL;
  std::shared_ptr <Species> right_child = NULL;

  Timestep extinction = -1; // Time step of extinction, or -1 if extant
  Timestep split = -1; // Time of speciation. parent->split is species
		  // origin time. -1 if not speciated
  Characteristics initial; // At species origin (i.e. split from parent)
  Characteristics latest; // Updated at each time step. Frozen after speciation.
  std::shared_ptr <DemeMap> demes; // Cells occupied by this species.

  float max_dispersal_radius = 1.0f;
  float dispersal_min = 0.2;
  DispersalKernel dk;

  Species(): demes(new(DemeMap)) {
    setup_dispersal();
  }

  static float dispersal_distance(int x, int y) {
    // Replace with some other distance metric if required.
    return sqrt(x*x + y*y);
  }

  void setup_dispersal() {
    // Calculate dispersal kernel
    // TODO only really have to store one quadrant (perhaps only one octant)
    int r = (int) (ceil(max_dispersal_radius) + 0.5);
    for (int i = -r; i <= r; ++i)
      for (int j = -r; j <= r; ++j) {
	float dd = dispersal_distance(i, j);
	if (dd <= max_dispersal_radius) {
	  dk.push_back(DispersalWeight {
	      i, j,
		1.0f - ((1.0f - dispersal_min) * dd/max_dispersal_radius)});
	}
      }
  }

  void print_kernel() {
    for (std::vector<DispersalWeight>::iterator it = dk.begin() ; it != dk.end(); ++it)
      std::cout <<  it->x << ", " << it->y <<  " " << it->weight << std::endl;
  }

};

struct EnvChange {
  double ramp = 0.0; // linear environment change per time step
  float sine_period = 4.0f; // wavelength of sinusoidal environment change in time steps
  float sine_offset = 0.0f; // shift sinusoidal change by this fraction of the
		     // wavelength. Use 0.25 for cos()
  float sine_amplitude = 0.0f; // maximum swing of sinusoidal environment change
};

struct Config {
  // Parameters for a simulation run.
  int debug = 0;
  int env_dims = 1; // must be <= max_env_dims
  int genetic_dims = max_genetic_dims; // <= max_genetic_dims
  EnvChange env_change[max_env_dims];
};

class DreadDs::Impl {
public:

  Impl(EnvIndex rows, EnvIndex cols): env(boost::extents[rows][cols]) {
    // FIXME WIP
  }

  ~Impl() {
    //FIXME WIP
  }

  std::shared_ptr<DemeMap> disperse(Species &species);

  Config conf;
  EnvMatrix env;
  float env_delta[max_env_dims] = {0.0f};
  std::vector <Species> roots; // Initial species
  std::vector <Species> tips; // extant leaf species

  void update_environment(int time_step) {
    // Set env_delta for the current time step. env_delta gets applied
    // to every cell in env
    float *ed  = env_delta;
    EnvChange *c = conf.env_change;
    for (int i = 0; i < conf.env_dims;
	 ++i, ++c, ++ed) {
      *ed = (time_step * c->ramp) +
	(c->sine_amplitude * sin(2 * M_PI *
				 ((double)c->sine_offset +
				  (double)time_step/c->sine_period)));
      if (conf.debug > 3)
	std::cout << time_step << " env " << *ed << std::endl;
    }
  }

  EnvCell get_env(const Location &loc) {
    EnvCell ec;
    const EnvCell &base_env =  env[loc.y][loc.x];
    for (int i =0; i < conf.env_dims; i++)
      ec.v[i] =  base_env.v[i] + env_delta[i];
    return ec;
  }

  static float suitability(float env_value, float niche_centre, float niche_tolerance) {
    // The suitability function is cos() from -pi to pi, scaled to the
    // range 0…1.0. This gives 1.0 at niche_centre, and 0 and dy/dx ==
    // 0 at niche_centre±niche_tolerance
    return
      (fabs(niche_centre - env_value) > niche_tolerance)?
      0.0f:
      0.5f + 0.5f * cos(M_PI * (env_value - niche_centre) / niche_tolerance);
  }

  float niche_suitability(const EnvCell &cell, const Deme &d) {
    float v = 1.0f;
    for (int i = 0; i < conf.env_dims; i++)
      v *= suitability(cell.v[i], d.niche_centre[i], d.niche_tolerance[i]);
    return v;
  }

};


std::shared_ptr<DemeMap> DreadDs::Impl::disperse(Species &species) {
  auto target = std::make_shared <DemeMap>();
  // Iterate over all cells where this species occurs
  for (DemeMap::iterator dm_it = species.demes->begin();
       dm_it != species.demes->end();
       ++dm_it) {
    const Location &loc = dm_it->first;
    std::vector <Deme> &dv =  dm_it->second;
    EnvCell &&source_env = get_env(loc);

    // Iterate over all demes (sort of "sub species") in the cell
    for (auto deme_itr =  dv.begin();
	 deme_itr != dv.end();
	 ++deme_itr) {
      float abundance = niche_suitability(source_env, *deme_itr);
      // Disperse into the area around this cell
      for (auto k = species.dk.begin();
	   k != species.dk.end();
	   ++k) {
	Location new_loc;
	new_loc.x = loc.x + k->x;
	new_loc.y = loc.y + k->y;
	if (new_loc.x < 0 || new_loc.y < 0 ||
	    new_loc.x >= env.shape()[1] || new_loc.y >= env.shape()[0])
	  continue;

	EnvCell &&target_env = get_env(new_loc);
	float target_abundance = abundance * niche_suitability(target_env, *deme_itr) * k->weight;
	std::vector <Deme> &target_dv = (*target)[new_loc]; // new or existing vector at location
	// Disperse to new deme. Everything except amount in the source deme gets
	// copied into the target deme.
	target_dv.emplace_back(*deme_itr, target_abundance);
      }
    }
  }
  return target;
}


DreadDs::DreadDs(int cols, int rows): // FIXME rows cols should come from env grid
  current_step(0), impl(new Impl(cols, rows)) {

  // TODO load config
  impl->conf.debug = 4;

  // TODO load env

  // TODO initialise env

  // TODO initialise roots, tips
  {
    // FIXME STUB
    impl->roots.push_back(Species {});
    if (impl->conf.debug > 3) {
      impl->roots.back().print_kernel();
    }
  }


}

DreadDs::~DreadDs() = default;


int DreadDs::step() {
  impl->update_environment(current_step);

  for (std::vector<Species>::iterator s_it = impl->tips.begin();
       s_it != impl->tips.end();
       ++s_it) {
    auto target = impl->disperse(*s_it);


    // TODO collapse demes in cell within genetic tolerance.
    // TODO? can do this row-lagged in the dispersal loop, providing we
    // stay far enough behind the dispersal area
    // TODO?? do this in another thread?

    // Finished with source demes now - replace with target;
    s_it->demes = target;

    // TODO niche evolution. Remove demes with 0 abundance

    // TODO check for extinction. Remove from DemeMap

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
