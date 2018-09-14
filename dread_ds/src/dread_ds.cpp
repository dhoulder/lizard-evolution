/**
   Implementation of DREaD_ds model. See ./README.md
*/

#include <vector>
#include <map>

#include <Eigen/Dense>

#include "dread_ds.h"

using Eigen::MatrixXd;


struct Niche {
  /**
     Describes a niche on an environmental variable for a species.
  */
  // mean and sd of niche position of all demes of this species
  float position_mean;
  float position_sd;
  float breadth_mean;
  float breadth_sd;
  // max and min  values  of  position - (breadth  /2)
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


typedef int Timestep;

struct Range {
  int cell_count; // number of demes (occupied cells).
  float population; // total population across all occupied cells
};


struct Spacies;

struct Species {
  /**
     Describes species and their phylogeny
   */
  Species *parent;
  Timestep extinction; // Time step of extinction, or -1 if extant
  Timestep split; // Time of speciation. parent->split is species
		  // origin time. -1 if not speciated
  Range initial; // At species origin (i.e. split from parent)
  Range latest; // Updated at each time step. Frozen after speciation.
		// Extinction implies empty range
  // Niches derived from all demes of this species.
  std::vector <Niche> initial_niche;
  std::vector <Niche> latest_niche; // Updated at each time step. frozen after speciation
  std::vector <Genetics> initial_genetics;
  std::vector <Genetics> latest_genetics;
};

struct Deme {
  /**
     Holds location and characteristics of a population.
  */
  Species *species;
  std::vector <float> niche_position;
  std::vector <float> niche_sd;
  float amount; // population per cell
  float genetic_position[]; // genetic position in n-dimensional space. See struct Genetics
};


class DreadDs::Impl {
public:

  Impl(int cols, int rows): env(MatrixXd::Random(cols, rows)) {
    // FIXME WIP
  }

  ~Impl() {
    //FIXME WIP
  }


  MatrixXd env;
  Species *root;
  std::vector <Species> tips; // extant leaf species
  std::multimap <int, Deme> demes; // indexed by flattened cell index of env
};


DreadDs::DreadDs(int cols, int rows): impl(new Impl(cols, rows)) {
  // FIXME WIP
}


DreadDs::~DreadDs() = default;
