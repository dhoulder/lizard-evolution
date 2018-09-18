/**
   Implementation of DREaD_ds model. See ./README.md
*/

#include <memory>
#include <vector>
#include <map>

#include <Eigen/Dense>

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

struct Range {
  int cell_count; // number of demes (occupied cells).
  float population; // total population across all occupied cells
};

struct Characteristics {
  // FIXME std::vector may be overkill (constant lengths)
  std::vector <Niche> niche;   // Niches derived from all demes of this species.
  std::vector <Genetics> genetics;
  Range range;
};

struct Deme {
  /**
     Holds characteristics of a population in a cell.
  */
  // FIXME std::vector may be overkill (constant lengths)
  std::vector <float> niche_position;
  std::vector <float> niche_sd;
  float amount; // population per cell
  float genetic_position[]; // genetic position in n-dimensional space. See struct Genetics
};



struct Spacies;

struct Species {
  /**
     Describes species and their phylogeny
   */
  Species *parent;
  std::shared_ptr <Species> left_child = NULL;
  std::shared_ptr <Species> right_child = NULL;

  Timestep extinction; // Time step of extinction, or -1 if extant
  Timestep split; // Time of speciation. parent->split is species
		  // origin time. -1 if not speciated
  Characteristics initial; // At species origin (i.e. split from parent)
  Characteristics latest; // Updated at each time step. Frozen after speciation.
  std::map <int, Deme> demes; // Cells occupied by this
			      // species. indexed by flattened cell
			      // index of env
};

class DreadDs::Impl {
public:

  Impl(int cols, int rows): env(MatrixXd(cols, rows)) {
    // FIXME WIP
  }

  ~Impl() {
    //FIXME WIP
  }

  MatrixXd env;
  std::vector <Species> roots; // Initial species
  std::vector <Species> tips; // extant leaf species
};

DreadDs::DreadDs(int cols, int rows): step(0), impl(new Impl(cols, rows)) {
  // FIXME WIP
}

DreadDs::~DreadDs() = default;

int DreadDs::run(int n_steps) {
  for (int final = step + n_steps; step < final; step++) {
    // FIXME STUB
  }
  return step;
}
