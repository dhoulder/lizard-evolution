// -*- coding: utf-8 -*-

/**
   Implementation of DREaD_ds model. See ./README.md
*/
#include <memory>
#include <vector>
#include <list>
#include <map>
#include <math.h>


#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

// FIXME for debugging
#include <iostream>
#include <iomanip>

#include "dread_ds.h"

static const int max_env_dims = 2; // Max environment layers (e.g. 2
				   // for temperature, precipitation)
static const int max_genetic_dims= 3; // Max number of abstract genetic axes


struct EnvChange {
  double ramp = 0.0; // linear environment change per time step
  float sine_period = 4.0f; // wavelength of sinusoidal environment
			    // change in time steps
  float sine_offset = 0.0f; // shift sinusoidal change by this
			    // fraction of the wavelength. Use 0.25
			    // for cos()
  float sine_amplitude = 0.0f; // maximum swing of sinusoidal environment change
};

typedef boost::mt19937 rng_eng_t;
typedef boost::uniform_real<float> gene_flow_distr_t;
typedef boost::variate_generator<rng_eng_t&, gene_flow_distr_t> gene_flow_vg_t;

static rng_eng_t rng; // FIXME seed?


struct Config {
  // Parameters for a simulation run.
  int debug = 0;
  int env_dims = 1; // must be <= max_env_dims
  int genetic_dims = max_genetic_dims; // <= max_genetic_dims
  EnvChange env_change[max_env_dims];

  float gene_flow_threshold = 0.001f;
  float gene_flow_zero_distance = 5.0f;
};


struct EnvCell {
  // has to be a class|struct to keep boost::multi_array allocator happy
  float v[max_env_dims];
};


typedef boost::multi_array<EnvCell, 2> EnvMatrix;
typedef EnvMatrix::index EnvIndex;

typedef int Timestep;


struct Characteristics {
  // Describes a species

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


  Niche niche[max_env_dims];   // Niches derived from all demes of this species.
  Genetics genetics[max_genetic_dims];
  Range range;
};



class Deme {
  /**
     Describes a genetically homogeneous population in a cell.
  */
public:

  struct Genetics {
    float niche_centre[max_env_dims];
    float niche_tolerance[max_env_dims];
    float genetic_position[max_genetic_dims]; // genetic position in n-dimensional space. See struct Genetics
  };

  Genetics genetics;
  float amount; // population per cell
  bool is_primary; // indicates incumbency in a cell during dispersal

  Deme(): amount(0), is_primary(false) {}

  Deme(const Deme &from, float new_amount, bool new_primary):
    Deme(from) {
    // FIXME only need to do env_dims, not max_env_dims. Maybe use float xxxx[nnnn] = {0.0f} just to be sure? or pass in nnnn
    amount = new_amount;
    is_primary = new_primary;
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
typedef std::list <Deme> DemeList;
typedef std::map <Location, DemeList> DemeMap;

struct DispersalWeight {
  // Describes dispersal propensity for (x,y) offset from origin cell
  // at (0,0) due to distance cost.
  int x;
  int y;
  float weight; // 0 to 1.0
};

typedef std::vector <DispersalWeight> DispersalKernel;


class DemeMixer {
  /**
   * For computing weighted average of genetic and niche positions
   */

private:
  const Config &c;
  Deme::Genetics g = {{0.0f}, {0.0f}, {0.0f}};
  float total_abundance = 0.0f;

public:
  DemeMixer(const Config &conf, const Deme &primary):
    c(conf) {
    contribute(primary);
  }

  void contribute(const Deme &d) {
    for (int i = 0; i < c.genetic_dims; ++i)
      g.genetic_position[i] += d.genetics.genetic_position[i] * d.amount;
    for (int i = 0; i < c.env_dims; ++i) {
      g.niche_centre[i] += d.genetics.niche_centre[i] * d.amount;
      g.niche_tolerance[i] += d.genetics.niche_tolerance[i] * d.amount; // FIXME CHECK
    }
    total_abundance += d.amount;
  }

  void mix_to(Deme *d, float amount) {
    // Update d with weighted average of contributing demes
    for (int i = 0; i < c.genetic_dims; ++i)
      d->genetics.genetic_position[i] = g.genetic_position[i] / total_abundance;
    for (int i = 0; i < c.env_dims; ++i) {
      d->genetics.niche_centre[i] = g.niche_centre[i] /  total_abundance;
      d->genetics.niche_tolerance[i] = g.niche_tolerance[i] / total_abundance;  // FIXME CHECK
    }
    d->amount = amount;
    d->is_primary = true;
  }

};


class Species {
public:
  /**
     Describes a species and its phylogeny
   */
  const Config &conf;
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

  gene_flow_distr_t gene_flow_distr;
  gene_flow_vg_t gene_flow_random;
  // MT TODO mutex around RNG stuff http://www.bnikolic.co.uk/blog/cpp-boost-rand-normal.html
  // FIXME make rng stuff configurable so it can be deterministic for testing

  Species(const Config &c):
    conf(c),
    demes(new(DemeMap)),
    // Generate random values between [gene_flow_threshold,
    // 1-gene_flow_threshold] to effectively apply high and low
    // cutoffs when comparing against probability below.
    gene_flow_distr(
		    c.gene_flow_threshold,
		    1.0f - c.gene_flow_threshold),
    gene_flow_random(rng, gene_flow_distr)
  {
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
	if (!(i || j))
	  // dispersal into same-cell is a special case
	  continue;
	float dd = dispersal_distance(i, j);
	if (dd <= max_dispersal_radius) {
	  dk.push_back(DispersalWeight {
	      i, j,
		1.0f - ((1.0f - dispersal_min) * dd/max_dispersal_radius)});
	}
      }
  }

  void print_kernel() {
    for (auto &&v: dk)
      std::cout <<  v.x << ", " << v.y <<  " " << v.weight << std::endl;
  }


  Deme *choose_primary(DemeList &deme_list) {
    // Find "primary" deme at random by probability weighted by abundance
    return &deme_list.back(); // FIXME STUB
  }


  float gene_flow_probability(float distance) {
    // ~1.0 for distance == 0, approaches 0 for distance >= gene_flow_zero_distance
    const float A = 14.0f;
    const float B = 0.5f;
    return (1.0f / (1.0f + exp(((distance / conf.gene_flow_zero_distance) - B) * A)));
  }


  float genetic_distance(const Deme &d1, const Deme &d2) {
    // Euclidean distance in gene space
    float sum = 0.0f;
    for (int i = 0; i < conf.genetic_dims; ++i) {
      const float d = d2.genetics.genetic_position[i] - d1.genetics.genetic_position[i];
      sum += d*d;
    }
    return sqrt(sum);
  }


  bool gene_flow_occurs(const Deme &d1, const Deme &d2) {
    // see 3.4.1 "Does  gene  flow  occur?"
    return (gene_flow_probability(genetic_distance(d1, d2)) >
	    gene_flow_random());
  }


  void merge(std::shared_ptr<DemeMap> dm) {

    // Merge demes where gene flow occurs

    for (auto &&deme_cell: *dm) {
      auto &&deme_list = deme_cell.second;
      if (deme_list.size() < 2)
	// FIXME adjust abundance? anything else?
	continue;

      // Have at least two demes in this cell, so check for gene flow
      // and merge if required.
      const Location &loc = deme_cell.first;
      // Any incumbent primary deme will be at the front, otherwise
      // we have to choose one
      Deme *first_deme =  &deme_list.front();

      // NOTE: we assume there is at most one primary deme. In the
      // current model, all demes of a species in a cell mix down into
      // one final primary deme, so this assumption holds.
      // FIXME CHECK: max 1 primary deme after gene flow?

      Deme *primary = first_deme->is_primary?
	first_deme:
	choose_primary(deme_list);

      DemeMixer mixer(conf, *primary);
      for (auto &&deme: deme_list) {
	if (&deme == primary)
	  // don't merge with self
	  continue;

	if (gene_flow_occurs(deme, *primary))
	  mixer.contribute(deme);
	// otherwise "excluded by competition
      }

      // mix down to first deme
      float pa = primary->amount;
      deme_list.resize(1);
      mixer.mix_to(first_deme, pa); // FIXME CHECK:  amount after merge
    }
  }

};

class DreadDs::Model {
public:

  Model(EnvIndex rows, EnvIndex cols): env(boost::extents[rows][cols]) {
    // FIXME WIP
  }

  ~Model() {
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

  float niche_suitability(const EnvCell &cell, const Deme::Genetics &g) {
    float v = 1.0f;
    for (int i = 0; i < conf.env_dims; i++)
      v *= suitability(cell.v[i],
		       g.niche_centre[i],
		       g.niche_tolerance[i]);
    return v;
  }


  static inline float dispersal_abundance(float source_abundance,
					  float suitability,
					  float distance_weight) {
    // See 3.3 Dispersal, equation 1
    return source_abundance * suitability * distance_weight;
  }
};


std::shared_ptr<DemeMap> DreadDs::Model::disperse(Species &species) {
  /**
   * Iterate over all extant demes of a species and disperse to
   * neighbouring cells, weighted by environmental niche
   * suitability. This can result in several demes per cell.
   */
  auto target = std::make_shared <DemeMap>();
  // Iterate over all cells where this species occurs
  for (auto &&deme_cell: *(species.demes)) {

    const Location &loc = deme_cell.first;
    const EnvCell &&source_env = get_env(loc);

    // Iterate over all demes (sort of "sub species") in the cell
    for (auto &&deme:  deme_cell.second) {
      float abundance = niche_suitability(source_env, deme.genetics); // FIXME CHECK: or do this at end of step in merge??
      if (abundance <= 0.0)
	// FIXME CHECK extinction??? check this against spec.
	continue;

      // "disperse" into the same cell. This becomes a primary deme
      // after dispersal due to incumbency. Note that this inserts at
      // the front of the list. Any primary demes will be at the front.
      (*target)[loc].emplace_front(deme,
				   dispersal_abundance(abundance,
						       abundance, // FIXME CHECK
						       1.0), // same cell, no travel cost
				   true);
      // Disperse into the area around this cell
      for (auto &&k: species.dk) {
	Location new_loc;
	new_loc.x = loc.x + k.x;
	new_loc.y = loc.y + k.y;
	if (new_loc.x < 0 || new_loc.y < 0 ||
	    new_loc.x >= env.shape()[1] || new_loc.y >= env.shape()[0])
	  continue;
	// FIXME CHECK check target_abundance for extinction?? (<=0.0). check against spec
	(*target)[new_loc].emplace_back(deme,
					dispersal_abundance(abundance,
							    niche_suitability(get_env(new_loc),
									      deme.genetics),
							    k.weight),
					false);
      }
    }
  }
  return target;
}


DreadDs::DreadDs(int cols, int rows): // FIXME rows cols should come from env grid
  current_step(0), model(new Model(cols, rows)) {

  // TODO load config
  model->conf.debug = 4;

  // TODO load env

  // TODO initialise env

  // TODO initialise roots
  {
    // FIXME STUB init tips
    model->tips.push_back(Species(model->conf));
    if (model->conf.debug > 3) {
      model->tips.back().print_kernel();
    }
  }

}

DreadDs::~DreadDs() = default;


int DreadDs::step() {
  model->update_environment(current_step);

  for (auto && species: model->tips) {
    auto target = model->disperse(species);

    // TODO handle range contraction (extinction) here ???? see "3.3.2 Range contraction"

    // merge demes in each cell that are within genetic tolerance.
    species.merge(target);

    // TODO? can do this row-lagged in the dispersal loop, providing we
    // stay far enough behind the dispersal area
    // TODO?? do this in another thread?

    // Finished with source demes now - replace with target;
    species.demes = target;

    // TODO competition/co-occurrence. (see 3.5)

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
