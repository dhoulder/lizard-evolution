// -*- coding: utf-8 -*-

/**
 *  Implementation of DREaD_ds model. See ./README.md
 */
#include <memory>
#include <math.h>

// FIXME for debugging
#include <iostream>
#include <iomanip>

#include "model-config.h"
#include "environment.h"
#include "model.h"

using namespace DreadDs;

namespace  DreadDs  {

  static rng_eng_t rng; // FIXME seed?

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
	g.niche_tolerance[i] += d.genetics.niche_tolerance[i] * d.amount;
      }
      total_abundance += d.amount;
    }

    void mix_to(Deme *d) {
      // Update d with weighted average of contributing demes
      // FIXME handle / 0.0
      for (int i = 0; i < c.genetic_dims; ++i)
	d->genetics.genetic_position[i] = g.genetic_position[i] / total_abundance;
      for (int i = 0; i < c.env_dims; ++i) {
	d->genetics.niche_centre[i] = g.niche_centre[i] /  total_abundance;
	d->genetics.niche_tolerance[i] = g.niche_tolerance[i] / total_abundance;
      }
      d->is_primary = true;
    }
  };
}

Model::Model(const char *config_path,
	     const char *output_path_arg):
  Model(Config(config_path),
	output_path_arg) {}


Model::Model(const Config &conf_,
	     const char *output_path_arg):
  conf(conf_),
  env(Environment(conf_.env_params)),
  // Generate random values between [gene_flow_clip,
  // 1-gene_flow_clip] to effectively apply high and low
  // cutoffs when comparing against probability below.
  gene_flow_distr(
		  conf_.gene_flow_clip,
		  1.0f - conf_.gene_flow_clip),
  gene_flow_random(rng, gene_flow_distr),
  deme_choice_distr(0.0f, 1.0f),
  gene_drift_distr(0.0f,
		   conf_.gene_drift_sd),
  gene_drift_random(rng, gene_drift_distr)
{
  for (const auto &sp: conf.initial_species) {
    tips.push_back(std::shared_ptr <Species>(new Species(conf, sp, env)));
  }
  roots = tips;
  output_path = output_path_arg;
}

// MT TODO mutex around RNG stuff http://www.bnikolic.co.uk/blog/cpp-boost-rand-normal.html
// FIXME make rng stuff configurable so it can be deterministic for testing

void Model::update_environment(int time_step) {
  // Set env_delta for the current time step. env_delta gets applied
  // to every cell in env
  float *ed  = env_delta;
  auto p = conf.env_params.begin();
  for (int i = 0; i < conf.env_dims;
       ++i, ++p, ++ed) {
    *ed = (time_step * p->ramp) +
      (p->sine_amplitude * sin(2 * M_PI *
			       ((double)p->sine_offset +
				(double)time_step/p->sine_period)));
    if (conf.debug > 3)
      std::cout << "env " << i << " step " << time_step << " delta " << *ed << std::endl;
  }
}

inline EnvCell Model::get_env(const Location &loc) {
  EnvCell ec;
  const EnvCell &base_env =  (env.values)[loc.y][loc.x];
  for (int i =0; i < conf.env_dims; i++)
    ec.v[i] =  base_env.v[i] + env_delta[i];
  return ec;
}

static float suitability(float env_value, float niche_centre, float niche_tolerance) {
  // The suitability function is cos() from -π to π, scaled to the
  // range 0…1.0. This gives 1.0 at niche_centre, and 0 and dy/dx ==
  // 0 at niche_centre±niche_tolerance
    return
      (fabs(niche_centre - env_value) > niche_tolerance)?
      0.0f:
      0.5f + 0.5f * cos(M_PI * (env_value - niche_centre) / niche_tolerance);
}

float Model::niche_suitability(const EnvCell &cell, const Deme::Genetics &g) {
  // compute geometric mean of all niche suitabilities
  float v = 1.0f;
  for (int i = 0; i < conf.env_dims; i++)
    v *= suitability(cell.v[i],
		     g.niche_centre[i],
		     g.niche_tolerance[i]);
  return pow(v,
	     1.0f / (float)conf.env_dims);
}


static inline float dispersal_abundance(float source_abundance,
					float suitability,
					float distance_weight) {
  // See 3.3 Dispersal, equation 1
  return source_abundance * suitability * distance_weight;
}


void Model::evolve_towards_niche(Deme &deme, const EnvCell &ec) {
  // apply niche selection pressure
  float *nc =  deme.genetics.niche_centre;
  for (int i=0; i < conf.env_dims; nc++, i++)
    *nc += (ec.v[i] - *nc) * conf.niche_evolution_rate;
}

void Model::do_genetc_drift(Deme &deme) {
  float * gp = deme.genetics.genetic_position;
  for (int i=0; i < conf.genetic_dims; gp++, i++)
    *gp += gene_drift_random();
}


std::shared_ptr<DemeMap> Model::disperse(Species &species) {
  /**
   * Iterate over all extant demes of a species and disperse to
   * neighbouring cells, weighted by environmental niche
   * suitability. This can result in several demes per cell.
   */
  auto target = std::make_shared <DemeMap>(); // FIXME not C++11
  auto &&env_shape = env.values.shape();

  // Iterate over all cells where this species occurs
  for (auto &&deme_cell: *(species.demes)) {
    const Location &loc = deme_cell.first;
    const EnvCell &&source_env = get_env(loc);

    // Iterate over all demes (sort of "sub species") in the cell
    for (auto &&deme:  deme_cell.second) {
      evolve_towards_niche(deme, source_env);
      // TODO check for extinction. Remove from DemeMap
      // CHECK update abundance?
      do_genetc_drift(deme);

      // "disperse" into the same cell. This becomes a primary deme
      // after dispersal due to incumbency. Note that this inserts at
      // the front of the list. Any primary demes will be at the front.
      (*target)[loc].emplace_front(deme,
				   dispersal_abundance(deme.amount,
						       niche_suitability(source_env,
									 deme.genetics),
						       1.0), // same cell, no travel cost
				   true);
      // Disperse into the area around this cell
      for (auto &&k: species.dk) {
	Location new_loc;
	new_loc.x = loc.x + k.x;
	new_loc.y = loc.y + k.y;
	if (new_loc.x < 0 || new_loc.y < 0 ||
	    new_loc.x >= env_shape[1] || new_loc.y >= env_shape[0])
	  continue;
	// FIXME CHECK check target_abundance for extinction?? (<=0.0). check against spec
	(*target)[new_loc].emplace_back(deme,
					dispersal_abundance(deme.amount,
							    niche_suitability(get_env(new_loc),
									      deme.genetics),
							    k.weight),
					false);
      }
    }
  }
  return target;
}


Deme *Model::choose_primary(DemeList &deme_list) {
  // Find "primary" deme at random by abundance-weighted probability
  // Do not call on empty list
  float sum = 0.0f;;
  for (auto &&d: deme_list)
    sum += d.amount;
  float rv = deme_choice_distr(rng) * sum;
  for (auto &&d: deme_list) {
    if (rv <= d.amount )
      return &d;
    rv -= d.amount;
  }
  // fp precision effects may leave us here
  return &deme_list.back();
}


float Model::gene_flow_probability(float distance) {
  // ~1.0 for distance == 0, approaches 0 for distance >= gene_flow_zero_distance
  const float A = 14.0f;
  const float B = 0.5f;
  return (1.0f / (1.0f + exp(((distance / conf.gene_flow_zero_distance) - B) * A)));
}


float Model::genetic_distance(const Deme &d1, const Deme &d2) {
  // Euclidean distance in gene space
  float sum = 0.0f;
  for (int i = 0; i < conf.genetic_dims; ++i) {
    const float d = d2.genetics.genetic_position[i] - d1.genetics.genetic_position[i];
    sum += d*d;
  }
  return sqrt(sum);
}


bool Model::gene_flow_occurs(const Deme &d1, const Deme &d2) {
  // see 3.4.1 "Does  gene  flow  occur?"
  return (gene_flow_probability(genetic_distance(d1, d2)) >
	  gene_flow_random());
}


void Model::merge(DemeMap &dm) {
  // Merge demes where gene flow occurs

  for (auto &&deme_cell: dm) {
    auto &&deme_list = deme_cell.second;
    if (deme_list.size() < 1)
      continue;

    const Location &loc = deme_cell.first;
    Deme *first_deme =  &deme_list.front();

    if (deme_list.size() >1) {
      // Have at least two demes in this cell, so check for gene flow
      // and merge if required.

      // Any incumbent primary deme will be at the front, otherwise
      // we have to choose one

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
      deme_list.resize(1);
      mixer.mix_to(first_deme);
    }
    // update abundance according to current environment
    first_deme->amount = niche_suitability(get_env(loc),
					   first_deme->genetics);

    // FIXME CHECK extinction??? check this against spec.
  }
}


int Model::do_step() {
  /**
   * Execute one time step of the model
   */

  update_environment(step);
  for (auto && species: tips) {
    auto target = disperse(*species);
    // TODO handle range contraction (extinction) here ???? see "3.3.2 Range contraction"

    // merge demes in each cell that are within genetic tolerance.
    merge(*target);

    // TODO? can do this row-lagged in the dispersal loop, providing we
    // stay far enough behind the dispersal area
    // TODO?? do this in another thread?

    // Finished with source demes now - replace with target;
    species->demes = target;

    // TODO competition/co-occurrence. (see 3.5)

    // TODO speciate.
    // Make sure iterator doesn't see this.
  }

  ++step;
  return step;
}
