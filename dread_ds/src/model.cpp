// -*- coding: utf-8 -*-

/**
 * Implementation of DREaD_ds model. See ./README.md
 */

// For C open(), fdopen() etc.
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>

#include <memory>
#include <math.h>
#include <iostream>
#include <string>

#include "model-config.h"
#include "environment.h"
#include "deme.h"
#include "species.h"
#include "model.h"
#include "exceptions.h"

using namespace DreadDs;

namespace DreadDs {

  static rng_eng_t rng;

  class DemeMixer {
    /**
     * For computing weighted average of genetic and niche positions
     */
  private:
    const Config &c;
    Deme::Genetics g;
    float total_abundance;

  public:
    DemeMixer(const Config &conf):
      c {conf},
      g {},
      total_abundance {0.0f}
    {}

    void contribute(const Deme &d) {
      for (int i = 0; i < c.genetic_dims; ++i)
	g.genetic_position[i] += d.genetics.genetic_position[i] * d.amount;
      for (int i = 0; i < c.env_dims; ++i) {
	g.niche_centre[i] += d.genetics.niche_centre[i] * d.amount;
	g.niche_tolerance[i] += d.genetics.niche_tolerance[i] * d.amount;
      }
      total_abundance += d.amount;
    }

    bool mix_to(Deme *d) {
      // If possible, update d with weighted average of contributing demes
      if (total_abundance <= 0.0f)
	return false;
      for (int i = 0; i < c.genetic_dims; ++i)
	d->genetics.genetic_position[i] = g.genetic_position[i] / total_abundance;
      for (int i = 0; i < c.env_dims; ++i) {
	d->genetics.niche_centre[i] = g.niche_centre[i] / total_abundance;
	d->genetics.niche_tolerance[i] = g.niche_tolerance[i] / total_abundance;
      }
      d->is_primary = true;
      return true;
    }
  };
}

Model::Model(const Config &c):
  conf(c),
  step(0),
  env(Environment(conf)), // must use conf member, not arg

  // Generate random values between [gene_flow_clip,
  // 1-gene_flow_clip] to effectively apply high and low
  // cutoffs when comparing against probability below.
  gene_flow_distr(
		  conf.gene_flow_clip,
		  1.0f - conf.gene_flow_clip),
  gene_flow_random(rng, gene_flow_distr),
  deme_choice_distr(0.0f, 1.0f),
  gene_drift_distr(0.0f,
		   conf.gene_drift_sd),
  gene_drift_random(rng, gene_drift_distr)
{
  env.update(step);
  for (const auto &sp: conf.initial_species) {
    tips.push_back(std::shared_ptr <Species>(new Species(conf, sp, env)));
  }
  roots = tips;
}

// MT TODO mutex around RNG stuff http://www.bnikolic.co.uk/blog/cpp-boost-rand-normal.html
// FIXME make rng stuff configurable so it can be deterministic for testing

static inline float dispersal_abundance(float source_abundance,
					float suitability,
					float distance_weight) {
  // See 3.3 Dispersal, equation 1
  return source_abundance * suitability * distance_weight;
}


void Model::evolve_towards_niche(Deme &deme, const EnvCell &ec) {
  // apply niche selection pressure
  float *nc = deme.genetics.niche_centre;
  for (int i=0; i < conf.env_dims; ++nc, ++i)
    *nc += (ec.v[i] - *nc) * conf.niche_evolution_rate;
}

void Model::do_genetc_drift(Deme &deme) {
  float * gp = deme.genetics.genetic_position;
  for (int i=0; i < conf.genetic_dims; gp++, i++)
    *gp += gene_drift_random();
}


std::shared_ptr<DemeMap> Model::evolve_and_disperse(Species &species) {
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
    const EnvCell &&source_env = env.get(loc);

    // Iterate over all demes (sort of "sub species") in the cell
    for (auto &&deme: deme_cell.second) {
      evolve_towards_niche(deme, source_env);
      // FIXME MAYBE: update abundance?
      do_genetc_drift(deme);

      // "disperse" into the same cell. This becomes a primary deme
      // after dispersal due to incumbency. Note that this inserts at
      // the front of the list. Any primary demes will be at the front.
      (*target)[loc].emplace_front(
	      deme,
	      dispersal_abundance(
		      deme.amount,
		      deme.niche_suitability(conf, source_env),
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
	(*target)[new_loc].emplace_back(
		deme,
		dispersal_abundance(
			deme.amount,
			deme.niche_suitability(conf, env.get(new_loc)),
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
  // see 3.4.1 "Does gene flow occur?"
  return (gene_flow_probability(genetic_distance(d1, d2)) >
	  gene_flow_random());
}


void Model::merge(DemeMap &dm) {
  // Merge demes where gene flow occurs

  for (auto cell_itr = dm.begin(); cell_itr != dm.end(); ) {

    auto &&deme_list = cell_itr->second;
    if (deme_list.size() < 1) {
      // no demes, so get rid of the cell
      cell_itr = dm.erase(cell_itr);
      continue;
    }

    const Location &loc = cell_itr->first;
    Deme *first_deme = &deme_list.front();

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

      DemeMixer mixer(conf);
      for (auto &&deme: deme_list)
	if ((&deme == primary) || gene_flow_occurs(deme, *primary))
	  mixer.contribute(deme);
	// otherwise "excluded by competition

      // mix down to first deme
      if (!mixer.mix_to(first_deme)) {
	// extinct in this cell, so drop it
	  cell_itr = dm.erase(cell_itr);
	  continue;
      }
      deme_list.resize(1);
    }
    // update abundance according to current environment
    first_deme->amount = first_deme->niche_suitability(
	    conf, env.get(loc));
    if (first_deme->amount <= 0.0f)
      // another extinction case
      cell_itr = dm.erase(cell_itr);
    else
      ++cell_itr;
  }
}


void Model::save() {
  // Using plain old C file IO because C++ doesn't provide an
  // equivalent to O_CREAT|O_EXCL.  Could also use
  // https://www.boost.org/doc/libs/1_67_0/libs/iostreams/doc/classes/file_descriptor.html#file_descriptor_sink

  std::string output_filename = conf.output_dir + "/" +
    conf.output_file_prefix + std::to_string(step) + ".csv";
  int fd = open(output_filename.c_str(),
		// atomically create and open, fail if exists
		O_WRONLY | O_CREAT | O_EXCL,
		0664);
  FILE *of = NULL;
  if (fd > -1)
    of = fdopen(fd, "w");
  if (NULL == of) {
    char msg_buf[200] = "";
#if (_POSIX_C_SOURCE >= 200112L) && !  _GNU_SOURCE
    strerror_r(errno, msg_buf, sizeof(msg_buf));
    char *mb = msg_buf;
#else
    char *mb = strerror_r(errno, msg_buf, sizeof(msg_buf));
#endif
    throw ApplicationException("Could not create " +
			       output_filename + ": " + mb);
  }

  fprintf(of, "species, row, column, amount");
  for (int i=0; i < conf.env_dims; ++i)
    fprintf(of, ", env_%d, niche_centre_%d, niche_breadth_%d", i, i, i);
  for (int i=0; i < conf.genetic_dims; i++)
    fprintf(of, ", genetic_position_%d", i);
  fprintf(of, "\n");

  for (auto &&species: tips) {
    for (auto &&deme_cell: *(species->demes)) {
      const Location &loc = deme_cell.first;
      const EnvCell &&ec = env.get(loc);
      for (auto &&deme: deme_cell.second) {
	auto &&g = deme.genetics;
	fprintf(of,
		"%d, %d, %d, %f",
		species->id, loc.y, loc.x, deme.amount);
	for (int i=0; i < conf.env_dims; ++i)
	  fprintf(of,
		  ", %f, %f, %f",
		  ec.v[i],
		  g.niche_centre[i], g.niche_tolerance[i] *2.0f);
	for (int i=0; i < conf.genetic_dims; i++)
	  fprintf(of,
		  ", %f",
		  g.genetic_position[i]);
	fprintf(of, "\n");
      }
    }
  }

  fclose(of);

}


int Model::do_step() {
  /**
   * Execute one time step of the model
   */
  ++step;
  env.update(step);
  int n_occupied = 0;
  for (auto && species: tips) {
    auto target = evolve_and_disperse(*species);
    // TODO handle range contraction (extinction) here ???? see "3.3.2 Range contraction"

    // merge demes in each cell that are within genetic tolerance.
    merge(*target);

    // TODO? can do this row-lagged in the dispersal loop, providing we
    // stay far enough behind the dispersal area
    // TODO?? do this in another thread?

    // Finished with source demes now - replace with target;
    species->demes = target;

    // TODO competition/co-occurrence. (see 3.5)

    // TODO speciate. (make sure auto && species iterator doesn't visit new species)

    species->update_stats(species->latest_stats, step);
    n_occupied += species->latest_stats.cell_count;

  }

  save();

  if (conf.verbosity > 1)
    std::cout << "End of step " << step << std::endl << std::endl;

  return n_occupied;
}
