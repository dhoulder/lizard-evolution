// -*- coding: utf-8 -*-

/**
 * Implementation of DREaD_ds model. See ./README.md
 */

#include <stdio.h>

#include <ctime>
#include <memory>
#include <math.h>
#include <iostream>
#include <string>

#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "deme.h"
#include "species.h"
#include "model.h"
#include "exceptions.h"
#include "output-file.h"

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
  gene_drift_distr(0.0f, gene_drift_sd),
  gene_drift_random(rng, gene_drift_distr)
{
  for (auto &&sp: conf.get_initial_species(env)) {
    tips.push_back(std::make_shared <Species>(conf, sp, env));
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


void Model::evolve_towards_niche(Deme &deme, const EnvCell ec) {
  // apply niche selection pressure
  float *nc = deme.genetics.niche_centre;
  for (int i=0; i < conf.env_dims; ++nc, ++i)
    *nc += (ec[i] - *nc) * conf.niche_evolution_rate;
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
  EnvCell source_env, target_env;
  bool no_source, no_target;
  auto target = std::make_shared <DemeMap>();
  auto &&env_shape = env.values.shape();

  // Iterate over all cells where this species occurs
  for (auto &&deme_cell: *(species.demes)) {
    const Location &loc = deme_cell.first;
    env.get(loc, source_env, &no_source);
    if (no_source)
      continue;

    // Iterate over all demes in the cell
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
      Location new_loc;
      for (auto &&k: species.dk) {
	new_loc.x = loc.x + k.x;
	new_loc.y = loc.y + k.y;
	if (new_loc.x < 0 || new_loc.y < 0 ||
	    new_loc.x >= env_shape[1] || new_loc.y >= env_shape[0])
	  continue;
	env.get(new_loc, target_env, &no_target);
	if (no_target)
	  continue;
	// FIXME CHECK check target_abundance for extinction?? (<=0.0). check against spec
	(*target)[new_loc].emplace_back(
		deme,
		dispersal_abundance(
			deme.amount,
			deme.niche_suitability(conf, target_env),
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
  return (1.0f / (1.0f + exp(((distance /
			       conf.gene_flow_zero_distance) - B) * A)));
}


float Model::genetic_distance(const Deme &d1, const Deme &d2) {
  // Euclidean distance in gene space
  float sum = 0.0f;
  for (int i = 0; i < conf.genetic_dims; ++i) {
    const float d = d2.genetics.genetic_position[i] -
      d1.genetics.genetic_position[i];
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
  EnvCell ec;
  for (auto cell_itr = dm.begin(); cell_itr != dm.end(); ) {

    auto &&deme_list = cell_itr->second;
    if (deme_list.size() < 1) {
      // no demes, so get rid of the cell // FIXME can this actually happen??
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
    env.get(loc, ec); // known location so no need to check for no_data
    first_deme->amount = first_deme->niche_suitability(conf, ec);
    if (first_deme->amount <= 0.0f)
      // another extinction case
      cell_itr = dm.erase(cell_itr);
    else
      ++cell_itr;
  }
}

void Model::write_phylogeny(FILE *of, const Species::Vec &sv) {
  for (auto &&species: sv) {
    species->as_yaml(of, step, env.current_delta);
    write_phylogeny(of, species->sub_species);
  }
}

void Model::save() {
  std::string ofn = std::to_string(step) + ".csv";
  FILE *of = open_output_file(conf, ofn);
  fprintf(of, "species, row, column, amount");
  for (int i=0; i < conf.env_dims; ++i)
    fprintf(of, ", env_%d, niche_centre_%d, niche_breadth_%d", i, i, i);
  for (int i=0; i < conf.genetic_dims; i++)
    fprintf(of, ", genetic_position_%d", i);
  fprintf(of, "\n");

  EnvCell ec;
  bool no_data;
  char s1[16], s2[16];
  auto write_f = [&](float f) {
    // Write the more compact representation of a float.
    // e.g. 1234 instead of 1.23e+03
    int n1 = snprintf(s1, sizeof s1, "%g", f);
    int n2 = snprintf(s2, sizeof s2, "%0.*g", conf.csv_precision, f);
    return fprintf(of, ",%s",  n1 < n2? s1: s2);
  };

  for (auto &&species: tips) {
    for (auto &&deme_cell: *(species->demes)) {
      const Location &loc = deme_cell.first;
      env.get(loc, ec, &no_data);
      for (auto &&deme: deme_cell.second) {
	auto &&g = deme.genetics;
	fprintf(of,
		"%d,%d,%d",
		species->id, loc.y, loc.x);
	write_f(deme.amount);
	for (int i=0; i < conf.env_dims; ++i) {
	  write_f(no_data? NAN : ec[i]);
	  write_f(g.niche_centre[i]);
	  write_f(g.niche_tolerance[i] *2.0f);
	}
	for (int i=0; i < conf.genetic_dims; i++)
	  write_f(g.genetic_position[i]);
	if (fprintf(of, "\n") <1)
	  throw ApplicationException("Error writing " + ofn);
      }
    }
  }
  if (fclose(of) != 0)
    throw ApplicationException("Error closing " + ofn);

  ofn = std::to_string(step) + "-stats.yml";
  of = open_output_file(conf, ofn);
  fprintf(of, "# Step %d\n", step);
  write_phylogeny(of, roots);
  if (fclose(of) != 0)
    throw ApplicationException("Error closing " + ofn);
}


int Model::do_step() {
  /**
   * Execute one time step of the model
   */

  env.update(step); // step is 0-based here
  ++step;
  // step is 1-based below

  if (conf.verbosity > 0) {
    std::time_t now = std::time(nullptr);
    char tstr[50];
    if (!std::strftime(tstr, sizeof(tstr),
		       "%Y %b %d %H:%M:%S",
		       std::localtime(&now)))
      tstr[0] = '\0';
    std::cout << (step>1? "\n": "") << "Starting step " <<
      step << " at " << tstr << std::endl;
  }

  int n_occupied = 0;
  Species::Vec new_tips;
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

    if (conf.check_speciation && (step % conf.check_speciation == 0))
      species->speciate(step);
    if (species->sub_species.empty()) {
      // no speciation
      int n = species->update_stats(species->latest_stats, step);
      if (n > 0) {
        n_occupied += n;
        new_tips.push_back(species);
      }
    } else {
      // species has split into sub_species
      for (auto &&s: species->sub_species) {
        n_occupied += s->update_stats(s->initial_stats, step);
        s->latest_stats = s->initial_stats;
        new_tips.push_back(s);
      }
    }
  }

  tips = new_tips;
  save();
  return n_occupied;
}
