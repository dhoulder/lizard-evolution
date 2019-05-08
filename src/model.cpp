// -*- coding: utf-8 -*-

/**
 * Implementation of DREaD_ds model. See ./README.md
 */

#include <stdio.h>

#include <ctime>
#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "deme.h"
#include "species.h"
#include "model.h"
#include "exceptions.h"
#include "output-file.h"

using namespace DreadDs;
using boost::algorithm::replace_all_copy;

namespace DreadDs {

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
      return true;
    }
  };
}

Model::Model(const Config &c):
  conf(c),
  step(0),
  demes(std::make_shared <DemeMap>()),
  species_id_counter(0),
  env(Environment(conf)), // must use conf member, not arg

  // Generate random values between [gene_flow_clip,
  // 1-gene_flow_clip] to effectively apply high and low
  // cutoffs when comparing against probability below.
  gene_flow_distr(
		  conf.gene_flow_clip,
		  1.0f - conf.gene_flow_clip),
  deme_choice_distr(0.0f, 1.0f),
  gene_drift_distr(0.0f, gene_drift_sd)
{
  // Seed our RNG from the already randomly seeded Config::rng
  rng.seed(c.rng());

  for (auto &&sp: conf.get_initial_species(env))
    load_initial_species(sp);
}


/**
 * Load initial species values, locations and abundance.
 */

void Model::load_initial_species(const SpeciesParameters &sp) {
  const auto &&species = std::make_shared<Species>(conf,
					     // must be loaded in increasing id order
					     ++species_id_counter);
  species->setup_dispersal(sp);

  // get initial bounding rectangle row and column limits.

  // sp.north…west have already been bounds checked in sp.get_initial_species()
  long n = env.row(sp.north);
  long s = env.row(sp.south);
  long e = env.col(sp.east);
  long w = env.col(sp.west);

  if (conf.verbosity > 1)
    std::cout << " Loading species. "
      "west=" << sp.west <<
      ", east=" << sp.east <<
      ", south=" << sp.south <<
      ", north=" << sp.north  << std::endl;

  // fill all suitable cells in initial-rectangle

  Deme d(sp);
  Location loc;
  EnvCell ec;
  bool no_data;
  bool absent = true;;
  Species::StatsAccumulator acc(conf);
  // roots and tips must hold species in ascending id order
  roots.push_back(species);

  for (loc.y = n; loc.y <= s; ++loc.y)
    for (loc.x = w; loc.x <= e; ++loc.x) {
      env.get(loc, ec, &no_data);
      if (no_data)
	continue;
      d.amount = d.niche_suitability(conf, ec);
      if (d.amount > 0.0f) {
	absent = false;
	// push_front() so they end up in descending id order. Note
	// that species is now in roots so we can safely store a plain
	// pointer to it for the duration of the model.
	auto &&sp = SpeciesPresence(&(*species), d);
	acc.accumulate(sp);
	(*demes)[loc].push_front(std::move(sp)); // these end up in descending id order.
      }
    }

  species->set_initial_stats(acc);
  if (absent)
    species->extinction = 0; // Dead on arrival

  if (species->extinction < 0)
    tips.push_back(species);
}


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
    *gp += gene_drift_distr(rng);
}


static inline SpeciesPresence &get_target_presence(SpeciesPresenceList &target_cell,
                                                   SpeciesPresenceList::iterator &target_presence_itr,
                                                   Species *species) {
  // Insert or find list element. They're in descending species.id order in target_cell
  for (;;) {
    auto itr = target_presence_itr; // note: may be before_begin(), so not safe to dereference
    ++target_presence_itr;

    if (target_presence_itr == target_cell.end() || target_presence_itr->species->id < species->id) {
      // append to list or insert before next species
      target_presence_itr = target_cell.emplace_after(itr, species);
      break;
    }
    if (target_presence_itr->species->id == species->id)
      break; // found existing SpeciesPresence
  }

  return *target_presence_itr;
}


void Model::evolve_and_disperse() {
  /**
   * Iterate over all extant demes of all species and disperse to
   * neighbouring cells, weighted by environmental niche
   * suitability.
   */
  EnvCell source_env, target_env;
  bool no_data, no_target;
  auto target = std::make_shared <DemeMap>();
  auto &&env_shape = env.values.shape();

  // Collect all the dispersal kernels into a map so that we can
  // iterate over the dispersal values for all species at a given
  // location.
  std::map<
    Location,
    std::vector<SpeciesDispersalWeight> // must be in ascending species id order.
    >  aggregated_dk;
  for (const auto &s: tips)
    for (const auto &dw: s->dk)
      aggregated_dk[Location(dw.x, dw.y)].emplace_back(s->id,
						      dw.weight);
  // TODO? collapse aggregated_dk to vector for faster iteration?

  // Iterate over all cells
  for (auto &dememap_entry: *demes) {
    const Location &loc = dememap_entry.first;
    env.get(loc, source_env, &no_data);
    if (no_data)
      continue;

    // "disperse" into the same cell. These become primary demes after
    // dispersal due to incumbency.
    SpeciesPresenceList &source_cell = dememap_entry.second;
    SpeciesPresenceList &target_cell = (*target)[loc];
    auto target_presence_itr = target_cell.before_begin();

    for (SpeciesPresence &presence: source_cell) {
      // evolve and disperse to same location. source_cell holds species in decreasing id order
      auto &&source_deme = presence.incumbent;
      evolve_towards_niche(source_deme, source_env);
      do_genetc_drift(source_deme);

      Deme &target_deme = get_target_presence(target_cell,
                                              target_presence_itr,
                                              presence.species).incumbent;
      target_deme = source_deme;
      target_deme.amount = dispersal_abundance(
		    source_deme.amount,
		    source_deme.niche_suitability(conf, source_env),
		    1.0); // same cell, no travel cost
    }

    // Loop over locations covered by dispersal kernels and disperse each species

    for (const auto &adk_entry: aggregated_dk) {
      const Location &offset = adk_entry.first;
      Location new_loc(loc.x + offset.x, loc.y + offset.y);
      if (new_loc.x < 0 || new_loc.y < 0 ||
	  new_loc.x >= env_shape[1] || new_loc.y >= env_shape[0])
	continue;
      env.get(new_loc, target_env, &no_target);
      if (no_target)
	continue;

      auto dw_itr = adk_entry.second.crbegin();  // iterating backwards, largest id first
      SpeciesPresenceList &target_cell = (*target)[new_loc];
      auto target_presence_itr = target_cell.before_begin();
      for (SpeciesPresence &presence: source_cell) {
	// Not every species in adk_entry may be present in this cell,
	// and not every species in the cell may be present in
	// adk_entry. Species are ordered by ascending id in adk_entry and by
	// descending id in source_cell, so skip ahead until we find
	// the dispersal weight for this species.
        for (;dw_itr !=  adk_entry.second.crend(); ++dw_itr) {
          if (dw_itr->species_id == presence.species->id) {
            // Got a dispersal weight for this species in this
            // cell. Disperse to this species in the target cell
            auto &source_deme = presence.incumbent;
            auto &immigrants = get_target_presence(target_cell,
                                                   target_presence_itr,
                                                   presence.species).immigrants;
            if (immigrants.size() == 0) // first time
              immigrants.reserve(presence.species->dk.size());
            immigrants.emplace_back(
		  source_deme,
		  dispersal_abundance(
		       source_deme.amount,
		       source_deme.niche_suitability(conf, target_env),
		       dw_itr->weight));
            ++dw_itr; // to disperse next species
            break;
          }
          if (dw_itr->species_id < presence.species->id)
	    // No dispersal for this species at this location.
	    // Skip to next species
            break;
	  // Otherwise the species for this dispersal weight isn't
	  // present in this cell. Loop to next dispersal weight
        }
      }
    }
  }

  // Finished dispersing so replace source cells
  demes = target;
}


Deme &Model::choose_primary(std::vector<Deme> &immigrants) {
  // Choose a "primary" deme at random from immigrant demes in a cell using
  // abundance-weighted probability. Do not call with an empty vector.
  float sum = 0.0f;;
  for (auto &&d: immigrants)
    sum += d.amount;
  float rv = deme_choice_distr(rng) * sum;
  for (auto &&d: immigrants) {
    if (rv <= d.amount )
      return d;
    rv -= d.amount;
  }
  // Very occasionally fp precision effects may leave us here
  return immigrants.back();
}


float Model::gene_flow_probability(float distance) {
  // ~1.0 for distance == 0, approaches 0 for distance >= gene_flow_zero_distance
  const float A = 14.0f;
  const float B = 0.5f;
  return (1.0f / (1.0f + std::exp(((distance /
				    conf.gene_flow_zero_distance) - B) * A)));
}


bool Model::gene_flow_occurs(const Deme &d1, const Deme &d2) {
  // see 3.4.1 "Does gene flow occur?"
  return (gene_flow_probability(d1.genetic_distance(conf, d2)) >
	  gene_flow_distr(rng));
}


/**
 * Merge immigrant demes into incumbents.
 */
MergeResultVector Model::merge(bool do_speciation) {

  // Merge demes where gene flow occurs
  EnvCell ec;
  MergeResultVector mrv(get_species_count(), conf);

  for (auto cell_itr = demes->begin(); cell_itr != demes->end(); ) {
    const Location &loc = cell_itr->first;
    env.get(loc, ec); // known location so no need to check for no_data

    SpeciesPresenceList &spl = cell_itr->second;
    for (auto sp_itr = spl.begin(), prev = spl.before_begin();
	 sp_itr != spl.end(); ) {
      Deme &incumbent = sp_itr->incumbent;

      if (!sp_itr->immigrants.empty()) {
	// Have at least two demes in this cell, so check for gene flow
	// and merge if required.

	// Any incumbent population will be the primary deme, otherwise
	// we have to choose one
	Deme &primary = (incumbent.amount < 0.0f) ?
	  choose_primary(sp_itr->immigrants):
	  incumbent;

	DemeMixer mixer(conf);
	mixer.contribute(primary);
	for (auto &&deme: sp_itr->immigrants)
	  if ((&deme != &primary) && gene_flow_occurs(deme, primary))
	    mixer.contribute(deme);
        // otherwise excluded by competition

	// mix down to primary deme
	if (!mixer.mix_to(&incumbent)) {
	  // extinct in this cell, so drop it
	  sp_itr = spl.erase_after(prev);
	  continue;
	}
	sp_itr->immigrants.clear();
      }
      // update abundance according to current environment
      incumbent.amount = incumbent.niche_suitability(conf, ec);


      /**
       * TODO: Implement species competition here. May need to do all
       * the immigration first and then loop over all the resulting
       * incumbents in order to simulate competition.  That will
       * probably involve some refactoring of the enclosing loops and
       * the code below.
       */


      if (incumbent.amount > 0.0f) {
	// Present.
	Species &species = *sp_itr->species;
	MergeResult &mr = mrv[species.id - 1];
	mr.acc.accumulate(*sp_itr); // contribute to statistics
	if (do_speciation)
	  mr.sp_candidates.add(spl, *sp_itr);
	// Advance to next species in cell
	prev = sp_itr;
	++sp_itr;
      } else
	// Species extinct here. drop from list
	sp_itr = spl.erase_after(prev);
    }

    if (spl.empty())
      // all species extinct in this cell, so get rid of it
      cell_itr = demes->erase(cell_itr);
    else
      ++cell_itr;
  }
  return mrv;
}


void Model::save() {
  std::string ofn = std::to_string(step) + ".csv";
  FILE *of = open_output_file(conf, ofn);
  fprintf(of, "species, row, column, amount");
  for (int i=0; i < conf.env_dims; ++i)
    fprintf(of, ", env_%d, niche_centre_%d, niche_breadth_%d", i, i, i);
  for (int i=0; i < conf.genetic_dims; i++)
    fprintf(of, ", genetic_position_%d", i);
  if (fprintf(of, "\n") < 1)
    throw ApplicationException("Error writing " + ofn);

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

  for (const auto &entry: *demes) {
    const Location &loc = entry.first;
    const SpeciesPresenceList &spl = entry.second;
    env.get(loc, ec, &no_data);
    for (const SpeciesPresence &presence: spl) {
      auto &&d = presence.incumbent;
      auto &&g = d.genetics;
      fprintf(of,
	      "%d,%d,%d",
	      presence.species->id, loc.y, loc.x);
      write_f(d.amount);
      for (int i=0; i < conf.env_dims; ++i) {
	write_f(no_data? NAN : ec[i]);
	write_f(g.niche_centre[i]);
	write_f(g.niche_tolerance[i] *2.0f);
      }
      for (int i=0; i < conf.genetic_dims; i++)
	write_f(g.genetic_position[i]);
      if (fprintf(of, "\n") < 1)
	throw ApplicationException("Error writing " + ofn);
    }
  }
  if (fclose(of) != 0)
    throw ApplicationException("Error closing " + ofn);

  ofn = std::to_string(step) + "-stats.yml";
  of = open_output_file(conf, ofn);
  fprintf(of,
	  "step: %d\n"
	  "input_environment:\n",
	  step);
  for (int i=0; i < conf.env_dims; ++i)
    fprintf(of,
            " - filename: '%s'\n"
            "   delta: %f\n",
            replace_all_copy(conf.env_params[i]->get_filename(step-1),
                             "'", "''").c_str(),
            env.current_delta[i]);
  if (fprintf(of, "species:\n") <1)
    throw ApplicationException("Error writing " + ofn);
  for (auto &&species: roots)
    species->phylogeny_as_yaml(of, " - ");

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

  bool do_speciation = conf.check_speciation && (step % conf.check_speciation == 0);

  evolve_and_disperse();
  // merge demes in each cell that are within genetic tolerance.
  MergeResultVector &&mrv = merge(do_speciation);

  // compute statistics. do speciation if required
  int n_occupied = 0;
  Species::Vec new_tips;
  for (auto &species: tips) {
    MergeResult &mr = mrv[species->id -1];
    species->update(step, mr.acc);
    if (species->extinction > -1)
      // Just went extinct in this step
      continue;

    n_occupied += mr.acc.count;
    if (do_speciation)
      species->speciate(&species_id_counter,
			mr.sp_candidates);
    if (species->sub_species.empty()) {
      // No speciation.
      // We have to update_stats() after every step as we may
      // speciate on the next iteration and we want to leave the
      // parent species with accurate stats.
      species->update_latest_stats(mr.acc);
      new_tips.push_back(species);
    } else {
      // species has split into sub_species

      if (conf.verbosity > 0)
	std::cout << "Speciation occurred" << std::endl;
      for (auto &&s: species->sub_species)
	new_tips.push_back(s);
    }
  }
  // tips must be in ascending id order so that aggregated_dk gets built properly
  std::sort(new_tips.begin(), new_tips.end(),
	    [] (std::shared_ptr <Species> a,
		std::shared_ptr <Species> b) {
	      return a->id < b->id;
	    });
  tips = new_tips;
  return n_occupied;
}


static void traverse_species(Species::Vec &flat,
                             std::shared_ptr <Species> &s) {
  flat.push_back(s);
  for (auto &&child: s->sub_species)
    traverse_species(flat, child);
}


Species::Vec Model::get_all_species()  {
  Species::Vec v;
  for (auto &&species: roots)
    traverse_species(v, species);
  return v;
}
