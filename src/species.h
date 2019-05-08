// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_SPECIES_H
#define DREADDS_SPECIES_H
#include <stdio.h>

#include <memory>
#include <vector>
#include <string>
#include <forward_list>
#include <utility>

#include <boost/serialization/array_wrapper.hpp> // just for boost 1.64. see https://svn.boost.org/trac10/ticket/12982
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "deme.h"

namespace DreadDs {

  namespace ba = boost::accumulators;

  template <class ...Types>
  using Acc = ba::accumulator_set<
    float,
    boost::accumulators::stats<Types...>>;

  typedef int Timestep;

  class SpeciesPresence;
  // SpeciesPresenceList stores entries in descending species id
  // order. This allows new species from speciation to be pushed on
  // the front.
  typedef std::forward_list<SpeciesPresence> SpeciesPresenceList;

  struct DispersalWeight {
    // Describes dispersal propensity for (x,y) offset from origin cell
    // at (0,0) due to distance cost.
    int x;
    int y;
    float weight; // 0 to 1.0
  };

  typedef std::vector <DispersalWeight> DispersalKernel;


  typedef std::pair<SpeciesPresenceList &,
		    SpeciesPresence &> SpeciationItem;
  typedef std::vector<SpeciationItem> SpeciationItemVector;

  class SpeciationCandidates {
  private:
    int count = 0;
    std::forward_list<SpeciationItem> candidates;
  public:
    void add(SpeciesPresenceList &spl, SpeciesPresence &sp);
    std::vector<SpeciationItemVector> split(float distance);
  };


  class Species {
  public:
    /**
       Describes a species and its phylogeny
    */
    typedef std::vector<std::shared_ptr <Species>> Vec;
    class Characteristics {
    public:
      // Derived from demes of this species.
      struct NicheStats {
	/**
	 * Describes a niche on an environmental variable for a species.
	 */
	// mean and sd of niche position of all demes of this species
	float position_mean = 0.0f;
	float position_sd = 0.0f;
	float breadth_mean = 0.0f;
	float breadth_sd = 0.0f;
	// max and min values of position ± tolerance
	float max = 0.0f;
	float min = 0.0f;
      };

      struct GeneticStats {
	/**
	   Holds the genetic position (on an abstract genetic trait) and
	   variance of all the demes of a species.
	*/
	float mean = 0.0f;
	float sd = 0.0f;
      };

      NicheStats niche_stats[max_env_dims];
      GeneticStats genetic_stats[max_genetic_dims];
      int cell_count; // number of demes (occupied cells).
      float population; // total population across all occupied cells
      int step = -1; // time step of last stats update
    };

    class StatsAccumulator {
    public:
      int count = 0;
      double population = 0.0;

      StatsAccumulator(const Config &c):
	niche_pos_acc(c.env_dims),
	niche_tol_acc(c.env_dims),
	genetic_acc(c.genetic_dims),
	niche_min(c.env_dims),
	niche_max(c.env_dims)
      {}

      void accumulate(const SpeciesPresence &sp);
      int update_stats(Characteristics &ch, int step);
    private:
      std::vector <Acc <ba::tag::mean, ba::tag::variance>>
      niche_pos_acc, niche_tol_acc, genetic_acc;
      std::vector <Acc <ba::tag::min>> niche_min;
      std::vector <Acc <ba::tag::max>> niche_max;
    };

    const Config &conf;
    Vec sub_species;
    Species *parent = NULL; // no need for smart pointer here as
                            // sub-species will always have a parent
    Timestep extinction = -1; // Time step of extinction, or -1 if extant
    Timestep split = -1; // Time step of speciation. -1 if we have no
                         // sub-species. This implies parent->split is
                         // species origin time.
    Characteristics latest_stats;  // Updated after each time step.
    DispersalKernel dk;
    int id = -1; // unique id starting at 0
    int family_id = -1; // species id of top level ancestor
    int step = -1; // time step of most recent dispersal

    Species(const Config &conf, const int species_id);
    Species(Species *parent, const int species_id);

    void setup_dispersal(const SpeciesParameters &sp);
    void update(int step, StatsAccumulator &acc);
    void speciate(int *id_counter, SpeciationCandidates &candidates);

    void as_yaml(FILE *of,
                 const std::string &first_indent);
    void phylogeny_as_yaml(FILE *of,
                           const std::string &first_indent);
    std::string phylogeny_as_newick(); // See https://en.wikipedia.org/wiki/Newick_format
    std::string get_name();
    void log_summary_stats(const Characteristics &ch);

    inline int get_id() {
      return id +1; // For 1-based indexing in output files and R API
    }

    void set_initial_stats(StatsAccumulator &acc) {
      acc.update_stats(initial_stats, step);
      latest_stats = initial_stats;
      if (conf.verbosity > 1)
	log_summary_stats(initial_stats);
    }
    void update_latest_stats(StatsAccumulator &acc) {
      acc.update_stats(latest_stats, step);
      if (conf.verbosity > 1)
	log_summary_stats(latest_stats);
    }

  private:
    Characteristics initial_stats; // At species origin (i.e. split from parent)
    std::shared_ptr <Species> add_child(const int species_id);
    std::vector <SpeciationItemVector> get_clusters(float distance,
						    SpeciationCandidates &speciation_candidates);
  };


  /**
   * Describes the presence of a species in a cell
   */
  class SpeciesPresence {
  public:
    Species *species;
    Deme incumbent; // holds the population that exists due to prior occupation
    std::vector<Deme> immigrants; // During dispersal, several demes
				  // can occupy a cell

    SpeciesPresence(Species *s):
      incumbent(-1.0f), species(s) {}

    SpeciesPresence(Species *s, const Deme &d):
      incumbent(d), species(s) {}

  };

}
#endif
