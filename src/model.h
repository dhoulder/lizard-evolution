// -*- coding: utf-8; Mode: c++ -*-

/**
 * Internal model declarations
 */

#ifndef DREADDS_MODEL_H
#define DREADDS_MODEL_H

#include <memory>
#include <vector>
#include <string>
#include <map>

#define BOOST_DISABLE_ASSERTS

#include "random.h"
#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "deme.h"
#include "species.h"

namespace DreadDs {
  typedef std::map <Location, SpeciesPresenceList> DemeMap;
  // Note that compared to a simple 2D array, DemeMap is faster to
  // iterate through when the population is sparse, but slower to
  // index by location. Profiling indicates that around 20% to 30% of
  // the total model runtime is spent accessing or modifying
  // DemeMap()s. A 2D array with some kind of adaptive bounding box to
  // mask out many of the empty cells might be faster.

  // The Model() holds a single object that describes the dispersal
  // characteristics of all species in a way that can be iterated over
  // easily when we're doing dispersal.
  struct SpeciesDispersalWeight {
    float weight;
    int species_id;

    SpeciesDispersalWeight(int id, float w): species_id(id), weight(w) {}
  };

  class MergeResult {
  public:
    Species::StatsAccumulator acc;
    SpeciationCandidates sp_candidates;

    MergeResult(const Config &c):
      acc(c) {};
  };

  typedef std::vector<MergeResult> MergeResultVector;

  class Model {
  public:
    const Config conf; // must be first (initialisation order)
    int step = 0;
    Environment env;
    Species::Vec roots; // Initial species
    Species::Vec tips; // Extant leaf species
    std::shared_ptr <DemeMap> demes; // Cells occupied by all species.

    Model(const Config &conf);

    /**
     * Execute one time step of the model.
     * Returns: total number of occupied cells across all species.
     */
    int do_step();

    /**
     * Save current state to files in output directory
     */
    void save();

    /**
     * Return all species known to the model as a flattened list
     */
    Species::Vec get_all_species();

    inline int get_species_count() const {
      return species_counter;
    }

  private:
    int species_counter = 0;
    rng_eng_t rng;
    uniform_real_distr_t gene_flow_distr;
    uniform_real_distr_t deme_choice_distr;
    normal_distr_t gene_drift_distr;

    void load_initial_species(const SpeciesParameters &sp);
    void evolve_towards_niche(Deme &deme, const EnvCell env);
    void do_genetc_drift(Deme &deme);
    Deme &choose_primary(std::vector<Deme> &immigrants);
    float gene_flow_probability(float distance);
    bool gene_flow_occurs(const Deme &d1, const Deme &d2);
    void evolve_and_disperse();
    MergeResultVector merge(bool do_speciation);
  };
}
#endif
