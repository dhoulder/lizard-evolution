// -*- coding: utf-8; Mode: c++ -*-

/**
 * Internal model declarations
 */

#ifndef DREADDS_MODEL_H
#define DREADDS_MODEL_H

#include <memory>
#include <vector>
#include <string>

#define BOOST_DISABLE_ASSERTS

#include "random.h"
#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "deme.h"
#include "species.h"

namespace DreadDs {

  class Model {
  public:
    const Config conf; // must be first (initialisation order)
    int step = 0;
    int species_id_counter = 0;
    Environment env;
    Species::Vec roots; // Initial species
    Species::Vec tips; // Extant leaf species

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

  private:
    rng_eng_t rng;
    uniform_real_distr_t gene_flow_distr;
    uniform_real_distr_t deme_choice_distr;
    normal_distr_t gene_drift_distr;

    void evolve_towards_niche(Deme &deme, const EnvCell env);
    void do_genetc_drift(Deme &deme);
    Deme *choose_primary(DemeList &deme_list);
    float gene_flow_probability(float distance);
    bool gene_flow_occurs(const Deme &d1, const Deme &d2);
    std::shared_ptr<DemeMap> evolve_and_disperse(Species &species);
    void merge(DemeMap &dm);
  };
}
#endif
