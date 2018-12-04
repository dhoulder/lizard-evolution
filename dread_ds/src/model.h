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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "deme.h"
#include "species.h"

namespace DreadDs {

  typedef boost::random::mt19937 rng_eng_t;
  typedef boost::random::uniform_real_distribution<float> uniform_distr_t;
  typedef boost::random::variate_generator<rng_eng_t&, uniform_distr_t> uniform_vg_t;
  typedef boost::random::normal_distribution<float> normal_distr_t;
  typedef boost::random::variate_generator<rng_eng_t&, normal_distr_t> normal_vg_t;

  class Model {
  public:
    const Config conf; // must be first (initialisation order)
    int step = 0;
    Environment env;
    std::vector <std::shared_ptr <Species>> roots; // Initial species
    std::vector <std::shared_ptr <Species>> tips; // extant leaf species

    Model(const Config &conf);

    /**
     * Execute one time step of the model.
     * Returns: total number of occupied cells across all species.
     */
    int do_step();

  private:
    uniform_distr_t gene_flow_distr;
    uniform_vg_t gene_flow_random;
    uniform_distr_t deme_choice_distr;
    normal_distr_t gene_drift_distr;
    normal_vg_t gene_drift_random;

    void evolve_towards_niche(Deme &deme, const EnvCell &env);
    void do_genetc_drift(Deme &deme);
    Deme *choose_primary(DemeList &deme_list);
    float gene_flow_probability(float distance);
    float genetic_distance(const Deme &d1, const Deme &d2);
    bool gene_flow_occurs(const Deme &d1, const Deme &d2);
    std::shared_ptr<DemeMap> evolve_and_disperse(Species &species);
    void merge(DemeMap &dm);
    void save();
  };
}
#endif
