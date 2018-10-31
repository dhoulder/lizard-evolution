// -*- coding: utf-8; Mode: c++ -*-

/**
 * Internal model declarations
 */

#ifndef DREADDS_MODEL_H
#define DREADDS_MODEL_H

#define BOOST_DISABLE_ASSERTS

#include <memory>
#include <vector>
#include <string>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "model-limits.h"
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
    const Config conf;
    Environment env;
    std::vector <std::shared_ptr <Species>> roots; // Initial species
    std::vector <std::shared_ptr <Species>> tips; // extant leaf species
    std::string output_path;

    Model(const char *config_path,
	  const char *output_path);

    Model(const Config &conf,
	  const char *output_path);


    ~Model() {
    //FIXME WIP
    }

    /**
     * Execute one time step of the model.
     * Returns: time step just executed. First step is 1.
     */
    int do_step();

  private:
    int step = 0;
    float env_delta[max_env_dims] = {0.0f};
    uniform_distr_t gene_flow_distr;
    uniform_vg_t gene_flow_random;
    uniform_distr_t deme_choice_distr;
    normal_distr_t gene_drift_distr;
    normal_vg_t gene_drift_random;

    EnvCell get_env(const Location &loc);
    void evolve_towards_niche(Deme &deme, const EnvCell &env);
    void do_genetc_drift(Deme &deme);
    Deme *choose_primary(DemeList &deme_list);
    float gene_flow_probability(float distance);
    float genetic_distance(const Deme &d1, const Deme &d2);
    bool gene_flow_occurs(const Deme &d1, const Deme &d2);
    void update_environment(int time_step);
    std::shared_ptr<DemeMap> disperse(Species &species);
    void merge(DemeMap &dm);
  };
}
#endif
