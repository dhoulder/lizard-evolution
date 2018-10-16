// -*- coding: utf-8; Mode: c++  -*-

/**
 * Internal model declarations
 */

#ifndef DREADDS_MODEL_H
#define DREADDS_MODEL_H

#define BOOST_DISABLE_ASSERTS

#include <memory>
#include <vector>
#include <map>
#include <list>

#include "boost/multi_array.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "model-config.h"

namespace DreadDs {

  typedef boost::random::mt19937 rng_eng_t;
  typedef boost::random::uniform_real_distribution<float> uniform_distr_t;
  typedef boost::random::variate_generator<rng_eng_t&, uniform_distr_t> uniform_vg_t;
  typedef boost::random::normal_distribution<float> normal_distr_t;
  typedef boost::random::variate_generator<rng_eng_t&, normal_distr_t> normal_vg_t;

  struct EnvCell {
    // has to be a class|struct to keep boost::multi_array allocator happy
    float v[max_env_dims];
  };

  typedef int Timestep;
  typedef boost::multi_array<EnvCell, 2> EnvMatrix;
  typedef EnvMatrix::index EnvIndex;

  struct Location {
    int x;
    int y;

    // Used as key in a map, so needs an ordering
    friend bool operator< (const Location &a, const Location &b) {
      return (a.x < b.x) || (a.x == b.x && a.y < b.y);
    }
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

  class Species {
  public:
    /**
       Describes a species and its phylogeny
    */

    struct Range {
      int cell_count; // number of demes (occupied cells).
      float population; // total population across all occupied cells
    };

    const float max_dispersal_radius_;
    const float dispersal_min_;

    Species *parent = NULL;
    std::shared_ptr <Species> left_child = NULL;
    std::shared_ptr <Species> right_child = NULL;

    Timestep extinction = -1; // Time step of extinction, or -1 if extant
    Timestep split = -1; // Time of speciation. parent->split is species
    // origin time. -1 if not speciated
    struct {
      SpeciesParameters params;
      Range range;
    } initial, // At species origin (i.e. split from parent)
      latest; // Updated at each time step. Frozen after speciation.

    std::shared_ptr <DemeMap> demes; // Cells occupied by this species.

    DispersalKernel dk;

    Species(float max_dispersal_radius,
	    float dispersal_min = 0.2f);

    void print_kernel() { // FIXME debugging
      for (auto &&v: dk)
	std::cout <<  v.x << ", " << v.y <<  " " << v.weight << std::endl;
    }

  private:
    void setup_dispersal();
  };


  class Model {
  public:
    const Config conf;
    EnvMatrix env;
    std::vector <Species> roots; // Initial species
    std::vector <Species> tips; // extant leaf species

    Model(const char *config_file);

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
    float niche_suitability(const EnvCell &cell, const Deme::Genetics &g);
    void evolve_towards_niche(Deme &deme, const EnvCell &env);
    void do_genetc_drift(Deme &deme);
    Deme *choose_primary(DemeList &deme_list);
    float gene_flow_probability(float distance);
    float genetic_distance(const Deme &d1, const Deme &d2);
    bool gene_flow_occurs(const Deme &d1, const Deme &d2);
    void update_environment(int time_step);
    std::shared_ptr<DemeMap> disperse(Species &species);
    void merge(std::shared_ptr<DemeMap> dm);
  };
}
#endif
