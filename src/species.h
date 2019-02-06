// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_SPECIES_H
#define DREADDS_SPECIES_H
#include <stdio.h>

#include <memory>
#include <vector>
#include <string>

#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "deme.h"


namespace DreadDs {

  typedef int Timestep;

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
    typedef std::vector <std::shared_ptr <Species>> Vec;

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

    const Config &conf;
    Vec sub_species;

    Species *parent = NULL; // no need for smart pointer here as
                            // sub-species will always have a parent
    Timestep extinction = -1; // Time step of extinction, or -1 if extant
    Timestep split = -1; // Time step of speciation. -1 if we have no
                         // sub-species. This implies parent->split is
                         // species origin time.
    Characteristics latest_stats;  // Updated after each time step.
    std::shared_ptr <DemeMap> demes; // Cells occupied by this species.
    DispersalKernel dk;
    int id = 0;
    int step = -1; // time step of most recent dispersal

    Species(const Config &conf);
    Species(const Config &conf,
            const SpeciesParameters &sp,
            const Environment &env);
    void set_initial_stats();
    void update(const std::shared_ptr <DemeMap> &d, int step);
    void speciate();
    int update_stats(Characteristics &ch);
    void as_yaml(FILE *of,
                 const std::string &first_indent);
    void phylogeny_as_yaml(FILE *of,
                           const std::string &first_indent);

  private:
    Characteristics initial_stats; // At species origin (i.e. split from parent)
    void setup_dispersal(const SpeciesParameters &sp);
    void load_initial(const SpeciesParameters &sp, const Environment &env);
  };

}
#endif
