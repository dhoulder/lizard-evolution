// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_SPECIES_H
#define DREADDS_SPECIES_H

#include <memory>
#include <vector>

#include <iostream> // FIXME for debugging

#include "model-limits.h"
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

    class Characteristics {
    public:
      // Derived from demes of this species.

      struct Range {
	int cell_count; // number of demes (occupied cells).
	float population; // total population across all occupied cells
      };

      struct Niche {
	/**
	 * Describes a niche on an environmental variable for a species.
	 */
	// mean and sd of niche position of all demes of this species
	float position_mean = 0.0f;
	float position_sd = 0.0f;
	float breadth_mean = 0.0f;
	float breadth_sd = 0.0f;
	// max and min values of position - (breadth /2)
	float max = 0.0f;
	float min = 0.0f;
      };

      struct Genetics {
	/**
	   Holds the genetic position (on an abstract genetic trait) and
	   variance of all the demes of a species.
	*/
	float position = 0.0f;
	float variance = 0.0f;
      };

      Niche niche[max_env_dims];
      Genetics genetics[max_genetic_dims];
      Range range;

      void update(const Config &conf, const DemeMap &demes);
    };

    std::shared_ptr <Species> left_child = NULL;
    std::shared_ptr <Species> right_child = NULL;
    std::weak_ptr <Species> parent;

    Timestep extinction = -1; // Time step of extinction, or -1 if extant
    Timestep split = -1; // Time of speciation. parent->split is species
    // origin time. -1 if not speciated

    Characteristics initial_stats; // At species origin (i.e. split from parent)
    Characteristics latest_stats;  // Updated after each time step.

    std::shared_ptr <DemeMap> demes; // Cells occupied by this species.

    DispersalKernel dk;

    Species(const Config &conf, const SpeciesParameters &sp, const Environment &env);

    void print_kernel() { // FIXME debugging
      std::cout << "Dispersal kernel" << std::endl;
      for (auto &&v: dk)
	std::cout << v.x << ", " << v.y << " " << v.weight << std::endl;
    }

    inline void update_stats(const Config &conf) {
      latest_stats.update(conf, *demes);
    }


  private:
    void setup_dispersal(float dispersal_min, const SpeciesParameters &sp);
    void load_initial(const Config &conf, const SpeciesParameters &sp, const Environment &env);
  };

}
#endif
