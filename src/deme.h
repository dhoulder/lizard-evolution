// -*- coding: utf-8; Mode: c++ -*-
#ifndef DREADDS_DEME_H
#define DREADDS_DEME_H

/**
 * A "deme" is a population in a cell
 */
#include <cmath>

#include "constants.h"
#include "model-config.h"
#include "environment.h"

namespace DreadDs {

  class Deme {
    /**
     *  Describes a genetically homogeneous population in a cell.
     */
  public:

    class Genetics {
    public:
      float niche_centre[max_env_dims];
      float niche_tolerance[max_env_dims];
      float genetic_position[max_genetic_dims]; // genetic position in
      // n-dimensional space.

      Genetics():
        // zero everything
        niche_centre {},
        niche_tolerance {},
        genetic_position {} {}
      Genetics(const SpeciesParameters &sp);
    };

    Genetics genetics;
    float amount; // population per cell

    Deme(float _amount = 0.0f): amount(_amount) {}

    Deme(const SpeciesParameters &sp):
      amount(0.0f),
      genetics(sp) {}

    Deme(const Deme &from, float new_amount):
      // Delegate to implicit copy constructor to copy genetics
      Deme(from) {
      amount = new_amount;
    }

    float niche_suitability(const Config &conf, const EnvCell env);

    inline float genetic_distance_sq(const Config &conf, const Deme &d2) const {
      // Square of Euclidean distance in gene space
      float sum = 0.0f;
      for (int i = 0; i < conf.genetic_dims; ++i) {
	const float d = d2.genetics.genetic_position[i] -
	  genetics.genetic_position[i];
	sum += d*d;
      }
      return sum;
    }

    inline float genetic_distance(const Config &conf, const Deme &d2) const {
      return std::sqrt(genetic_distance_sq(conf, d2));
    }

  };
}
#endif
