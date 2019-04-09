// -*- coding: utf-8; Mode: c++ -*-
#ifndef DREADDS_DEME_H
#define DREADDS_DEME_H

/**
 * A "deme" is a population in a cell
 */
#include <cmath>
#include <list>
#include <map>

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
    bool is_primary; // indicates incumbency in a cell during dispersal

    Deme(): amount(0.0f), is_primary(false) {}

    Deme(const SpeciesParameters &sp):
      amount(0.0f), is_primary(false),
      genetics(sp) {}

    Deme(const Deme &from, float new_amount, bool new_primary):
      // Delegate to implicit copy constructor to copy genetics
      Deme(from) {
      amount = new_amount;
      is_primary = new_primary;
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

  // Cells are occupied by demes of a species, Several demes can
  // occupy a cell, hence std::list
  typedef std::list <Deme> DemeList;
  // TODO DemeList could probably just be a forward_list and use
  // emplace_after(whatever.begin(), …)  and emplace_front(). Or maybe
  // struct {Deme primary; std::forward_list <Deme> others}, or maybe
  // just use a vector and reserve(worst_case_length). Do we need to
  // handle more than 1 primary deme?
  typedef std::map <Location, DemeList> DemeMap;
  // Note that compared to a simple 2D array, DemeMap is faster to
  // iterate through when the population is sparse, but slower to
  // index by location. Profiling indicates that around 20% to 30% of
  // the total model runtime is spent accessing or modifying
  // DemeMap()s. A 2D array with some kind of adaptive bounding box to
  // mask out many of the empty cells might be faster.
  typedef std::pair<const Location, DemeList> DemeMapEntry;
}
#endif
