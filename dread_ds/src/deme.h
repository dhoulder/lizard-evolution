// -*- coding: utf-8; Mode: c++ -*-
#ifndef DREADDS_DEME_H
#define DREADDS_DEME_H

/**
 * A "deme" is a population in a cell
 */

#include <list>
#include <map>

#include "model-limits.h"
#include "model-config.h"
#include "environment.h"

namespace DreadDs {

  class Deme {
    /**
       Describes a genetically homogeneous population in a cell.
    */
  public:

    class Genetics {
    public:
      float niche_centre[max_env_dims];
      float niche_tolerance[max_env_dims];
      float genetic_position[max_genetic_dims]; // genetic position in n-dimensional space. See struct Genetics

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
      Deme(from) {
      // FIXME only need to do env_dims, not max_env_dims. Maybe use float xxxx[nnnn] = {0.0f} just to be sure? or pass in nnnn
      amount = new_amount;
      is_primary = new_primary;
    }

    float niche_suitability(const Config &conf, const EnvCell &env);

  };

  // Cells occupied by demes of a species, Several demes can occupy a
  // cell, hence std::list
  typedef std::list <Deme> DemeList;
  typedef std::map <Location, DemeList> DemeMap;
}
#endif
