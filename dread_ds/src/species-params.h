// -*- coding: utf-8; Mode: c++  -*-

#ifndef DREADDS_SPECIES_PARAMS_H
#define DREADDS_SPECIES_PARAMS_H

#include "model-limits.h"

namespace DreadDs {

  struct SpeciesParameters {
    struct Niche {
      /**
	 Describes a niche on an environmental variable for a species.
      */
      // mean and sd of niche position of all demes of this species
      float position_mean = 0.0f;
      float position_sd = 0.0f;
      float breadth_mean = 0.0f;
      float breadth_sd = 0.0f;
      // max and min values of position - (breadth  /2)
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

    Niche niche[max_env_dims];   // Niches derived from all demes of this species.
    Genetics genetics[max_genetic_dims];
  };
}

#endif
