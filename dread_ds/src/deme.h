// -*- coding: utf-8; Mode: c++  -*-
#ifndef DREADDS_DEME_H
#define DREADDS_DEME_H

/**
 * A "deme" is a population in a cell
 */


#include <iostream> // FIXME debugging

#include <list>
#include <map>

#include "model-limits.h"
#include "model-config.h"

namespace DreadDs {

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

    class Genetics {
    public:
      float niche_centre[max_env_dims];
      float niche_tolerance[max_env_dims];
      float genetic_position[max_genetic_dims]; // genetic position in n-dimensional space. See struct Genetics

      Genetics() {}
      Genetics(float v):
	niche_centre {v},
	niche_tolerance {0.0f},
	genetic_position {0.0f} {

	  // FIXME debugging - just checking brace init ok for these arrays
	  std::cout << "\n\ngenetics:\n";
	  for (int i=0; i<max_env_dims; ++i)
	    std::cout << niche_centre[i] << " " << niche_tolerance[i] << " " << genetic_position[i] << std::endl;


	}
      Genetics(const SpeciesParameters &sp);
      float niche_suitability(const Config &conf, const float *env);
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
}
#endif
