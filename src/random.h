// -*- coding: utf-8; Mode: c++ -*-

/**
 * Types for random number generation. These typedefs allow the random
 * number generation to be replaced with a deterministic alternative
 * at compile time if required for testing purposes.
 */

#ifndef DREADDS_RANDOM_H
#define DREADDS_RANDOM_H
#include <random>

namespace DreadDs {
  typedef std::mt19937 rng_eng_t;
  typedef std::uniform_real_distribution<float> uniform_real_distr_t;
  typedef std::normal_distribution<float> normal_distr_t;
  typedef std::uniform_int_distribution<> uniform_int_distr_t;
}

#endif
