// -*- coding: utf-8; Mode: c++ -*-

/**
 * Convenience types for Boost random stuff.
 *
 * Example usage:
 *   unsigned int seed = something;
 *   rng_eng_t rng(seed);
 *   uniform_distr_t distr(from_low, to_high);
 *
 *   float my_random_number = distr(rng);
 *   // or
 *   uniform_vg_t vg(rng, distr);
 *   float my_random_number = vg();
 */

#ifndef DREADDS_RANDOM_H
#define DREADDS_RANDOM_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace DreadDs {

  typedef boost::random::mt19937 rng_eng_t;

  typedef boost::random::uniform_real_distribution<float> uniform_real_distr_t;
  typedef boost::random::uniform_int_distribution<int> uniform_int_distr_t;
  typedef boost::random::variate_generator<rng_eng_t&, uniform_real_distr_t>  uniform_vg_t;
  typedef boost::random::normal_distribution<float> normal_distr_t;
  typedef boost::random::variate_generator<rng_eng_t&, normal_distr_t> normal_vg_t;
}

#endif
