// -*- coding: utf-8 -*-

#include <math.h>

#include "model-config.h"
#include "environment.h"
#include "deme.h"

using namespace DreadDs;

static float suitability(float env_value, float niche_centre, float niche_tolerance) {
  // The suitability function is cos() from -π to π, scaled to the
  // range 0…1.0. This gives 1.0 at niche_centre, and 0 and dy/dx ==
  // 0 at niche_centre ± niche_tolerance
  return
    (fabs(niche_centre - env_value) > niche_tolerance)?
    0.0f:
    (0.5f + 0.5f * cos(M_PI *
                       (env_value - niche_centre) / niche_tolerance));
}

float Deme::niche_suitability(const Config &conf, const EnvCell env) {
  // compute geometric mean of all niche suitabilities
  float v = 1.0f;
  for (int i = 0; i < conf.env_dims; ++i)
    v *= suitability(env[i],
                     genetics.niche_centre[i],
                     genetics.niche_tolerance[i]);
  return pow(v,
             1.0f / (float)conf.env_dims);
}


Deme::Genetics::Genetics(const SpeciesParameters &sp) {
  int i = 0;
  for (auto &&n: sp.niche) {
    niche_centre[i] = n.centre;
    niche_tolerance[i] = n.tolerance;
    i++;
  }
  for(int i=0; i<max_genetic_dims; i++)
    genetic_position[i] = sp.genetics[i];
}
