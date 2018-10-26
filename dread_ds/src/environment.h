// -*- coding: utf-8; Mode: c++  -*-

#ifndef DREADDS_ENVIRONMENT_H
#define DREADDS_ENVIRONMENT_H

#include "boost/multi_array.hpp"

#include "model-limits.h"
#include "model-config.h"


namespace DreadDs {

  struct EnvCell {
    // has to be a class|struct to keep boost::multi_array allocator happy
    float v[max_env_dims];
  };
  typedef boost::multi_array<EnvCell, 2> EnvMatrix;

  class Environment {
  public:
    EnvMatrix values;
    // See https://en.wikipedia.org/wiki/Esri_grid
    float xllcorner;
    float yllcorner;
    float cellsize;
    float nodata_value = NAN;

    Environment(const EnvParamsVec &env_inputs);
  };
}
#endif
