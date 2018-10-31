// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_ENVIRONMENT_H
#define DREADDS_ENVIRONMENT_H

#include "boost/multi_array.hpp"

#include "model-limits.h"
#include "model-config.h"


namespace DreadDs {

  struct EnvCell {
    // array wrapped in struct so we we can return it by value
    float v[max_env_dims];
  };
  typedef boost::multi_array<float, 3> EnvMatrix; // [rows][columns][layers]

  class Environment {
  public:
    EnvMatrix values;
    double geo_transform[6]; // See https://www.gdal.org/classGDALDataset.html#a5101119705f5fa2bc1344ab26f66fd1d

    Environment(const EnvParamsVec &env_inputs);
  };
}
#endif
