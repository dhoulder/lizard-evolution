// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_ENVIRONMENT_H
#define DREADDS_ENVIRONMENT_H

#include "boost/multi_array.hpp"

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

  struct EnvCell {
    // array wrapped in struct so we we can return it by value
    float v[max_env_dims];
  };
  typedef boost::multi_array<float, 3> EnvMatrix; // [rows][columns][layers]

  class Environment {
  public:
    EnvMatrix values;
    float env_delta[max_env_dims] = {0.0f};

    double geo_transform[6];
    // https://www.gdal.org/gdal_tutorial.html

    // In the particular, but common, case of a "north up" image
    // without any rotation or shearing, the georeferencing transform
    // takes the following form :

    // geo_transform[0] /* top left x */
    // geo_transform[1] /* w-e pixel resolution */
    // geo_transform[2] /* 0 */
    // geo_transform[3] /* top left y */
    // geo_transform[4] /* 0 */
    // geo_transform[5] /* n-s pixel resolution (negative value) */

    Environment(const EnvParamsVec &env_inputs);
    void update(const Config &conf, int time_step);

    inline long row(double ns) const {
      return  long((geo_transform[3] - ns) / -geo_transform[5]);
    }

    inline long col(double ew) const {
      return long((ew - geo_transform[0]) / geo_transform[1]);
    }

    inline EnvCell get(long row, long col) const {
      EnvCell ec;
      const auto &&base_env = values[row][col];
      for (int i =0; i < values.shape()[2]; i++)
	ec.v[i] = base_env[i] + env_delta[i];
      return ec;
    }

    inline EnvCell get(const Location &loc) const {
      return get(loc.y, loc.x);
    }

  };
}
#endif
