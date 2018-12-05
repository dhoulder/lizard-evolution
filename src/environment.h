// -*- coding: utf-8; Mode: c++ -*-

#ifndef DREADDS_ENVIRONMENT_H
#define DREADDS_ENVIRONMENT_H

#include <string>

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"
// GDAL
#include "gdal_priv.h"

#include "constants.h"
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


  class EnvReader {
    /**
     * Read raster input data. Understands many common GIS formats.
     */
  public:
    int ncol = 0;
    int nrow = 0;
    float *row_buffer = NULL;
    const std::string grid_filename;

    EnvReader(const std::string &filename);
    // have to prohibit copying due to simple management of GDAl stuff
    // in ~EnvReader()
    EnvReader(const EnvReader &er) = delete;
    EnvReader& operator=(const EnvReader &er) = delete;

    void get_coordinates(double *buff6);
    void read_row(int row);
    ~EnvReader();

  private:
    GDALDataset *dataset = NULL;
    GDALRasterBand *band = NULL;
  };

  class Environment {
  public:
    const Config &conf;
    EnvMatrix values;
    float env_delta[max_env_dims] = {0.0f};
    int current_step_offset = 0;

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

    Environment(const Config &conf);
    void update(int step_offset); // step_offset is 0-based time step

    inline long row(double ns) const {
      return  long((ns - geo_transform[3]) / geo_transform[5]);
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

    void check_coordinates(EnvReader &er);
    void load(EnvReader &er, int layer);
  };
}
#endif
