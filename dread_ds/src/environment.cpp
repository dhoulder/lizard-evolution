// -*- coding: utf-8 -*-


/**
 * Read environment values from grid files using GDAL See
 * https://www.gdal.org
 * Each grid has to have the same dimensions, resolution and location.
 */


#include <string>

#include "boost/multi_array.hpp"
#include <boost/math/special_functions/next.hpp>

// GDAL
#include "gdal_priv.h"
#include "cpl_conv.h"

#include "model-limits.h"
#include "model-config.h"
#include "environment.h"

static bool gdal_reg = false;

using namespace DreadDs;

class Reader {
public:

  GDALDataset *dataset = NULL;
  GDALRasterBand *band = NULL;
  int ncol = 0;
  int nrow = 0;
  float *row_buffer = NULL;

  Reader(const std::string &filename) {
    if (!gdal_reg) {
      GDALAllRegister();
      gdal_reg = true;
    }
    dataset = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
    if (dataset == NULL)
      // ~GDALDataset() should clean up everything
      throw ConfigError("Cannot read environment input file " + filename);

    band = dataset->GetRasterBand(1);
    assert(band);
    ncol = band->GetXSize();
    nrow =  band->GetYSize();
    row_buffer = (float *) CPLMalloc(sizeof(float)*ncol);
    assert(row_buffer);
  }

  void get_coordinates(double *buff6) {
    if (dataset->GetGeoTransform(buff6) != CE_None)
      throw ConfigError("Cannot determine coordinates of cells in grid file");
  }

  void read_row() {
    if (band->RasterIO(GF_Read, 0, 0, ncol, 1,
		       row_buffer, ncol, 1, GDT_Float32,
		       0, 0 ) != CE_None)
      throw ConfigError("Error reading row data from grid file");
  }

  ~Reader() {
    if (dataset)
      GDALClose(dataset);
    if (row_buffer)
      CPLFree(row_buffer);
  }

};


Environment::Environment(const EnvParamsVec &env_inputs) {

  int layer = 0;
  for (const auto &ep: env_inputs) {

    Reader env_reader(ep.grid_filename);

    if (0 == layer) {
      // First file defines bounding box, allocates space and fills in v[0]
      values.resize(boost::extents[env_reader.nrow][env_reader.ncol]);
      env_reader.get_coordinates(geo_transform);
    } else {
      // subsequent files fill in v[1], v[2], ...
      assert(layer < max_env_dims);
      if (env_reader.nrow != values.shape()[0] || env_reader.ncol != values.shape()[1])
	throw ConfigError("Dimensions of " + ep.grid_filename + "don't match those of first grid file");
      double gt[6];
      env_reader.get_coordinates(gt);
      for (int j=0; j<6; j++)
	if (boost::math::float_distance(gt[j], geo_transform[j]) > 2.0)
	  throw ConfigError("Coordinates of " +  ep.grid_filename + " don't match those of first grid file");
    }

    for (int row=0; row < env_reader.nrow; ++row) {
      env_reader.read_row();
      for (int col=0; col < env_reader.ncol; ++col)
	values[row][col].v[layer] = env_reader.row_buffer[col];
    }

    ++layer;
  }
}
