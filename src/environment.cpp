// -*- coding: utf-8 -*-

/**
 * Read environment values from grid files using GDAL See
 * https://www.gdal.org
 * Each grid has to have the same dimensions, resolution and location.
 */


#include <string>
#include <iostream>
#include <cassert>
#include <cmath>

#define BOOST_DISABLE_ASSERTS

#include "boost/multi_array.hpp"
#include <boost/math/special_functions/next.hpp>

// GDAL
#include "gdal_priv.h"
#include "cpl_conv.h"

#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "exceptions.h"

static bool gdal_reg = false;

using namespace DreadDs;
namespace bm = boost::math;

EnvReader:: EnvReader(const std::string &filename):
  grid_filename(filename)
{
  if (!gdal_reg) {
    GDALAllRegister();
    gdal_reg = true;
  }
  dataset = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
  if (dataset == NULL)
    throw ConfigError("Cannot read environment input file " + filename);

  band = dataset->GetRasterBand(1);
  assert(band);
  ncol = band->GetXSize();
  nrow = band->GetYSize();
  row_buffer = (float *) CPLMalloc(sizeof(float)*ncol);
  assert(row_buffer);
  int has_nodata_value;
  double ndv = band->GetNoDataValue(&has_nodata_value);
  if (has_nodata_value)
    nodata_value = (float)ndv; // otherwise leave it is NAN
}

void EnvReader::get_coordinates(double *buff6) {
  if (dataset->GetGeoTransform(buff6) != CE_None)
    throw ConfigError("Cannot determine coordinates of cells in grid file");
}

void EnvReader::read_row(int row) {
  if (band->RasterIO(GF_Read,
		     0, row, ncol, 1,
		     row_buffer, ncol, 1,
		     GDT_Float32,
		     0, 0 ) != CE_None)
    throw ConfigError("Error reading row data from grid file");
}

bool EnvReader::is_nodata(float v) {
  return (!std::isnan(nodata_value)) &&
    (fabs(bm::float_distance(v, nodata_value)) <= 2.0f);
}

EnvReader::~EnvReader() {
  if (dataset)
    GDALClose(dataset);
  if (row_buffer)
    CPLFree(row_buffer);
}


Environment::Environment(const Config &conf_):
  current_step_offset(0),
  conf(conf_)
{
  int layer = 0;
  for (const auto &ep: conf.env_params) {
    std::string &&grid_filename = ep->get_filename(0);

    EnvReader env_reader(grid_filename);
    if (0 == layer) {
      // First file defines bounding box, allocates space and fills in first layer
      values.resize(boost::extents
		    [env_reader.nrow]
		    [env_reader.ncol]
		    [conf.env_params.size()]);
      env_reader.get_coordinates(geo_transform);
    } else {
      // subsequent files fill in 2nd, 3rd, … layers
      assert(layer < max_env_dims);
      check_coordinates(env_reader);
    }
    load(env_reader, layer);
    ++layer;
  }
  update(0); // initial species determination needs this
}

void Environment::check_coordinates(EnvReader &er) {
  // assumes values and geo_transform have been set
  if (er.nrow != values.shape()[0] || er.ncol != values.shape()[1])
    throw ConfigError("Dimensions of " + er.grid_filename +
		      "don't match those of first grid file");
  double gt[6];
  er.get_coordinates(gt);
  for (int j=0; j<6; j++)
    if (fabs(bm::float_distance(gt[j], geo_transform[j])) > 2.0)
      throw ConfigError("Coordinates of " + er.grid_filename +
			" don't match those of first grid file");
}

void Environment::load(EnvReader &er, int layer) {
  if (conf.verbosity > 1)
    std::cout << "Loading environment variable " << layer <<
      " from " << er.grid_filename << std::endl;
  for (int row=0; row < er.nrow; ++row) {
    er.read_row(row);
    for (int col=0; col < er.ncol; ++col) {
      const float &v = er.row_buffer[col];
      values[row][col][layer] = er.is_nodata(v)? NAN : v;
      if  (conf.verbosity > 2)
	std::cout <<
	  " row=" << row << " col=" << col <<
	  " value=" <<  values[row][col][layer] << std::endl;
    }
  }
}


void Environment::update(int step_offset) {
  /**
   * For extrapolated environment layers, set env_delta for the
   * current time step. For time-series layers, load the corresponding
   * files
   */
  float *ed = env_delta;
  auto p = conf.env_params.begin();
  for (int i = 0; i < conf.env_dims;
       ++i, ++p, ++ed)
    *ed = (*p)->update_environment(this, step_offset, i);
  current_step_offset = step_offset;
}
