// -*- coding: utf-8 -*-

#include <fstream>
#include <string>
#include <iostream> // FIXME debugging

#include "boost/multi_array.hpp"

#include "model-limits.h"
#include "model-config.h"
#include "environment.h"

// GDAL
#include "gdal_priv.h"
#include "cpl_conv.h"


static bool gdal_reg = false;

using namespace DreadDs;

class Reader {
public:

  GDALDataset *poDataset;

  Reader(const std::string &filename) {
    if (!gdal_reg) {
      GDALAllRegister();
      gdal_reg = true;
    }
    poDataset = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
    if (poDataset == NULL)
      // ~GDALDataset() should clean up everyhting
      throw ConfigError(std::string("Cannot read environment input file ") + filename);
  }


  ~Reader() {
    if (poDataset)
      GDALClose(poDataset); // FIXME necessary?
  }


};


Environment::Environment(const EnvParamsVec &env_inputs) {

  bool first = true;
  for (const auto &ep: env_inputs) {
    // FIXME WIP load cells from file

    // first file defines bb. subsequent files must match
    // first file creates env and v[0], subsequent files fill in v[1], v[2], ...


    Reader er(ep.grid_filename);


    if (first) {
      first = false;

      values.resize(boost::extents
		    [er.poDataset->GetRasterYSize()]
		    [er.poDataset->GetRasterXSize()]);

      std::cout << "Env dimensions rows=" << values.shape()[0] <<", cols=" << values.shape()[1] << std::endl; // FIXME

      if (er.poDataset->GetGeoTransform(adfGeoTransform) != CE_None)
	throw ConfigError("Cannot determine coordinates of cells in " + ep.grid_filename);



    } else {
      // FIXME WIP
    }


    // FIXME WIP see https://gis.stackexchange.com/questions/186190/reading-1-3-arcsecond-dem-data-using-c

  }

  // FIXME WIP fail if empty (first still true)

}
