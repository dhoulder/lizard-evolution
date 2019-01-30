// -*- coding: utf-8 -*-

/**
 * Convert grids to file of tab separated values
 *
 */

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include <boost/math/special_functions/next.hpp>


#include "environment.h"
#include "exceptions.h"


using namespace DreadDs;

int main(int argc, const char *argv[]) {
  std::vector <std::unique_ptr<EnvReader>> erv;
  double geo_transform[6];

  try {
    for (int i=1; i < argc; ++i) {
      double gt[6];
      EnvReader er(argv[i]);
      if (i < 2)
        er.get_coordinates(geo_transform);
      else {
        er.get_coordinates(gt);
        for (int j=0; j<6; j++)
          if (fabs(boost::math::float_distance(gt[j], geo_transform[j])) > 2.0) {
            std::cerr << "Grids have different sizes or resolutions" << std::endl;
            std::exit(2);
          }
      }
      erv.push_back(std::unique_ptr<EnvReader>(new EnvReader(argv[i])));
    }
  }
  catch (const ApplicationException &ae) {
    std::cerr << ae.what()  << std::endl;
    std::exit(2);
  }

  if (erv.size() < 1) {
    std::cerr << "Usage: " << argv[0] << " grid-file [ grid-file ...] > output.tsv" << std::endl;
  std:;exit(1);
  }

  std::cout << "Long\tLat";
  for (int i=1; i < argc; ++i)
    std::cout << "\t" << argv[i];
  std::cout << std::endl;

  int nr = erv[0]->nrow;
  int nc = erv[0]->ncol;

  for (int r=0; r < nr; ++r) {
    for (auto &erp: erv)
      erp->read_row(r);
    float y = geo_transform[3] + (geo_transform[5] * (r+0.5f));
    for (int c=0; c < nc; ++c) {
      bool nd = true;
      for (auto &erp: erv)
        nd = nd && erp->is_nodata(erp->row_buffer[c]);
      if (nd)
        continue;
      float x = geo_transform[0] + (geo_transform[1] * (c+0.5f));
      std::cout << x << "\t" << y;
      for (auto &erp: erv) {
        float v = erp->row_buffer[c];
        std::cout << "\t" << (erp->is_nodata(v)? NAN : v);
      }
      std::cout << std::endl;

    }
  }

}
