// -*- coding: utf-8 -*-

/**
 * Convert CSV of points to grid file.
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <algorithm>
#include <cassert>
#include <map>
#include <cstdlib>
#include <stdexcept>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

// GDAL
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_string.h"

using namespace std;
namespace ba = boost::algorithm;

typedef boost::tokenizer<boost::escaped_list_separator<char>> CsvTokenizer;

int main(int argc, const char *argv[]) {
    if (argc != 11) {
      cerr << "Usage: " << argv[0] <<
        " input.csv data_column_label "
        "xmin xmax ymin ymax ncol nrow nodata_value output_file.tif"
           << endl;
      exit(1);
    }

    ifstream in(argv[1]);
    if (!in.is_open())
      throw runtime_error(string("Cannot open ") + argv[1]);

    string data_col(argv[2]);

    // These can throw std::invalid_argument, std::out_of_range
    double xmin = stod(argv[3]);
    double xmax = stod(argv[4]);
    double ymin = stod(argv[5]);
    double ymax = stod(argv[6]);
    int nXsize = stoi(argv[7]);
    int nYsize = stoi(argv[8]);
    float nodata_value = stof(argv[9]);
    if (nXsize < 1 ||  nYsize < 1)
      throw runtime_error("nrow and ncol must be greater than 0");
    if (xmin >= xmax)
      throw runtime_error("xmax must be greater than xmin");
    if (ymin >= ymax)
      throw runtime_error("ymax must be greater than ymin");

    const char *output_file = argv[10];

    double we_res = (xmax - xmin) / (double)nXsize;
    double ns_res = (ymax - ymin) / (double)nYsize;

    map<string, int> header;
    string line;

    if (!getline(in, line))
      throw runtime_error(string("No header line in ") + argv[1]);

    CsvTokenizer tok(line);
    int ncol = 0, xcol = -1, ycol = -1;
    for (auto &s: tok) {
      string sl = ba::to_lower_copy(s);
      if (sl == "x" || sl == "lon" || sl == "long" || sl == "longitude") {
        if (xcol != -1)
          throw runtime_error("Ambiguous or repeated X column");
        xcol = ncol;
      } else if (sl == "y" || sl == "lat" || sl == "latitude") {
        if (ycol != -1)
          throw runtime_error("Ambiguous or repeated Y column");
        ycol = ncol;
      }

      header[s] = ncol;  // map column names to column index
      ncol++;
    }
    if (xcol == -1 || ycol == -1)
      throw runtime_error("Need columns for X (labelled  x, lon, long or longitude), "
                         "and Y (labelled y, lat or latitude)");

    auto hitr = header.find(data_col);
    if (hitr == header.end())
      throw runtime_error("No column labelled " + data_col);

    int dc = hitr->second;

    GDALAllRegister();
    float *raster = (float *) CPLMalloc(sizeof(float)* nXsize * nYsize);
    assert(raster);

    GDALDataset *poDstDS;

    const char *pszFormat = "GTiff";
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    assert(poDriver);
    poDstDS = poDriver->Create(output_file,
                               nXsize, nYsize, 1,
                               GDT_Float32,
                               NULL); // papszOptions	list of driver specific control parameters.
    if (!poDstDS)
      throw runtime_error(string("Could not open ") + output_file + " for writing");

    // In the particular, but common, case of a "north up" image
    // without any rotation or shearing, the georeferencing transform
    // takes the following form :
    // geo_transform[0] /* top left x */
    // geo_transform[1] /* w-e pixel resolution */
    // geo_transform[2] /* 0 */
    // geo_transform[3] /* top left y */
    // geo_transform[4] /* 0 */
    // geo_transform[5] /* n-s pixel resolution (negative value) */
    double adfGeoTransform[6] = {xmin, we_res, 0, ymax, 0, -ns_res };
    poDstDS->SetGeoTransform(adfGeoTransform);

    // Can set projection info too if required. e.g.:
    // char *pszSRS_WKT = NULL;
    // OGRSpatialReference oSRS;
    // oSRS.SetUTM( 11, TRUE );
    // oSRS.SetWellKnownGeogCS( "NAD27" );
    // oSRS.exportToWkt( &pszSRS_WKT );
    // poDstDS->SetProjection( pszSRS_WKT );
    // CPLFree( pszSRS_WKT );

    fill(raster, raster + (nXsize * nYsize), nodata_value);
    float value;
    int x, y, lc = 0;
    while (getline(in, line)) {
      lc++;
      tok.assign(line);
      int i = 0, found = 0;
      for (auto &s: tok) {
        if (xcol == i) {
          x = int((stod(s) - xmin) / we_res);
          found |= 1;
        } else if (ycol == i) {
          y = int((ymax - stod(s)) / ns_res);
          found |= 2;
        } else if (dc == i) {
          value = stod(s);
          found |= 4;
        }

        if (found == 7)
          break;
        i++;
      }
      if (found != 7) {
        cerr << "Skipping line " << lc << ". Not enough values." << endl;
        continue;
      }
      if (x < 0 || y < 0 || x >= nXsize || y >= nYsize) {
        cerr << "Skipping line " << lc << ". Point out of bounds." << endl;
        continue;
      }

      raster[(y * nXsize) + x] = value;
    }

    GDALRasterBand *poBand = poDstDS->GetRasterBand(1);
    assert(poBand);
    poBand->SetNoDataValue(nodata_value);
    int err = poBand->RasterIO(GF_Write,
                               0, 0, nXsize, nYsize,
                               raster,
                               nXsize, nYsize,
                               GDT_Float32,
                               0, 0 );
    if (err != CE_None)
      cerr << "Raster write failed" << endl;

    GDALClose((GDALDatasetH) poDstDS );
    CPLFree(raster);
}
