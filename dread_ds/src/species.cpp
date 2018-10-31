// -*- coding: utf-8 -*-

#include <cmath>

#include "model-limits.h"
#include "model-config.h"
#include "species.h"
//#include "model.h"

using namespace DreadDs;

static float dispersal_distance(int x, int y) {
  // Replace with some other distance metric if required.
  return sqrt(x*x + y*y);
}

void Species::setup_dispersal(float dispersal_min, const SpeciesParameters &sp) {
  // Calculate dispersal kernel
  // TODO only really have to store one quadrant (perhaps only one octant)
  int r = (int) (ceil(sp.max_dispersal_radius) + 0.5);
  for (int i = -r; i <= r; ++i)
    for (int j = -r; j <= r; ++j) {
      if (!(i || j))
	// dispersal into same-cell is a special case
	continue;
      float dd = dispersal_distance(i, j);
      if (dd <= sp.max_dispersal_radius) {
	dk.push_back(DispersalWeight {
	    i, j,
	      1.0f - ((1.0f - dispersal_min) * dd/sp.max_dispersal_radius)});
      }
    }
}


void Species::load_initial(const Config &conf, const SpeciesParameters &sp, const Environment &env) {

  // FIXME STUB load from file (initial species and locations)

}


Species::Species(const Config &conf, const SpeciesParameters &sp, const Environment &env):
  demes(new(DemeMap))

  /**
   * Load initial species values, locations and abundance.
   */

{

  setup_dispersal(conf.dispersal_min, sp);
  load_initial(conf, sp, env);

  // FIXME get niche params, bb, genetics, from sp

  // FIXME create demes according to suitability.


}
