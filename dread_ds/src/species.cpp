// -*- coding: utf-8 -*-

#include "model.h"

using namespace DreadDs;

static float dispersal_distance(int x, int y) {
  // Replace with some other distance metric if required.
  return sqrt(x*x + y*y);
}

void Species::setup_dispersal() {
  // Calculate dispersal kernel
  // TODO only really have to store one quadrant (perhaps only one octant)
  int r = (int) (ceil(max_dispersal_radius) + 0.5);
  for (int i = -r; i <= r; ++i)
    for (int j = -r; j <= r; ++j) {
      if (!(i || j))
	// dispersal into same-cell is a special case
	continue;
      float dd = dispersal_distance(i, j);
      if (dd <= max_dispersal_radius) {
	dk.push_back(DispersalWeight {
	    i, j,
	      1.0f - ((1.0f - dispersal_min) * dd/max_dispersal_radius)});
      }
    }
}


Species::Species(const char *filename, EnvMatrix *env):
  demes(new(DemeMap))
{

  max_dispersal_radius = 1.0f;
  dispersal_min = 0.2f;

  setup_dispersal();

  // FIXME STUB load from file (initial species and  locations)

}
