// -*- coding: utf-8 -*-

#include <cmath>
#include <algorithm>
#include <iostream> // FIXME debugging

#include "model-limits.h"
#include "model-config.h"
#include "species.h"

using namespace DreadDs;

static float dispersal_distance(int x, int y) {
  // Replace with some other distance metric if required.
  return sqrt(x*x + y*y);
}

void Species::setup_dispersal(float dispersal_min, const SpeciesParameters &sp) {
  // Calculate dispersal kernel (TODO? we could probably store just
  // one quadrant (perhaps only one octant). Probably little gain
  // though)
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


void Species::load_initial(const Config &conf,
			   const SpeciesParameters &sp,
			   const Environment &env) {

  auto &&env_shape = env.values.shape();

  // FIXME check overflow, zero cases

  std::cout << "species cell bounds (NSEW) " << // FIXME debugging
    env.row(sp.north) << " " << sp.north << ", " <<
    env.row(sp.south) <<  " " << sp.south << ", " <<
    env.col(sp.east) << " " << sp.east << ", " <<
    env.col(sp.west) << " " << sp.west << std::endl;

  long n = std::max(0L, env.row(sp.north));
  long s = std::min((long)(env_shape[0] - 1),
		    env.row(sp.south));
  long e = std::min((long)(env_shape[1] -1),
		    env.col(sp.east));
  long w = std::max(0L, env.col(sp.west));
  // FIXME enforce n <s, e < w

  Deme d(sp);
  Location loc;
  std::cout << "Loading species " << n << " " << s <<  " " << e << " " << w << std::endl;
  for (loc.y = n; loc.y <= s; ++loc.y)
    for (loc.x = w; loc.x <= e; ++loc.x) {
      d.amount = d.genetics.niche_suitability(conf,
					      env.get(loc));
      if (d.amount > 0.0f)
	(*demes)[loc] = DemeList(1, d);
    }
}

void Species::Characteristics::update(const Config &conf, const DemeMap &demes) {
  range.cell_count = demes.size(); // FIXME assumes no demes with amount == 0

  // FIXME WIP STUB

  std::cout << "Species cell count " << range.cell_count << std::endl;
}



Species::Species(const Config &conf,
		 const SpeciesParameters &sp,
		 const Environment &env):
  demes(new(DemeMap))
  /**
   * Load initial species values, locations and abundance.
   */
{
  setup_dispersal(conf.dispersal_min, sp);
  load_initial(conf, sp, env);

  initial_stats.update(conf, *demes);
  latest_stats = initial_stats;
}
