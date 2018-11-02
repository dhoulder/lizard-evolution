// -*- coding: utf-8 -*-

#include <cmath>
#include <algorithm>
#include <iostream> // FIXME debugging
#include <cassert>

#include "model-limits.h"
#include "model-config.h"
#include "species.h"

using namespace DreadDs;

static float dispersal_distance(int x, int y) {
  // Replace with some other distance metric if required.
  return sqrt(x*x + y*y);
}

void Species::setup_dispersal(const SpeciesParameters &sp) {
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
	      1.0f - ((1.0f - conf.dispersal_min) * dd/sp.max_dispersal_radius)});
      }
    }
}


static inline unsigned long clip(long v, unsigned long hi) {
  assert(hi);
  return std::min(hi-1, (unsigned long) std::max(0L, v));
}

void Species::load_initial(const SpeciesParameters &sp,
			   const Environment &env) {

  auto &&env_shape = env.values.shape();
  long n = clip(env.row(sp.north), env_shape[0]);
  long s = clip(env.row(sp.south), env_shape[0]);
  long e = clip(env.col(sp.east), env_shape[1]);
  long w = clip(env.col(sp.west),  env_shape[1]);
  assert(n <= s && w <= e);

  Deme d(sp);
  Location loc;

  if (conf.debug > 3)
    std::cout << "Loading species. north=" <<
      sp.north << " row " << n <<
      ", south=" << sp.south << " row " << s <<
      ", west=" << sp.west << " col " << w <<
      ", east=" << sp.east << " col " << e <<
      std::endl;

  for (loc.y = n; loc.y <= s; ++loc.y)
    for (loc.x = w; loc.x <= e; ++loc.x) {
      d.amount = d.niche_suitability(conf,
				     env.get(loc));
      if (d.amount > 0.0f)
	(*demes)[loc] = DemeList(1, d);
    }
}

void Species::update_stats(Characteristics &ch) {
  ch.range.cell_count = demes->size(); // FIXME assumes no demes with amount == 0

  // FIXME WIP STUB

  std::cout << "Species cell count " << ch.range.cell_count << std::endl;
}

Species::Species(const Config &conf_,
		 const SpeciesParameters &sp,
		 const Environment &env):
  conf(conf_),
  demes(new(DemeMap))
  /**
   * Load initial species values, locations and abundance.
   */
{
  setup_dispersal(sp);
  load_initial(sp, env);

  update_stats(initial_stats);
  latest_stats = initial_stats;
}
