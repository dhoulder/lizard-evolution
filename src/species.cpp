// -*- coding: utf-8 -*-

#include <cmath>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <atomic>

#include <boost/serialization/array_wrapper.hpp> // just for boost 1.64. see https://svn.boost.org/trac10/ticket/12982
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "constants.h"
#include "model-config.h"
#include "species.h"
#include "exceptions.h"

using namespace DreadDs;
namespace ba = boost::accumulators;

static std::atomic_int id_hwm(0);

Species::Species(const Config &c,
		 const SpeciesParameters &sp,
		 const Environment &env):
  conf(c),
  demes(new(DemeMap))
  /**
   * Load initial species values, locations and abundance.
   */
{
  ++id_hwm;
  id = id_hwm;

  if (conf.verbosity > 1)
    std::cout << "Creating species " << id << std::endl;

  setup_dispersal(sp);
  load_initial(sp, env);

  update_stats(initial_stats,
	       // use step=0 to indicate before first time step
	       0);
  latest_stats = initial_stats;
}


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
  if (conf.verbosity > 2) {
    std::cout << " Dispersal kernel" << std::endl;
    for (auto &&v: dk)
      std::cout << " x=" << v.x << ", y=" << v.y <<
	" weight=" << v.weight << std::endl;
  }
}


void Species::load_initial(const SpeciesParameters &sp,
			   const Environment &env) {

  // get initail bounding rectangle row and column limits.

  // sp.north…west have alreday been bounds checked in sp.get_initial_species()
  long n = env.row(sp.north);
  long s = env.row(sp.south);
  long e = env.col(sp.east);
  long w = env.col(sp.west);

  if (conf.verbosity > 1)
    std::cout << " Loading species. "
      "west=" << sp.west <<
      ", east=" << sp.east <<
      ", south=" << sp.south <<
      ", north=" << sp.north  << std::endl;

  // fill all suitable cells in initial-rectangle

  Deme d(sp);
  Location loc;
  EnvCell ec;
  bool no_data;
  for (loc.y = n; loc.y <= s; ++loc.y)
    for (loc.x = w; loc.x <= e; ++loc.x) {
      env.get(loc, ec, &no_data);
      if (no_data)
        continue;
      d.amount = d.niche_suitability(conf, ec);
      if (d.amount > 0.0f)
        (*demes)[loc] = DemeList(1, d);
    }
}


template <class ...Types>
using Acc = ba::accumulator_set<float,
			    ba::stats<Types...>>;

int Species::update_stats(Characteristics &ch, int current_step) {
  /**
   * Update species statistics.
   * Returns: cell count
   */

  if (extinction > -1)
    // went extinct in earlier time step. Leave the stats alone
    return 0;

  Acc <ba::tag::mean, ba::tag::variance>
    niche_pos_acc[conf.env_dims],
    niche_tol_acc[conf.env_dims],
    genetic_acc[conf.genetic_dims];
  Acc <ba::tag::min> niche_min[conf.env_dims];
  Acc <ba::tag::max> niche_max[conf.env_dims];

  ch.cell_count = 0;
  ch.population = 0.0f;

  // FIXME do this in merge() instead of separate loops here?
  for (auto &&kv: *demes)
    for (auto &&d: kv.second) {
      if (d.amount <= 0.0f)
	continue; // FIXME flag these cases???
      ++ch.cell_count;
      ch.population += d.amount;

      for (int i=0; i < conf.env_dims; i++) {
	float nc = d.genetics.niche_centre[i];
	float nt = d.genetics.niche_tolerance[i];
	niche_min[i]( nc - nt);
	niche_max[i](nc + nt);
	niche_pos_acc[i](nc);
	niche_tol_acc[i](nt);
      }
      for (int i=0; i < conf.genetic_dims; i++)
	genetic_acc[i](d.genetics.genetic_position[i]);
    }

  if (ch.cell_count == 0) {
    extinction = current_step;
    // leave other stats as they were on previous step
    if (conf.verbosity > 1)
      std::cout << "Species " << id << " extinct at step " <<
	current_step << std::endl;
    return 0;
  }

  for (int i=0; i < conf.env_dims; i++) {
    auto &ns = ch.niche_stats[i];
    ns.min = ba::min(niche_min[i]);
    ns.max = ba::max(niche_max[i]);
    ns.position_mean = ba::mean(niche_pos_acc[i]);
    ns.position_sd = sqrt(ba::variance(niche_pos_acc[i]));
    // breadth is 2 * tolerance
    ns.breadth_mean = 2.0 * ba::mean(niche_tol_acc[i]);
    ns.breadth_sd = 2.0 * sqrt(ba::variance(niche_tol_acc[i]));
  }

  for (int i=0; i < conf.genetic_dims; i++) {
    auto &gs = ch.genetic_stats[i];
    gs.mean = ba::mean(genetic_acc[i]);
    gs.sd = sqrt(ba::variance(genetic_acc[i]));
  }

  if (conf.verbosity > 1) {
    std::cout <<
      "Species " << id << " cell count " << ch.cell_count <<
      ", population = " << ch.population <<
      std::endl;
    for (int i=0; i < conf.env_dims; ++i) {
      const auto &ns = ch.niche_stats[i];
      std::cout << "Niche " << i <<
	" min=" << ns.min <<
	", max=" << ns.max <<
	", mean=" << ns.position_mean <<
	", sd=" << ns.position_sd <<
	"\n  breadth: mean=" << ns.breadth_mean <<
	" sd=" << ns.breadth_sd <<
	std::endl;
    }
  }
  return ch.cell_count;
}

void Species::speciate() {
  /**
   * Examine species for distinct clumps of demes in genetic space. If
   * only one is found, leave the species alone, otherwise split the
   * current species into distinct sub-species.
   */

  // FIXME STUB

}
