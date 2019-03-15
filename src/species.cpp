// -*- mode: C++; coding: utf-8; eval: (c-set-offset 'arglist-intro '+) -*-

#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <atomic>
#include <vector>
#include <forward_list>
#include <memory>
#include <utility>
#include <string>
#include <sstream>

#include <boost/serialization/array_wrapper.hpp> // just for boost 1.64. see https://svn.boost.org/trac10/ticket/12982
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "deme.h"
#include "species.h"
#include "exceptions.h"

using namespace DreadDs;
namespace ba = boost::accumulators;

static std::atomic_int id_hwm(0);

Species::Species(const Config &c):
  conf(c),
  step(0), // step=0 indicates before first time step
  demes(std::make_shared<DemeMap>())
{
  ++id_hwm;
  id = id_hwm;
  if (conf.verbosity > 1)
    std::cout << "Creating species " << id << std::endl;
}

/**
 * Load initial species values, locations and abundance.
 */
Species::Species(const Config &c,
		 const SpeciesParameters &sp,
		 const Environment &env):
  Species(c)
{
  setup_dispersal(sp);
  load_initial(sp, env);
  set_initial_stats();
  if (demes->size() == 0)
    extinction = 0; // Dead on arrival
}

void Species::set_initial_stats() {
  update_stats(initial_stats);
  latest_stats = initial_stats;
}

void Species::update(const std::shared_ptr <DemeMap> &d, int s) {
  demes = d;
  step = s;
  if (d->size() == 0 && extinction < 0) {
    update_stats(latest_stats); // TODO could probably short-circuit this
    extinction = step;
    if (conf.verbosity > 1)
      std::cout << "Species " << id << " extinct at step " <<
	extinction << std::endl;
  }
}

static float dispersal_distance(double x, double y) {
  // Replace with some other distance metric if required.
  return std::sqrt(x*x + y*y);
}

void Species::setup_dispersal(const SpeciesParameters &sp) {
  // Calculate dispersal kernel

  // Implements inverse-square falloff from 1.0 at distance 0 to
  // conf.dispersal_min at sp.max_dispersal_radius
  float dist_adjust =
    (std::sqrt(1.0f / conf.dispersal_min)
     - 1.0f) /
    sp.max_dispersal_radius;

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
	      1.0f / std::pow(((dd * dist_adjust) + 1.0f), 2.0f)});
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

  // get initial bounding rectangle row and column limits.

  // sp.north…west have already been bounds checked in sp.get_initial_species()
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

/**
 * Update species statistics.
 * Returns: cell count
 */
int Species::update_stats(Characteristics &ch) {
  if (extinction > -1 || split > -1)
    // went extinct or speciated in earlier time step. Leave the stats alone
    return 0;

  // no need to recalculate everything if nothing's changed
  if (step == ch.step)
    return ch.cell_count;
  ch.step = step;

  Acc <ba::tag::mean, ba::tag::variance>
    niche_pos_acc[conf.env_dims],
    niche_tol_acc[conf.env_dims],
    genetic_acc[conf.genetic_dims];
  Acc <ba::tag::min> niche_min[conf.env_dims];
  Acc <ba::tag::max> niche_max[conf.env_dims];

  ch.cell_count = 0;
  ch.population = 0.0f;

  // TODO? do this in merge() instead of separate loops here?
  for (auto &&kv: *demes)
    for (auto &&d: kv.second) {
      if (d.amount <= 0.0f)
	continue; // FIXME flag these cases???
      ++ch.cell_count;
      ch.population += d.amount;

      for (int i=0; i < conf.env_dims; i++) {
	float nc = d.genetics.niche_centre[i];
	float nt = d.genetics.niche_tolerance[i];
	niche_min[i](nc - nt);
	niche_max[i](nc + nt);
	niche_pos_acc[i](nc);
	niche_tol_acc[i](nt);
      }
      for (int i=0; i < conf.genetic_dims; i++)
	genetic_acc[i](d.genetics.genetic_position[i]);
    }

  if (ch.cell_count == 0) {
    // leave other stats as they were on previous step TODO think about this
    return 0;
  }

  for (int i=0; i < conf.env_dims; i++) {
    auto &ns = ch.niche_stats[i];
    ns.min = ba::min(niche_min[i]);
    ns.max = ba::max(niche_max[i]);
    ns.position_mean = ba::mean(niche_pos_acc[i]);
    ns.position_sd = std::sqrt(ba::variance(niche_pos_acc[i]));
    // breadth is 2 * tolerance
    ns.breadth_mean = 2.0 * ba::mean(niche_tol_acc[i]);
    ns.breadth_sd = 2.0 * std::sqrt(ba::variance(niche_tol_acc[i]));
  }

  for (int i=0; i < conf.genetic_dims; i++) {
    auto &gs = ch.genetic_stats[i];
    gs.mean = ba::mean(genetic_acc[i]);
    gs.sd = std::sqrt(ba::variance(genetic_acc[i]));
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


/**
 * Find clusters of points in genome space. Every deme in a cluster
 * can be reached by taking steps no larger than the threshold
 * distance (single-linkage clustering).
 */
std::vector <DemeMapEntryVec> Species::get_clusters(DemeMap &dm,
						    float distance) {
  std::forward_list<DemeMapEntry *> source;
  std::vector <DemeMapEntryVec> clusters;
  float dsq = distance*distance;
  for (auto &&kv: dm)
    source.push_front(&kv);

  int n = dm.size();
  while (!source.empty()) {
    // extract a cluster from source and add it to clusters

    DemeMapEntryVec
      edge, // Holds the set of points forming the "invasion front" of
	    // the search for the current cluster.
      cl;   // Holds points from edge after we've searched for their
	    // neighbours.
    cl.reserve(n);
    edge.reserve(n);

    // Choose any starting deme
    edge.push_back(source.front());
    source.pop_front();
    do {
      // Grab any edge point and find all unclustered points within
      // range
      DemeMapEntry *d = edge.back();
      edge.pop_back();
      cl.push_back(d);

      auto p = source.begin();
      auto prev = source.end();
      // We're doing this after Model::merge(), so we'll have a single
      // entry in DemeList in each cell.
      Deme &d1 = d->second.front();
      while (p != source.end()) {
	Deme &d2 = (*p)->second.front();
	if (d1.genetic_distance_sq(conf, d2) <= dsq) {
	  // Within range.
	  edge.push_back(*p);
	  // Pull the candidate out of source
	  if (prev == source.end()) { // at the front of the list
	    source.pop_front();
	    p = source.begin();
	  } else
	    p = source.erase_after(prev);
	  // p is now pointing to the next candidate, prev is the same.
	} else {
	  // Point too far away. Advance to the next candidate
	  prev = p;
	  ++p;
	}
      }
    } while (!edge.empty());

    // Got a cluster
    n -= cl.size();
    clusters.push_back(std::move(cl));
  }

  return clusters;
}

/**
 * Create a sub species with no demes
 */
std::shared_ptr <Species> Species::add_child() {
  auto s = std::make_shared <Species>(conf);
  s->dk = dk;
  s->step = step;
  s->parent = this;
  sub_species.push_back(s);
  return s;
}


/**
 * Examine species for distinct clumps of demes in genetic space. If
 * only one is found, leave the species alone, otherwise split the
 * current species into distinct sub-species.
 */
void Species::speciate() {
  std::vector<DemeMapEntryVec> clusters = get_clusters(
    *demes,
    conf.gene_flow_zero_distance);

  if (clusters.size() < 2)
    // All one species. Nothing to do
    return;

  // The species has split into two or more sub species
  split = step;

  for (auto &&dv: clusters) {
    auto s = add_child();
    // move demes into sub species
    for (auto &&d: dv) {
      (*(s->demes))[d->first] = std::move(d->second);
    }
  }
  demes->clear();
}


void Species::as_yaml(FILE *of,
                      const std::string &first_indent) {
  std::string spaces(first_indent.size(), ' ');
  const char *ind = spaces.c_str();
  if (fprintf(
        of,
        "%sid: %d\n"
        "%scell_count: %d\n"
        "%spopulation: %f\n"
        "%sextinction_time: %d\n"
        "%sspeciation_time: %d\n"
        "%sparent_species: %d\n"
        "%sniche:\n",
        first_indent.c_str(), id,
        ind, latest_stats.cell_count,
        ind, latest_stats.population,
        ind, extinction,
        ind, split,
        ind, parent? parent->id : -1,
        ind) < 1)
    throw ApplicationException("Error writing species YAML file");

  for (int i=0; i < conf.env_dims; ++i) {
    const auto &ns = latest_stats.niche_stats[i];
    fprintf(of,
            "%s - min: %f\n"
            "%s   max: %f\n"
            "%s   mean: %f\n"
            "%s   sd: %f\n"
            "%s   breadth_mean: %f\n"
            "%s   breadth_sd: %f\n",
            ind, ns.min,
            ind, ns.max,
            ind, ns.position_mean,
            ind, ns.position_sd,
            ind, ns.breadth_mean,
            ind, ns.breadth_sd);
  }

  fprintf(of,
          "%sgenetics:\n", ind);
  for (int i=0; i < conf.genetic_dims; ++i) {
    const auto &gs = latest_stats.genetic_stats[i];
    fprintf(of,
            "%s - mean: %f\n"
            "%s   sd: %f\n",
            ind, gs.mean,
            ind, gs.sd);
  }
}

void Species::phylogeny_as_yaml(FILE *of,
                              const std::string &first_indent) {
  as_yaml(of, first_indent);
  std::string spaces(first_indent.size(), ' ');

  if (sub_species.size() < 1) {
    fprintf(of,
            "%schildren: []\n",
            spaces.c_str());
    return;
  }
  fprintf(of,
	  "%schildren:\n",
          spaces.c_str());
  std::string list_indent = spaces + " - ";
  for (auto &&c: sub_species) {
    c->phylogeny_as_yaml(of, list_indent);
  }
}

std::string Species::get_name() {
  // There's currently no explicit species naming, but to ensure
  // consistent behaviour now and into the future, use this to get
  // a species' name or label.
  return "species_" + std::to_string(id);
}

static void traverse_newick(std::ostringstream &sstr, Species &s) {
  // See https://en.wikipedia.org/wiki/Newick_format
  if (s.sub_species.size() > 0) {
    sstr << "(";
    int i = 0;
    for (auto &&child: s.sub_species) {
      if (i == 0)
        ++i;
      else sstr << ",";
      traverse_newick(sstr, *child);
    }
    sstr << ")";
  }
  sstr << s.get_name();
  if (s.parent)
    // branch length if not root
    sstr << ":" << (s.step - s.parent->step);
}

std::string Species::phylogeny_as_newick() {
  std::ostringstream sstr;
  traverse_newick(sstr, *this);
  return sstr.str() + ";";
}
