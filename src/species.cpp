// -*- mode: C++; coding: utf-8; eval: (c-set-offset 'arglist-intro '+) -*-

#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <vector>
#include <forward_list>
#include <memory>
#include <utility>
#include <string>
#include <sstream>

#include "constants.h"
#include "model-config.h"
#include "environment.h"
#include "deme.h"
#include "species.h"
#include "exceptions.h"

using namespace DreadDs;
namespace ba = boost::accumulators;


Species::Species(const Config &c, const int species_id, const int creation_step):
  conf(c),
  step(creation_step), // step=0 indicates before first time step
  id(species_id)
{
  if (conf.verbosity > 1)
    std::cout << "Creating species " << get_id() << std::endl;
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

/**
 * Accumulate values to derive species statistics
 */

void Species::StatsAccumulator::accumulate(const SpeciesPresence &sp) {
  auto &&d =  sp.incumbent;

  ++count;
  population += d.amount;

  for (int i=0; i <  niche_min.size(); i++) {
    float nc = d.genetics.niche_centre[i];
    float nt = d.genetics.niche_tolerance[i];
    niche_min[i](nc - nt);
    niche_max[i](nc + nt);
    niche_pos_acc[i](nc);
    niche_tol_acc[i](nt);
  }
  for (int i=0; i < genetic_acc.size(); i++)
    genetic_acc[i](d.genetics.genetic_position[i]);
}


/**
 * Update species statistics.
 * Returns: cell count
 */
int Species::StatsAccumulator::update_stats(Characteristics &ch, int step) {
  ch.step = step;

  ch.cell_count = count;
  ch.population = population;

  if (ch.cell_count == 0) {
    // leave other stats as they were on previous step. TODO? or maybe reset them?
    return 0;
  }

  for (int i=0; i < niche_min.size(); i++) {
    auto &ns = ch.niche_stats[i];
    ns.min = ba::min(niche_min[i]);
    ns.max = ba::max(niche_max[i]);
    ns.position_mean = ba::mean(niche_pos_acc[i]);
    ns.position_sd = std::sqrt(ba::variance(niche_pos_acc[i]));
    // breadth is 2 * tolerance
    ns.breadth_mean = 2.0 * ba::mean(niche_tol_acc[i]);
    ns.breadth_sd = 2.0 * std::sqrt(ba::variance(niche_tol_acc[i]));
  }

  for (int i=0; i < genetic_acc.size(); i++) {
    auto &gs = ch.genetic_stats[i];
    gs.mean = ba::mean(genetic_acc[i]);
    gs.sd = std::sqrt(ba::variance(genetic_acc[i]));
  }

  return ch.cell_count;
}

void Species::log_summary_stats(const Characteristics &ch) {
  std::cout <<
    "Species " << get_id() << " cell count " << ch.cell_count <<
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

void Species::update(int s, StatsAccumulator &acc) {
  step = s;
  if (acc.count == 0  && extinction < 0) {
    // we just went extinct in this step
    update_latest_stats(acc);
    extinction = step;
    if (conf.verbosity > 1)
      std::cout << "Species " << get_id() << " extinct at step " <<
	extinction << std::endl;
  }
}

void SpeciationCandidates::add(SpeciesPresenceList &spl, SpeciesPresence &sp) {
  candidates.emplace_front(spl, sp);
  ++count;
}

/**
 * Find clusters of points in genome space. Every deme in a cluster
 * can be reached by taking steps no larger than the threshold
 * distance (single-linkage clustering).
 */
std::vector<SpeciationItemVector> SpeciationCandidates::split(float distance) {
  std::vector <SpeciationItemVector> clusters;
  float dsq = distance*distance;
  while (!candidates.empty()) {
    // extract a cluster from candidates and add it to clusters

    SpeciationItemVector
      edge, // Holds the set of points forming the "invasion front" of
	    // the search for the current cluster.
      cl;   // Holds points from edge after we've searched for their
	    // neighbours.
    cl.reserve(count);
    edge.reserve(count);

    // Choose any starting deme
    edge.push_back(std::move(candidates.front()));
    candidates.pop_front();
    do {
      // Grab any edge point and find all unclustered points within
      // range
      cl.push_back(std::move(edge.back()));
      edge.pop_back();
      SpeciationItem &si = cl.back();

      auto p = candidates.begin();
      auto prev = candidates.before_begin();
      // We're doing this after Model::merge(), so we have one deme in
      // SpeciesPresence::incumbent in each cell.
      Deme &d1 = si.second.incumbent;
      while (p != candidates.end()) {
	Deme &d2 = p->second.incumbent;
	if (d1.genetic_distance_sq(si.second.species->conf, d2) <= dsq) {
	  // Within range.
	  edge.push_back(*p);
	  // Pull the candidate out of candidates
	  p = candidates.erase_after(prev);
	  // p is now pointing to the next candidate, prev is the same.
	} else {
	  // Point too far away. Advance to the next candidate
	  prev = p;
	  ++p;
	}
      }
    } while (!edge.empty());

    // Got a cluster
    count -= cl.size();
    clusters.push_back(std::move(cl));
  }

  return clusters;
}

/**
 * Create a child species
 */
std::shared_ptr <Species> Species::add_child(const int species_id) {
  auto s = std::make_shared <Species>(conf, species_id, step);
  s->dk = dk;
  s->parent = this;
  sub_species.push_back(s);
  return s;
}


/**
 * Examine species for distinct clumps of demes in genetic space. If
 * only one is found, leave the species alone, otherwise split the
 * current species into distinct sub-species.
 */
void Species::speciate(int *id_counter, SpeciationCandidates &candidates) {
  auto &&clusters = candidates.split(conf.gene_flow_zero_distance);
  if (clusters.size() < 2)
    // All one species. Nothing to do
    return;

  // The species has split into two or more sub species
  split = step;
  for (auto &&c: clusters) {
    auto s = add_child((*id_counter)++);
    StatsAccumulator acc(conf);
    for (SpeciationItem &si: c) {
      SpeciesPresenceList &spl = si.first;
      SpeciesPresence &sp = si.second;
      // reclassify species presence as the new species
      sp.species = &(*s);
      acc.accumulate(sp); //  accumulate for initial stats

      // move it to the front because species have to be in descending
      // id order and our new species currently has the largest id.
      auto prev = spl.before_begin();
      for (auto &i: spl) {
	if (i.species == sp.species) {
	  spl.push_front(std::move(i));
	  spl.erase_after(prev);
	  break;
	}
	++prev;
      }
    }
    s->set_initial_stats(acc);
  }
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
        first_indent.c_str(), get_id(),
        ind, latest_stats.cell_count,
        ind, latest_stats.population,
        ind, extinction,
        ind, split,
        ind, parent? parent->get_id() : -1,
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
  return "species_" + std::to_string(get_id());
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
  sstr << s.get_name() << ":" <<
    (s.step - (s.parent? s.parent->step : 0));
}

std::string Species::phylogeny_as_newick() {
  std::ostringstream sstr;
  traverse_newick(sstr, *this);
  return "(" + sstr.str() + ");";
}
