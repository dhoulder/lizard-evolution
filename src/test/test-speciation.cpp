// -*- coding: utf-8 -*-

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <functional>
#include <vector>
#include <utility>
#include <algorithm>

#include "species.h"
#include "model-config.h"
#include "deme.h"

using namespace DreadDs;

typedef std::vector<SpeciesPresenceList> SplVec;

static void setup_candidates(SpeciationCandidates &candidates,
			     SplVec &ret,
			     int from, int to,
			     Species &s,
			     std::function<float(int)> position_callback) {
  // Populate speciation_candidates with distinct clumps of demes
  for (int i=from; i < to; ++i) {
    SpeciesPresenceList &spl = ret[i];
    Deme d;
    d.amount = 1.0;
    for (int j=0; j < s.conf.genetic_dims; j++)
      d.genetics.genetic_position[j] = position_callback(i);
    spl.emplace_front(&s, d);
    SpeciesPresence &sp = spl.front();
    candidates.add(spl, sp);
  }
}

static bool compare_spl(const SpeciesPresenceList &a, const SpeciesPresenceList &b) {
  return a.front().species->id < b.front().species->id;
}

static SpeciesPresence &get_sp(SpeciesPresenceList spl, const char *msg) {
  auto itr = spl.cbegin();
  ++itr;
  if (itr !=  spl.cend()) {
    printf("%s List too long.\n", msg);
    exit(1);
  }
  return spl.front();
}

int main(int argc, char *argv[]) {
  Config c;
  c.gene_flow_zero_distance = 1.0f;
  c.genetic_dims = 3;

  {
    // Test speciation into 100 species

    int species_id = 1;
    Species s(c, species_id);
    SplVec splv(100);
    SpeciationCandidates candidates;
    // 100 demes, each one should become a separate species
    setup_candidates(candidates, splv,
		     0, splv.size(),
		     s,
		     [] (int i) { return (float)i;});
    s.speciate(&species_id, candidates);

    if (s.sub_species.size() != 100) {
      printf("100-cluster test failed. Wrong species count: %lu\n",
             s.sub_species.size());
      exit(1);
    }

    // should now have species 2 to 102
    std::sort(splv.begin(), splv.end(), compare_spl);
    for (int i=0; i<splv.size(); ++i) {
      auto sp = get_sp(splv[i], "100-cluster test failed.");
      if (sp.species->id != i+2) {
	printf("100-cluster test failed. Bad species id: %i\n",
	       sp.species->id);
	exit(1);
      }
    }
  }

  {
    // Test speciation into 3 species

    // Use genetic distance threshold of 0.9 * Euclidean distance
    // between points
    int species_id = 1;
    c.gene_flow_zero_distance = 0.9 * sqrt(3.0f);
    Species s(c, species_id);
    // Create three distinct equally-sized clumps at genetic positions
    // (0, 0, 0), (1, 1, 1), (2, 2, 2)
    SplVec splv(99);
    SpeciationCandidates candidates;
    setup_candidates(candidates, splv,
		     0, splv.size(),
		     s,
		     [] (int i) { return (float)(i / 33);});
    s.step = 123;
    s.speciate(&species_id, candidates);
    if (s.sub_species.size() != 3) {
      printf("three-cluster test failed: %lu\n",
             s.sub_species.size());
      exit(2);
    }

    // should now have species 2, 3, 4
    std::sort(splv.begin(), splv.end(), compare_spl);
    for (int i=0; i<splv.size(); ++i) {
      auto sp = get_sp(splv[i], "three-cluster test failed.");
      int xid =  (i / 33)  +  2;
      if (sp.species->id != xid) {
	printf("three-cluster test failed. Bad species id: %i should be %i at %i\n",
	       sp.species->id, xid, i);
	exit(1);
      }
    }
    s.phylogeny_as_yaml(stdout, "");
    std::string n = s.phylogeny_as_newick();
    if (n != "((species_2:0,species_3:0,species_4:0)species_1:123);") {
      printf("phylogeny_as_newick() failed: Got %s\n",
             n.c_str());
      exit(3);
    }
  }

  {
    // Test cohesion of single species
    int species_id = 1;
    c.gene_flow_zero_distance =  0.9 * sqrt(3.0f);
    Species s(c, species_id);
    SplVec splv(101);
    SpeciationCandidates candidates;
    // same three distinct clumps as above
    setup_candidates(candidates, splv,
		     0, 99,
		     s,
		     [] (int i) { return (float) (i / 33);});
    // add two demes to bridge all demes into one cluster
    setup_candidates(candidates, splv,
		     99, 101,
		     s,
		     [] (int i) { return (i==99)? 0.5 : 1.5;});
    s.speciate(&species_id, candidates);
    if (s.sub_species.size() != 0) {
      printf("one-cluster test failed: %lu\n",
             s.sub_species.size());
      exit(4);
    }
    // check they're all still species 1
    for (int i=0; i<splv.size(); ++i) {
      auto sp = get_sp(splv[i], "one-cluster test failed.");
      if (sp.species->id != 1) {
	printf("one-cluster test failed. Bad species id: %i\n",
	       sp.species->id);
	exit(1);
      }
    }
  }
  printf("OK\n");
  exit(0);
}
