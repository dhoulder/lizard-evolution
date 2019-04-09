// -*- coding: utf-8 -*-

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include "species.h"
#include "model-config.h"
#include "deme.h"

using namespace DreadDs;

static void setup3(Species &s) {
  Location loc;
  Deme d;
  for (int i=0; i <100; ++i) {
    loc.x = loc.y = i;
    for (int j=0; j < s.conf.genetic_dims; j++)
      d.genetics.genetic_position[j] = (float)(i %3);
    (*s.demes)[loc].push_front(d);
  }
}

int main(int argc, char *argv[]) {
  Config c;
  c.gene_flow_zero_distance = 1.0f;
  c.genetic_dims = 3;
  Location loc;
  Deme d;
  {
    int species_id = 1;
    Species s(c, species_id);
    for (int i=0; i <100; ++i) {
      loc.x = loc.y = i;
      for (int j=0; j < c.genetic_dims; j++)
        d.genetics.genetic_position[j] = (float)i;
      (*s.demes)[loc].push_front(d);
    }
    s.speciate(&species_id);

    if (s.demes->size() != 0 || s.sub_species.size() != 100) {
      printf("100-cluster test failed: %lu %lu\n",
             s.demes->size(), s.sub_species.size());
      exit(1);
    }
  }

  {
    // use genetic distance threshold of 0.9 * diagonaal of 1x1x1
    int species_id = 1;
    c.gene_flow_zero_distance = 0.9 * sqrt(3.0f);
    Species s(c, species_id);
    setup3(s);
    s.step = 123;
    s.speciate(&species_id);
    if (s.demes->size() != 0 || s.sub_species.size() != 3) {
      printf("three-cluster test failed: %lu %lu\n",
             s.demes->size(), s.sub_species.size());
      exit(2);
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
    int species_id = 1;
    c.gene_flow_zero_distance =  0.9 * sqrt(3.0f);
    Species s(c, species_id);
    setup3(s);
    // add a couple of bridging demes to keep as one cluster
    float x = 0.5f;
    for (int i=100; i <102; x += 1.0f, ++i) {
      loc.x = loc.y = i;
      for (int j=0; j < c.genetic_dims; j++)
        d.genetics.genetic_position[j] = x;
      (*s.demes)[loc].push_front(d);
    }
    s.speciate(&species_id);
    if (s.demes->size() != 102 || s.sub_species.size() != 0) {
      printf("one-cluster test failed: %lu %lu\n",
             s.demes->size(), s.sub_species.size());
      exit(4);
    }
  }
  printf("OK\n");
  exit(0);
}
