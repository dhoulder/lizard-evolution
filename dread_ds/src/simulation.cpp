// -*- coding: utf-8 -*-

/**
 * Implementation of DREaD_ds model API. See ./README.md
 */

#include "model.h"
#include "simulation.h"

using DreadDs::Simulation;
using DreadDs::Model;

Simulation::Simulation(const char *config_path,
		       const char *output_path):
  model(new Model(config_path,
		  output_path))
{

  {
    // FIXME STUB
    if (model->conf.debug > 3) {
      model->tips.back()->print_kernel();
    }
  }

}

Simulation::~Simulation() = default;

int Simulation::run(int n_steps) {
  int n = 0;
  for (int i =0; i < n_steps; ++i) {
    // FIXME STUB. needs some way to stop when done
    n = model->do_step();
  }
  return n;
}
