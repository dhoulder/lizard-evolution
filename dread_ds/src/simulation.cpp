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
{}

Simulation::~Simulation() = default;

int Simulation::run(int n_steps) {
  for (int i =0; i < n_steps; ++i) {
    if (model->do_step() ==0)
      break;
  }
  return model->step;
}
