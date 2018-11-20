// -*- coding: utf-8 -*-

/**
 * Implementation of DREaD_ds model API. See ./README.md
 */

#include "model-config.h"
#include "model.h"
#include "simulation.h"

using DreadDs::Simulation;
using DreadDs::Model;
using DreadDs::Config;

Simulation::Simulation(const Config &conf):
  model(new Model(conf))
{}

Simulation::~Simulation() = default;

int Simulation::run(int n_steps) {
  for (int i =0; i < n_steps; ++i) {
    if (model->do_step() ==0)
      break;
  }
  return model->step;
}
