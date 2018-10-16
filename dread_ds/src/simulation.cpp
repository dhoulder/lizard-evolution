// -*- coding: utf-8 -*-

/**
 *  Implementation of DREaD_ds model API. See ./README.md
 */

#include <memory>

#include "model-config.h"
#include "model.h"
#include "simulation.h"

using DreadDs::Simulation;
using DreadDs::Config;
using DreadDs::Species;
using DreadDs::Model;

static const Config default_config;

Simulation::Simulation(const char *config_file) {

  // TODO load config

  // TODO load env

  // TODO load initial species and  locations


  model =  std::unique_ptr<Model>(new Model(default_config, 9, 9)); // FIXME STUB

  // TODO initialise roots, tips

  {
    // FIXME STUB init tips
    model->tips.push_back(Species(1.0f, 0.2f));
    if (model->conf.debug > 3) {
      model->tips.back().print_kernel();
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
