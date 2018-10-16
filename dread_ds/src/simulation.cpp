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

Simulation::Simulation(const char *config_file):
  current_step(0) {

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


int Simulation::step() {
  model->update_environment(current_step);

  for (auto && species: model->tips) {

    auto target = model->disperse(species);

    // TODO handle range contraction (extinction) here ???? see "3.3.2 Range contraction"

    // merge demes in each cell that are within genetic tolerance.
    model->merge(target);

    // TODO? can do this row-lagged in the dispersal loop, providing we
    // stay far enough behind the dispersal area
    // TODO?? do this in another thread?

    // Finished with source demes now - replace with target;
    species.demes = target;

    // TODO competition/co-occurrence. (see 3.5)

    // TODO speciate.
    // Make sure iterator doesn't see this.

  }

  ++current_step;
}


int Simulation::run(int n_steps) {
  for (int final = current_step + n_steps;
       current_step < final; ) {
    // FIXME STUB. needs some way to stop when done
    step();
  }
  return current_step;
}
