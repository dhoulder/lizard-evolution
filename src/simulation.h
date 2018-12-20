// -*- coding: utf-8; Mode: c++ -*-

/**
 * DREaD_ds model API. See ./README.md
 */

#ifndef DREADDS_SIMULATION_H
#define DREADDS_SIMULATION_H

#include <memory>
#include <vector>

namespace DreadDs {

  // for opaque pointers|references below
  class Model;
  class Config;

  // Represents a DREaD_ds evolution model. See ./README.md
  class Simulation {
  public:

    /**
     * Configure and initialise the simulation.
     */
    Simulation(const Config &conf);
    ~Simulation();

    /**
     * Run the simulation for a number of time steps.  Arguments:
     * n_steps: number of time steps to run for. Will stop prematurely
     * if all cells become empty.
     *
     * Returns: time step of final iteration.
     */
    int run(int n_steps);

  private:
    std::unique_ptr<Model> model;
  };
}

#endif
