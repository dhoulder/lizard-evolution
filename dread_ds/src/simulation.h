// -*-c++-*-

/**
   DREaD_ds model API
 */

#ifndef DREADDS_SIMULATION_H
#define DREADDS_SIMULATION_H

#include <memory>

namespace DreadDs {

  // Represents a DREaD_ds evolution model. See ./README.md
  class Simulation  {
  public:
    const char *const version = "0.01";

    // Configure and initialise the simulation.
    Simulation(int cols, int rows);
    ~Simulation();

    // Run the simulation for a number of time steps.
    //
    // Arguments: n_steps: number of time steps to run for.
    //
    // Returns: time step step of final iteration.
    int run(int n_steps);

  private:
    int step();
    class Model;
    std::unique_ptr<Model> model;
    int current_step;
  };
}

#endif
