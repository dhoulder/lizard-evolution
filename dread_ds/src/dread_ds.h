// -*-c++-*-

/**
   DREaD_ds model API
 */

#ifndef DREADDS_H
#define DREADDS_H

#include <memory>


// Represents a DREaD_ds evolution model. See ./README.md
class DreadDs  {
public:
  const char *const version = "0.01";

  // Configure and initialise the simulation.
  DreadDs(int cols, int rows);
  ~DreadDs();

  // Run the simulation for a number of time steps.
  //
  // Arguments: n_steps: number of time steps to run for.
  //
  // Returns: time step step of final iteration.
  int run(int n_steps);

private:
  int step();
  class Impl;
  std::unique_ptr<Impl> impl;
  int current_step;
};

#endif
