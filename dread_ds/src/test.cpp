#include <iostream>
#include "simulation.h"

// FIXME STUB

int main(int argc, char *argv[]) {
  DreadDs::filename_vec env_inputs(1, "JUNK");
  DreadDs::filename_vec species_inputs(1, "JUNK");

  DreadDs::Simulation sim("config.yml", env_inputs, species_inputs, "model-test-out.junk");

  int s = sim.run(4);
  
  std::cout << s << " done\n";
}
