#include <iostream>
#include "simulation.h"
#include <cassert>

// FIXME STUB

int main(int argc, char *argv[]) {
  DreadDs::filename_vec env_inputs(1, "JUNK");
  DreadDs::filename_vec species_inputs(1, "JUNK");
  assert(argc > 1);

  DreadDs::Simulation sim(argv[1], env_inputs, "model-test-out.junk");

  int s = sim.run(4);
  
  std::cout << s << " done\n";
}
