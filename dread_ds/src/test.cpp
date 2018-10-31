#include <iostream>
#include "simulation.h"
#include <cassert>

// FIXME STUB

int main(int argc, char *argv[]) {
  assert(argc > 1);

  DreadDs::Simulation sim(argv[1], "model-test-out.junk");

  int s = sim.run(4);
  
  std::cout << s << " done\n";
}
