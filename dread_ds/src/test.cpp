#include <string>
#include <stdexcept>
#include <iostream>
#include <cstdlib>

#include "exceptions.h"
#include "simulation.h"

static void usage(char *argv0) {
  std::cerr << "Usage: " << argv0 << "config-file n-iterations output-path" << std::endl;
  std::exit(1);
}

int main(int argc, char *argv[]) {
  int n =0;
  if (argc != 4)
    usage(argv[0]);
  try {
    n = std::stoi(argv[2]);
  }
  catch (const std::logic_error& oor) {
    usage(argv[0]);
  }

  try {
    DreadDs::Simulation sim(argv[1], argv[3]);
    int s = sim.run(n);
    std::cout << "Final step count: " << s << std::endl;
  }
  catch (const DreadDs::ApplicationError &ae) {
    std::cerr << ae.what()  << std::endl;
    std::exit(1);
  }
}
