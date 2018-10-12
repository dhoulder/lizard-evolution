#include <iostream>
#include "dread_ds.h"

// FIXME STUB

int main(int argc, char *argv[]) {
  DreadDs::Simulation d(5, 7);
  int s = d.run(4);
  
  std::cout << s << " done\n";
}
