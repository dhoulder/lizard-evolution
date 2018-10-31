
#include <Rcpp.h>

#include <simulation.h>

using namespace Rcpp;

// [[Rcpp::export]]
List dreadds() {
  // FIXME STUB

  DreadDs::Simulation sim("example.conf",
			  "model-test-out.junk");

  int final_step = sim.run(7);

  CharacterVector x = CharacterVector::create("dreadds STUB", sim.version);
  List z = List::create(final_step) ;

  return z ;
}
