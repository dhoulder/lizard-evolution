
#include <Rcpp.h>

#include <simulation.h>

using namespace Rcpp;

// [[Rcpp::export]]
List dreadds() {
  // FIXME STUB

  DreadDs::filename_vec env_inputs(1, "JUNK");
  DreadDs::filename_vec species_inputs(1, "JUNK");

  DreadDs::Simulation sim("config.yml",
			    env_inputs,  species_inputs,
			    "model-test-out.junk");

  int final_step = sim.run(7);

  CharacterVector x = CharacterVector::create("dreadds STUB", sim.version);
  List z            = List::create(final_step) ;

  return z ;
}
