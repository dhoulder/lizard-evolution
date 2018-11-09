/**
 * R interface for DreadDs model
 */

#include <Rcpp.h>

#include <simulation.h>

using namespace Rcpp;

// [[Rcpp::export]]

int dreadds(Rcpp::StringVector config_path_vec,
	    int n_steps,
	    Rcpp::StringVector output_path_vec) {

  if (config_path_vec.size() < 1 ||
      output_path_vec.size() <1)
    return -1;
  DreadDs::Simulation sim(config_path_vec[0],
			  output_path_vec[0]);
  return sim.run(n_steps);
}
