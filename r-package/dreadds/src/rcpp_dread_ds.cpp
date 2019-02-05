/**
 * R interface for DreadDs model
 */

#include <Rcpp.h>

#include <model-config.h>
#include <model.h>

using namespace Rcpp;

// [[Rcpp::export]]

int dreadds(Rcpp::StringVector config_path_vec) {
  // FIXME allow config options to be passed in as args

  if (config_path_vec.size() < 1)
    return -1;

  const char *argv[] = {"dreadds", "-c", config_path_vec[0]};
  DreadDs::Config conf(3, argv);
  DreadDs::Model model(conf);

  for (int i = 0; i < conf.n_iterations; ++i) {
    if (model.do_step() == 0)
      break;
  }
  return model.step;
}
