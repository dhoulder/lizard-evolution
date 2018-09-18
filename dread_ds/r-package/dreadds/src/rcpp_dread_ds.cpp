
#include <Rcpp.h>

#include <dread_ds.h>

using namespace Rcpp;

// [[Rcpp::export]]
List dreadds() {
  // FIXME STUB

  DreadDs model(3, 5);

  int final_step = model.run(7);

  CharacterVector x = CharacterVector::create("dreadds STUB", model.version);
  List z            = List::create(final_step) ;

  return z ;
}
