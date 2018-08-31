
#include <Rcpp.h>

#include <dread_ds.h>

using namespace Rcpp;

// [[Rcpp::export]]
List dreadds() {
  // FIXME STUB

  DreadDs model;

  CharacterVector x = CharacterVector::create("dreadds STUB", model.version);
  List z            = List::create(x ) ;

  return z ;
}
