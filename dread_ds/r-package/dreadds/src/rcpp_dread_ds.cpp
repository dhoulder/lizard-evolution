
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List dreadds() {
  // FIXME STUB
    CharacterVector x = CharacterVector::create("dreadds STUB");
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}
