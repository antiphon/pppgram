#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;


double taper_sqExp(NumericVector xy, double a){
  return exp( -a * pow(xy(0)/2,2) ) * exp( -a * pow(xy(1)/2,2) );
}

double taper1(NumericVector xy, double a){
  return 1.0;
}


// [[Rcpp::export]]
NumericVector c_iso_sum(NumericMatrix x, NumericVector t, NumericVector sl, double taper_a) {
  
  int nt = t.length();
  NumericVector out( nt );
  
  int n = x.nrow(), i, j, k;
  
  double d, dx, dy;
  
  double hi, hj;
  
  //double maxd = sqrt(2) * 2 * taper_a * sl(0) * sl(1);
  
  double (*taper)(NumericVector, double) = &(taper1);
  if(taper_a > 0) taper = &taper_sqExp;
  //   fh1d <- function(x, a, l) exp( -a * x^2 / (2*l)^2)
  //fry_h2d <- function(x, a, l) fh1d(x[,1], a, l[1])*fh1d(x[,2], a, l[2]) / prod(l)
    
  
  for(i = 0; i < n - 1 ; i++ ) {
    hi = taper(x(i,_)/sl, taper_a );
    for(j = i + 1; j < n ; j++) {
      dx = x(i,0)-x(j,0);
      dy = x(i,1)-x(j,1);
      d = sqrt( dx*dx + dy*dy);
        hj = taper(x(i,_)/sl(0), taper_a ) * hi;
        for(k = 0; k < nt; k++)
            out(k) += R::bessel_j(2 * PI * d * t(k), 0.0)* hj;
    }
  }
  return out;
}
