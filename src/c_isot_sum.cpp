#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;


double taper_sqExp(double x1, double x2, double a){
  return exp( -a/4.0 * (pow(x1,2) + pow(x2,2)) );
}

double taper1(double x1, double x2, double a){
  return 1.0;
}


// [[Rcpp::export]]
NumericVector c_iso_sum(NumericMatrix x, NumericVector t, NumericVector sl, double taper_a) {
  
  int nt = t.length();
  NumericVector out( nt );
  
  int n = x.nrow(), i, j, k;
  
  double d, dx, dy;
  
  double hij;
  
  double (*taper)(double, double, double) = &(taper1);
  if(taper_a > 0) taper = &taper_sqExp;

  for(i = 0; i < n - 1 ; i++ ) {
    Rcpp::checkUserInterrupt();
    for(j = i + 1; j < n ; j++) {
      dx = x(i,0)-x(j,0);
      dy = x(i,1)-x(j,1);
      d = sqrt( dx*dx + dy*dy) * M_2PI;
      hij = taper(dx/sl(0), dy/sl(1), taper_a );
      for(k = 0; k < nt; k++)
          out(k) += R::bessel_j(d * t(k), 0.0) * hij;
    }
  }
  return out;
}
