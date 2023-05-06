#include "cquadpack.h"

/* function to perform log(exp(lna) - exp(lnb)) maintaining numerical precision */
double logsubexp(const double x, const double y){
  double tmp = x - y;
  if ( x > y && fabs(tmp) > 1e3*GSL_DBL_EPSILON ){ /* numbers smaller than this can just give numerical noise in the gsl_sf_log_1plusx function */
    return x + gsl_sf_log_1plusx(-exp(-tmp));
  }
  else{
    return -INFINITY;
  }
}