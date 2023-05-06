#include "cquadpack.h"

double logsumexp(const double x, const double y){
  double tmp = x - y;
  if ( x == y || fabs(tmp) < 1e3*GSL_DBL_EPSILON ){ return x + M_LN2; } /* require the x == y to deal with cases when x and y are both -inf */
  else{
    if ( tmp > 0. ){
      return x + gsl_sf_log_1plusx(exp(-tmp));
    }
    else if ( tmp <= 0. ){
      return y + gsl_sf_log_1plusx(exp(tmp));
    }
    else{
      return tmp;
    }
  }
}