#include "cquadpack.h"

/*  DQAWS - Approximation to integral with algebraic and/or logarithmic
 *          singularities.
 *
 *  PARAMETERS:
 *
 *      f() - double precision function to be integrated.
 *
 *      a   - double lower limit of integration.
 *
 *      b   - upper limit of integration.
 *
 *      alfa - parameter in the weight function.
 *
 *      beta - parameter in the weight function.
 *
 *      wgtfunc - indicates which weight function is to be used.
 *                  = 1:    (x-a)^alfa * (b-x)^beta
 *                  = 2:    (x-a)^alfa * (b-x)^beta * log(x-a)
 *                  = 3:    (x-a)^alfa * (b-x)^beta * log(b-x)
 *                  = 4:    (x-a)^alfa * (b-x)^beta * log(x-a) * log(b-x)
 *
 *      epsabs  - absolute accuracy requested.
 *
 *      epsrel  - relative accuracy requested.
 *
 */
double dqaws(dq_function_type f,double a,double b,double alfa,double beta,
        int wgtfunc,double epsabs,double epsrel,double *abserr,
        int *neval,int *ier, void* user_data)
{
    double result;

    result = dqawse(f,a,b,alfa,beta,wgtfunc,epsabs,epsrel,abserr,
                neval,ier, user_data);
    return result;
}
