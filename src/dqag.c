/* dqag.c -- modified version of QUADPACK routine DQAG.
 * (C)1999, C. Bond. All right reserved.
 *
 * There are no changes to the basic computational method. Only
 * the temporary storage strategy is changed to utilize the
 * local stack at the appropriate level. This reduces the
 * need for memory allocation of arrays at higher levels and
 * the resulting passing of memory pointers down the line.
 *
 */
#include "cquadpack.h"

/* DQAG - Approximation to definite integral. (From QUADPACK)
 *
 *  Calls DQAGE with appropriate parameters assigned.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 *
 *    irule - integration rule to be used as follows:
 *        irule = 1 -- G_K 7-15
 *        irule = 2 -- G_K 10-21
 *        irule = 3 -- G_K 15-31
 *        irule = 4 -- G_K 20-41
 *        irule = 5 -- G_K 25-51
 *        irule = 6 -- G_K 30-61
 */
double dqag(dq_function_type f,double a,double b,double epsabs,
    double epsrel,int irule,double *abserr,int *neval,int *ier, void* user_data)
{
    double result;
    int last;

    result = dqage(f,a,b,epsabs,epsrel,irule,abserr,neval,ier,&last, user_data);

    return result;
}
