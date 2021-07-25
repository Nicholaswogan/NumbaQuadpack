#include "cquadpack.h"

/*  DQAWCE - Computation of Cauchy principal value
 *
 *  PARAMETERS:
 *
 *      f() -   double precision function defining the integrand.
 *
 *      a   -   lower limit of integration
 *
 *      b   -   upper limit of integration
 *
 *      c   -   parameter in the weight function
 *
 *      epsabs  -   absolute accuracy requested
 *
 *      epsrel  -   relative accuracy requested
 *
 *      abserr  -   estimate of the modulus of the absolute error
 *
 *      neval   -   number of function evaluations
 *
 *      ier     -   error code
 */
double dqawce(dq_function_type f,double a,double b,double c,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
    double aa,area,area1,area2,area12,a1,a2,bb,b1,b2;
    double errbnd,errmax,error1,error2,erro12,errsum,result;
    double alist[LIMIT],blist[LIMIT],rlist[LIMIT],elist[LIMIT];

    int iord[LIMIT],iroff1,iroff2,k,krule,last,maxerr,nrmax,nev;
    int limit;

    limit = LIMIT - 1;
    *ier = 6;
    *neval = 0;
    last = 0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    result = 0.0;
    *abserr = 0.0;
    if ((c==a) || (c==b) || ((epsabs < 0.0) && (epsrel < 0.0))) goto _999;

/*  First approximation to the integral.    */
    aa = a;
    bb = b;
    if (a <= b) goto _10;
    aa = b;
    bb = a;
_10:
    *ier = 0;
    krule = 1;
    result = dqc25c(f,aa,bb,c,abserr,&krule,neval, user_data);
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;
    alist[0] = a;
    blist[0] = b;

/*  Test on accuracy.   */
    errbnd = max(epsabs,epsrel * fabs(result));
    if (limit == 0) *ier = 1;
    if ((*abserr < min(1.0e-2 * fabs(result),errbnd))  || (*ier == 1))
        goto _70;

/*  Initialization. */
    alist[0] = aa;
    blist[0] = bb;
    rlist[0] = result;
    errmax = *abserr;
    maxerr = 0;
    area = result;
    errsum = *abserr;
    nrmax = 0;
    iroff1 = 0;
    iroff2 = 0;

/*  Main loop.  */
    for (last = 1;last < limit;last++) {
/* Bisect the subinterval with nrmax-th largest error estimate.    */
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr]+blist[maxerr]);
        b2 = blist[maxerr];
        if ((c <= b1) && (c > a1)) b1 = 0.5 * (c + b2);
        if ((c >  b1) && (c < b2)) b1 = 0.5 * (a1 + c);
        a2 = b1;
        krule = 2;
        area1 = dqc25c(f,a1,b1,c,&error1,&krule,&nev, user_data);
        *neval = *neval + nev;
        area2 = dqc25c(f,a2,b2,c,&error2,&krule,&nev, user_data);
        *neval = *neval + nev;

/*  Improve previous approximations to integral and error and
 *  test for accuracy.
 */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum += (erro12-errmax);
        area += (area12-rlist[maxerr]);
        if (fabs(rlist[maxerr]-area12) < (1.0e-5*fabs(area12)) &&
            (erro12 >= 0.99 * errmax) && (krule == 0)) iroff1++;
        if ((last > 10) && (erro12 > errmax) && (krule == 0))
            iroff2++;
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel*fabs(area));
        if (errsum <= errbnd) goto _15;

/*  Test for roundoff error and eventually set error flag. */
        if ((iroff1 >= 6) && (iroff2 > 20)) *ier = 2;

/*  Set error flag in the case that number of interval bisections
 *  exceeds limit.
 */
        if (last == limit) *ier = 1;

/*  Set error flag in the case of bad integrand behavior at a point
 *  of the integration range.
 */
        if (max(fabs(a1),fabs(b2)) <= (1.0 + 1.0e3*epmach)*
            (fabs(a2) + 1.0e3 * uflow)) *ier = 3;
/* Append newly created intervals to list. */
_15:
        if (error2 <= error1) {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = error1;
            elist[last] = error2;
        }
        else {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = area2;
            rlist[last] = area1;
            elist[maxerr] = error2;
            elist[last] = error1;
        }

/*  Call subroutine dqsort to maintain descending ordering in the list. */
        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);

/*  Jump out of loop.   */
        if ((*ier != 0) || (errsum <= errbnd)) goto _50;
    }

/*  Compute final result.   */
_50:
    result = 0.0;
    for (k=0;k<=last;k++) {
        result += rlist[k];
    }
    *abserr = errsum;
_70:
    if (aa == b) result = -result;
_999:
    return result;
}

