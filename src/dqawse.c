#include "cquadpack.h"

/*  DQAWSE - Approximation to integral with algebraic and/or logarithmic
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
double dqawse(dq_function_type f,double a,double b,double alfa,double beta,
        int wgtfunc,double epsabs,double epsrel,double *abserr,
        int *neval,int *ier, void* user_data)
{
    double alist[LIMIT],blist[LIMIT],rlist[LIMIT],elist[LIMIT];
    double ri[25],rj[25],rh[25],rg[25];
    double area,area1,area12,area2,a1,a2,b1,b2,centre;
    double errbnd,errmax,error1,erro12,error2,errsum;
    double resas1,resas2,result;

    int iord[LIMIT],iroff1,iroff2,k,last,limit,maxerr,nev,nrmax;

    limit = LIMIT - 1;
/*  Test on validity of parameters. */
    *ier = 6;
    *neval = 0;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    result = 0.0;
    *abserr = 0.0;
    if ((b <= a) || ((epsabs < 0.0) && (epsrel < 0.0)) ||
        (alfa <= -1.0) || (beta <= -1.0) || (wgtfunc < 1) ||
        (wgtfunc > 4) || (limit < 1)) goto _999;
    *ier = 0;

/*  Compute the modified Chebyshev moments. */
    dqmomo(alfa,beta,ri,rj,rg,rh,wgtfunc);

/*  Integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b). */
    centre = 0.5 * (a+b);
    area1 = dqc25s(f,a,b,a,centre,alfa,beta,ri,rj,rg,rh,&error1,
        &resas1,wgtfunc,&nev, user_data);
    *neval = *neval + nev;
    area2 = dqc25s(f,a,b,centre,b,alfa,beta,ri,rj,rg,rh,&error2,
        &resas2,wgtfunc,&nev, user_data);
    *neval = *neval + nev;
    result = area1 + area2;
    *abserr = error1 + error2;

/* Test on accuracy. */
    errbnd = max(epsabs,epsrel * fabs(result));

/*  Initialization. */
    if (error1 >= error2) {
        alist[0] = a;
        alist[1] = centre;
        blist[0] = centre;
        blist[1] = b;
        rlist[0] = area1;
        rlist[1] = area2;
        elist[0] = error1;
        elist[1] = error2;
    }
    else {
        alist[0] = centre;
        alist[1] = a;
        blist[0] = b;
        blist[1] = centre;
        rlist[0] = area2;
        rlist[1] = area1;
        elist[0] = error2;
        elist[1] = error1;
    }
    iord[0] = 0;
    iord[1] = 1;
    if (limit == 1) *ier = 1;
    if ((*abserr <= errbnd) || (*ier == 1)) goto _999;
    errmax = elist[0];
    maxerr = 0;
    nrmax = 0;
    area = result;
    errsum = maxerr;
    iroff1 = 0;
    iroff2 = 0;

/*  Main loop. */
    for (last = 2;last < limit;last++) {
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr]+blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];

        area1 = dqc25s(f,a,b,a1,b1,alfa,beta,ri,rj,rg,rh,&error1,
                &resas1,wgtfunc,&nev, user_data);
        *neval = *neval + nev;
        area2 = dqc25s(f,a,b,a2,b2,alfa,beta,ri,rj,rg,rh,&error2,
                &resas2,wgtfunc,&nev, user_data);
        *neval = *neval + nev;

/*  Improve previous approximation and error test for accuracy. */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum += (erro12 - errmax);
        area += (area12-rlist[maxerr]);
        if ((a == a1) || (b == b2)) goto _30;
        if ((resas1 == error1) || (resas2 == error2)) goto _30;

/*  Test for roundoff error. */
        if ((fabs(rlist[maxerr]-area12) < (1.0e-5 * fabs(area12))) &&
            (erro12 >= (0.99 *errmax))) iroff1++;
        if ((last > 9) && (erro12 > errmax)) iroff2++;
_30:
        rlist[maxerr] = area1;
        rlist[last] = area2;

/*  Test on accuracy. */
        errbnd = max(epsabs,epsrel*fabs(area));
        if (errsum <= errbnd) goto _35;

/*  Set error flag in the case that number of intervals exceeds limit. */
        if (last == limit) *ier = 1;

/*  Set error flag in the case of roundoff error. */
        if ((iroff1 > 5) || (iroff2 > 19)) *ier = 2;

/*  Set error flag in case of bad integrand behavior at interior points. */
        if ( max(fabs(a1),fabs(b2)) <= ((1.0 + 1.0e3 * epmach) *
                (fabs(a2)+1.0e3 * uflow)) ) *ier = 3;
/*  Append the newly created intervals to the list. */
_35:
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

/*  Call subroutine qsort to maintain the descending ordering. */
        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);

/*  Jump out of loop. */
    if ((*ier != 0) || (errsum <= errbnd)) break;
    }
    result = 0.0;
    for (k=0;k<=last;k++) {
        result += rlist[k];
    }
    *abserr = errsum;
_999:
    return result;
}
