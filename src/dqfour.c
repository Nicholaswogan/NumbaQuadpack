#include "cquadpack.h"

/* DQFOUR - Computation of oscillatory integrals.
 *
 *    Calculates an approximation to a given definite integral.
 *        I = integral of F(X) * W(X) over (A,B)
 *            where W(X) = COS(OMEGA * X)
 *               or W(X) = SIN(OMEGA * X)
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    omega - parameter in the integrand weight function.
 *
 *    sincos - indicates which weight function to use:
 *        sincos = COSINE,    W(X) = COS(OMEGA * X)
 *        sincos = SINE,        W(X) = SIN(OMEGA * X)
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 *
 *    limit - upper bound on number of of subdivisions in
 *        the partition of (A,B)
 *
 */

double dqfour(dq_function_type f,double a,double b,double omega,
    int sincos,double epsabs,double epsrel,
    int icall,int maxp1,double *abserr,
    int *neval,int *ier,
    int *momcom,double **chebmo, void* user_data)
{
    double abseps,area,area1,area12,area2;
    double a1,a2,b1,b2,correc,defabs,defab1;
    double defab2,domega,dres,erlarg,erlast,errbnd;
    double errmax,error1,error2,erro12,errsum,ertest;
    double resabs,reseps,result,res3la[3];
    double alist[LIMIT],blist[LIMIT],rlist[LIMIT];
    double elist[LIMIT],rlist2[52],small,width;

    int id,ierro,iroff1,iroff2,iroff3,jupbnd,k,ksgn,limit;
    int ktmin,last,maxerr,nev,nres,nrmax,nrmom,numrl2;
    int extrap,noext,extall,iord[LIMIT],nnlog[LIMIT];

    limit = LIMIT - 1;
/* Test validity of parameters. */
    *ier = 0;
    *neval = 0;
//    last = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    nnlog[0] = 0;
    if (((sincos != COSINE) && (sincos != SINE)) || ((epsabs < 0.0) &&
        (epsrel < 0.0)) || (icall < 1) || (maxp1 < 1)) *ier = 6;
    if (*ier == 6) return result;

/* First approximation to the integral. */
    domega = fabs(omega);
    nrmom = 0;
    if (icall <= 1)
        *momcom = 0;
_5:
    result = dqc25o(f,a,b,domega,sincos,nrmom,maxp1,0,
        abserr,neval,&defabs,&resabs,momcom,chebmo, user_data);
/* Test on accuracy. */
    dres = fabs(result);
    errbnd = max(epsabs,epsrel*dres);
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;
    if ((*abserr <= 100.0 * epmach * defabs) && (*abserr > errbnd))
        *ier = 2;
    if (limit == 0) *ier = 1;
    if ((*ier != 0) || (*abserr <= errbnd) && (*abserr != resabs) ||
        (*abserr == 0.0)) goto _200;

/* Initialization. */
    errmax = *abserr;
    maxerr = 0;             /* maxerr = 1 */
    area = result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 0;
    extrap = FALSE;
    noext = FALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ktmin = 0;
    small = fabs(b-a) * 0.75;
    nres = 0;
    numrl2 = -1;
    extall = FALSE;
    if ((0.5 * fabs(b-a) * domega) > 2.0)
        goto _10;
    numrl2 = 0;
    extall = TRUE;
    rlist2[0] = result;
_10:
    if ((0.25 * fabs(b-a) * domega) <= 2.0)
        extall = TRUE;
    ksgn = -1;
    if (dres > (1.0 - 50.0 * epmach) * defabs)
        ksgn = 1;

/* Main loop. */
    for (last = 1; last < limit; last++) {

/* Bisect the interval with the nrmax-th largest error estimate. */
        nrmom = nnlog[maxerr] + 1;
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        area1 = dqc25o(f,a1,b1,domega,sincos,nrmom,maxp1,0,
            &error1,&nev,&resabs,&defab1,momcom,chebmo, user_data);
        *neval += nev;
        area2 = dqc25o(f,a2,b2,domega,sincos,nrmom,maxp1,1,
            &error2,&nev,&resabs,&defab2,momcom,chebmo, user_data);
        *neval += nev;

/* Improve previous approximations to integral and error
      and test for accuracy. */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 == error1) || (defab2 == error2)) goto _25;
        if ((fabs(rlist[maxerr] - area12) > 1.0e-5 * fabs(area12))
            || (erro12 < .99 * errmax)) goto _20;
        if (extrap) iroff2++;
        else iroff1++;
_20:
        if ((last > 9) && (erro12 > errmax))    /* last > 10 */
            iroff3++;
_25:
        rlist[maxerr] = area1;
        rlist[last] = area2;
        nnlog[maxerr] = nrmom;
        nnlog[last] = nrmom;
        errbnd = max(epsabs,epsrel * fabs(area));

/* Test for roundoff error and eventually set error flag. */
        if (((iroff1 + iroff2) >= 10) || (iroff3 >= 20))
            *ier = 2;
        if (iroff2 >= 5)
            *ier = 3;

/* Set error flag in the case that the number of subintervals
    equals limit. */
        if (last == limit)    /* last == limit */
            *ier = 1;

/* Set error flag in the case of bad integrand behavior at some
    points in the integration range. */
        if (max(fabs(a1),fabs(b2)) <= (1.0 +1000.0 * epmach) *
            (fabs(a2) + 1000.0*uflow))
            *ier = 4;

/* Append the newly-created intervals to the list. */
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
/* Call dqsort to maintain the descending ordering in the list of error
    estimates and select the subinterval with nrmax-th largest
    error estimate (to be bisected next). */

        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);
        if (errsum <= errbnd) goto _170;
        if (*ier != 0) goto _150;
        if ((last == 1) && (extall)) goto _120;    /* last == 2 */
        if (noext) goto _140;
        if (!extall) goto _50;
        erlarg -= erlast;
        if (fabs(b1-a1) > small)
            erlarg += erro12;
        if (extrap) goto _70;

/* Test whether the interval to be bisected next is the smallest interval. */
_50:
        width = fabs(blist[maxerr] - alist[maxerr]);
        if (width > small)
            goto _140;
        if (extall)
            goto _60;

/* Test whether we can start with the extrapolation procedure (we do
 * this if we integrate over the next interval with use of a Gauss-
 * Kronrod rule) - see routine dqc25o. */
         small *= 0.5;
         if ((0.25 * width * domega) > 2.0)
             goto _140;
         extall = TRUE;
         goto _130;
_60:
        extrap = TRUE;
        nrmax = 1;        /* FORTRAN: nrmax = 2 */
_70:
        if ((ierro == 3) || (erlarg <= ertest))
            goto _90;

/* The smallest interval has the largest error. Before bisecting, decrease
    the sum of the erorrs over the larger intervals (erlarg) and
        perform extrapolation. */
        jupbnd = last;
        if (last > (2 + limit/2))
            jupbnd = limit + 3 - last;
        id = nrmax;
        for (k = id;k <= jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (fabs(blist[maxerr] - alist[maxerr]) > small)
                goto _140;
            nrmax++;
        }

/* Perform extrapolation. */
_90:
        numrl2++;
        rlist2[numrl2] = area;
        if (numrl2 < 2)
            goto _110;
        reseps = dqext(&numrl2,rlist2,&abseps,res3la,&nres);
        ktmin++;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) *ier = 5;
        if (abseps >= *abserr) goto _100;
        ktmin = 0;
        *abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs,epsrel * fabs(reseps));
        if (*abserr <= ertest) goto _150;

/* Prepare bisection of the smallest interval. */
_100:
        if (numrl2 == 0) noext = TRUE;
        if (*ier == 5) goto _150;
_110:
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = FALSE;
        small *= 0.5;
        erlarg = errsum;
        goto _140;
_120:
        small *= 0.5;
        numrl2++;
        rlist2[numrl2] = area;
_130:
        erlarg = errsum;
        ertest = errbnd;
_140:
        ;
    }

/* Set the final result. */
_150:
    if ((*abserr == oflow) || (nres == 0)) goto _170;
    if ((*ier + ierro) == 0) goto _165;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if ((result != 0.0) && (area != 0.0)) goto _160;
    if (*abserr > errsum) goto _170;
    if (area == 0.0) goto _190;
    goto _165;
_160:
    if (*abserr/fabs(result) > errsum/fabs(area)) goto _170;

/* Test on divergence. */
_165:
    if ((ksgn == -1) && (max(fabs(result),fabs(area)) <= defabs * .01))
        goto _190;
    if ((0.01 > result/area) || (result/area > 100.0) ||
        (errsum > fabs(area))) *ier = 6;
    goto _190;

/* Compute global integral. */
_170:
    result = 0.0;
    for (k = 0; k <= last; k++)
        result += rlist[k];
    *abserr = errsum;
_190:
    if (*ier > 2) (*ier)--;
_200:
    if ((sincos == SINE) && (omega < 0.0))
        result = - result;
    return result;
}
