#include "cquadpack.h"

/* DQAGI - Integration over (semi-) infinite intervals. (From QUADPACK)
 *
 *    Adaptive integration routine which handles functions
 *    to be integrated between -infinity to +infinity, or
 *    between either of those limits and some finite,
 *    real boundary.
 *
 *    The adaptive strategy compares results of integration
 *    over the interval with the sum of results obtained from
 *    integration of bisected interval. Since error estimates
 *    are available from each regional integration, the interval
 *    with the largest error is bisected and new results are
 *    computed. This bisection process is continued until the
 *    error is less than the prescribed limit or convergence
 *    failure is determined.
 *
 *    Note that bisection, in the sense used above, refers to
 *    bisection of the transformed interval.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    bound - optional finite bound on integral.
 *
 *    inf - specifies range of integration as follows:
 *        inf = -1 -- range is from -infinity to bound,
 *        inf =  1 -- range is from bound to +infinity,
 *        inf =  2 -- range is from -infinity to +infinity,
 *                (bound is immaterial in this case).
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 */
double dqagi(dq_function_type f,double bound,int inf,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
    double abseps,area,area1,area12,area2,a1,a2,b1,b2;
    double boun,correc,defabs,defab1,defab2,dres,erlarg;
    double erlast,errbnd,errmax,error1,error2,erro12;
    double errsum,ertest,resabs,reseps,result,res3la[3];
    double alist[LIMIT],blist[LIMIT],elist[LIMIT],rlist[LIMIT];
    double rlist2[52],small = 0; /* small will be initialized in _80 */

    int id,ierro,iord[LIMIT],iroff1,iroff2,iroff3,jupbnd,k,ksgn;
    int ktmin,last,maxerr,nres,nrmax,numrl2;
    int limit,extrap,noext;

    limit = LIMIT - 1;
/* Test validity of parameters. */
    *ier = 0;
    *neval = 0;
    last = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = 0.0;
    blist[0] = 1.0;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    if ((epsabs < 0.0) && (epsrel < 0.0)) *ier = 6;
    if ((inf != 1) && (inf != -1) && (inf != 2)) *ier = 6;
    if (*ier == 6) return result;

/* First approximation to the integral. */
    boun = bound;
    if (inf == 2) boun = 0.0;

    result = G_K15I(f,boun,inf,0.0,1.0,abserr,&defabs,&resabs, user_data);

/* Test on accuracy. */
    last = 0;
    rlist[0] = result;
    elist[0] = *abserr;
    iord[0] = 0;
    dres = fabs(result);
    errbnd = max(epsabs,epsrel*dres);
    if ((*abserr <= 100.0 * epmach * defabs) && (*abserr > errbnd))
        *ier = 2;
    if (limit == 0) *ier = 1;
    if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) ||
        (*abserr == 0.0)) goto _130;

/* Initialization for main loop. */
    rlist2[0] = result;
    errmax = *abserr;
    maxerr = 0;             /* maxerr = 1 */
    area = result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 0;
    nres = 0;          /* nres = 0 */
    ktmin = 0;
    numrl2 = 1;            /* numrl2 = 2 */
    extrap = FALSE;
    noext = FALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres > (1.0 - 50.0 * epmach) * defabs)
        ksgn = 1;

/* Main loop. */
    for (last = 1; last <= limit; last++) {
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        area1 = G_K15I(f,boun,inf,a1,b1,&error1,&resabs,&defab1, user_data);
        area2 = G_K15I(f,boun,inf,a2,b2,&error2,&resabs,&defab2, user_data);

/* Improve previous approxminations to integral and error
      and test for accuracy. */
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 == error1) || (defab2 == error2)) goto _15;
        if ((fabs(rlist[maxerr] - area12) > 1.0e-5 * fabs(area12))
            || (erro12 < .99 * errmax)) goto _10;
        if (extrap) iroff2++;
        else iroff1++;
_10:
        if ((last > 9) && (erro12 > errmax))    /* last > 10 */
            iroff3++;
_15:
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel * fabs(area));

/* Test for roundoff error and eventually set error flag. */
        if (((iroff1 + iroff2) >= 10) || (iroff3 >= 20))
            *ier = 2;
        if (iroff2 > 5)
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
        if (errsum <= errbnd) goto _115;
        if (*ier != 0) goto _100;
        if (last == 1) goto _80;    /* last == 2 */
        if (noext) continue;  //goto _90;
        erlarg -= erlast;
        if (fabs(b1-a1) > small)
            erlarg += erro12;
        if (extrap) goto _40;

/* Test whether the interval to be bisected next is the smallest interval. */
        if ((fabs(blist[maxerr] - alist[maxerr])) > small)
            goto _90;
        extrap = TRUE;
        nrmax = 1;        /* nrmax = 2 */
_40:
        if ((ierro == 3) || (erlarg <= ertest)) goto _60;

/* The smallest interval has the largest error. Before bisecting, decrease
    the sum of the erors over the larger intervals (erlarg) and
        perform extrapolation.) */
        id = nrmax;
        jupbnd = last;
        if (last > (2 + limit/2))
            jupbnd = limit + 3 - last;
        for (k = id;k <= jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (fabs(blist[maxerr] - alist[maxerr]) > small)
                goto _90;
            nrmax++;
        }

/* Perform extrapolation. */
_60:
        numrl2++;
        rlist2[numrl2] = area;
        reseps=dqext(&numrl2,rlist2,&abseps,res3la,&nres);
        ktmin++;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) *ier = 5;
        if (abseps >= *abserr) goto _70;
        ktmin = 0;
        *abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs,epsrel * fabs(reseps));
        if (*abserr <= ertest) goto _100;

/* Prepare bisection of the smallest interval. */
_70:
        if (numrl2 == 0) noext = TRUE;
        if (*ier == 5) goto _100;
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = FALSE;
        small = small * 0.5;
        erlarg = errsum;
        continue;
_80:
        small = .375;
        erlarg = errsum;
        ertest = errbnd;
        rlist2[1] = area;
_90:
        ;
    }                    /* 90: */
_100:
    if (*abserr == oflow) goto _115;
    if ((*ier + ierro) == 0) goto _110;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if ((result != 0.0) && (area != 0.0)) goto _105;
    if (*abserr > errsum) goto _115;
    if (area == 0.0) goto _130;
    goto _110;
_105:
    if (*abserr/fabs(result) > errsum/fabs(area)) goto _115;

/* Test on divergence. */
_110:
    if ((ksgn == -1) && (max(fabs(result),fabs(area)) <= defabs * .01))
        goto _130;
    if ((0.01 > result/area) || (result/area > 100.0) ||
        (errsum > fabs(area))) *ier = 6;
    goto _130;

/* Compute global integral. */
_115:
    result = 0.0;
    for (k = 0; k <= last; k++)
        result += rlist[k];
    *abserr = errsum;
_130:
    *neval = 30 * last + 15;
    if (inf == 2) *neval *= 2;
    if (*ier > 2) (*ier)--;
    return result;
}
