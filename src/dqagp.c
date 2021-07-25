#include "cquadpack.h"

/* DQAGP - Integration over finite intervals. (From QUADPACK)
 *       Accepts a list of known singularities.
 *
 *    Adaptive integration routine which handles functions
 *    to be integrated between two finite bounds.
 *
 *    The adaptive strategy compares results of integration
 *    over the given interval with the sum of results obtained
 *    from integration over a bisected interval. Since error
 *    estimates are available from each regional integration, the
 *    region with the largest error is bisected and new results
 *    are computed. This bisection process is continued until the
 *    error is less than the prescribed limit or convergence
 *    failure is determined.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    npts2 - number equal to 2 more than the number of sinularities.
 *
 *    points - vector of dimension npts2, the first (npts2-2) elements
 *         of which are the user provided interior break points.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 */
double dqagp(dq_function_type f,double a,double b,int npts2,double *points,
    double epsabs,double epsrel,double *abserr,int *neval,int *ier, void* user_data)
{
    double abseps,alist[LIMIT],area,area1,area12,area2;
    double a1,a2,blist[LIMIT],b1,b2,correc,defabs,defab1;
    double defab2,dres,elist[LIMIT],erlarg,erlast,errbnd;
    double errmax,error1,error2,erro12,errsum,ertest,ndin[40];
    double pts[40],resa,resabs,reseps,result,res3la[3];
    double rlist[LIMIT],rlist2[52],sign,temp;

    int i,id,ierro,ind1,ind2,ip1,iord[LIMIT],iroff1,iroff2,iroff3;
    int j,jlow,jupbnd,k,ksgn,ktmin,last,levcur,level[LIMIT],levmax;
    int maxerr,nint,nintp1,npts,nres,nrmax,numrl2,limit,extrap,noext;

    limit = LIMIT - 1;

/* Test validity of parameters. */
    *ier = 0;
    *neval = 0;
    last = 0;
    result = 0.0;
    *abserr = 0.0;
    alist[0] = a;
    blist[0] = b;
    rlist[0] = 0.0;
    elist[0] = 0.0;
    iord[0] = 0;
    level[0] = 0;
    npts = npts2-2;
    if ((npts2 < 2) || (limit < npts) || ((epsabs < 0.0) &&
        (epsrel < 0.0))) *ier = 6;
    if (*ier == 6) goto _999;

/* If any break points are provided, sort them into an ascending sequence. */
    sign = (a < b) ? 1.0 : -1.0;
    pts[0] = min(a,b);
    if (npts == 0) goto _15;
    for (i = 0; i < npts; i++)
        pts[i+1] = points[i];
_15:
    pts[npts+1] = max(a,b);
    nint = npts + 1;
    a1 = pts[0];
    if (npts == 0) goto _40;
    nintp1 = nint + 1;
    for (i = 0; i < nint; i++) {
        ip1 = i + 1;
        for (j = ip1; j < nintp1; j++) {
            if (pts[i] <= pts[j])
                goto _20;
            temp = pts[i];
            pts[i] = pts[j];
            pts[j] = temp;
_20:
            ;
        }
    }
    if ((pts[0] != min(a,b)) || (pts[nintp1-1] != max(a,b)))
        *ier = 6;
    if (*ier == 6)
        goto _999;

/* Compute first integral and error approximations. */
_40:
    resabs = 0.0;
    for (i = 0; i < nint; i++) {
        b1 = pts[i+1];
        area1 = G_K21(f,a1,b1,&error1,&defabs,&resa, user_data);
        *abserr = *abserr + error1;
        result = result + area1;
        ndin[i] = 0;
        if ((error1 == resa) && (error1 != 0.0))
            ndin[i] = 1;
        resabs += defabs;
        level[i] = 0;
        elist[i] = error1;
        alist[i] = a1;
        blist[i] = b1;
        rlist[i] = area1;
        iord[i] = i;
        a1 = b1;
    }
    errsum = 0.0;
    for (i = 0; i < nint; i++) {
        if (ndin[i] == 1)
            elist[i] = *abserr;
        errsum += elist[i];
    }

/* Test on accuracy. */
/*      last = nint; */
    *neval = 21 * nint;
    dres = fabs(result);
    errbnd = max(epsabs,epsrel*dres);
    if ((*abserr <= 100.0 * epmach * resabs) && (*abserr > errbnd))
        *ier = 2;
    if (nint == 0)
        goto _80;
    for (i = 0; i < npts; i++) {
        jlow = i + 1;
        ind1 = iord[i];
        for (j = jlow; j < nint; j++) {
            ind2 = iord[j];
            if (elist[ind1] > elist[ind2])
                goto _60; /* use continue after debugging */
            ind1 = ind2;
            k = j;
_60:
            ;
        }
        if (ind1 == iord[i])
            goto _70;
        iord[k] = iord[i];
        iord[i] = ind1;
_70:
        ;
    }
    if (limit < npts2)
        *ier = 1;
_80:
    if ((*ier != 0) || (*abserr <= errbnd))
        goto _999;

/* Initialization. */
    res3la[0] = 0.0;
    res3la[1] = 0.0;
    res3la[2] = 0.0;
    rlist2[0] = result;
    maxerr = iord[0];
    errmax = elist[maxerr];
    area = result;
    nrmax = 0;
    nrmax = 0;
    nres = -1;            /* nres = 0 */
    numrl2 = 0;            /* numrl2 = 1 */
    ktmin = 0;
    extrap = FALSE;
    noext = FALSE;
    erlarg = errsum;
    ertest = errbnd;
    levmax = 1;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ierro = 0;
    *abserr = oflow;
    ksgn = -1;
    if (dres > (1.0 - 50.0 * epmach) * resabs)
        ksgn = 1;

/* Main loop. */
    for (last = npts2; last <= limit; last++) {

/* Bisect the interval with the nrmax-th largest error estimate. */
        levcur = level[maxerr] + 1;
        a1 = alist[maxerr];
        b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        area1 = G_K21(f,a1,b1,&error1,&resa,&defab1, user_data);
        area2 = G_K21(f,a2,b2,&error2,&resa,&defab2, user_data);
/* Improve previous approximations to integral and error
      and test for accuracy. */
          *neval += 42;
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if ((defab1 == error1) || (defab2 == error2)) goto _95;
        if ((fabs(rlist[maxerr] - area12) > 1.0e-5 * fabs(area12))
            || (erro12 < .99 * errmax)) goto _90;
        if (extrap) iroff2++;
        else iroff1++;
_90:
        if ((last > 9) && (erro12 > errmax))    /* last > 10 */
            iroff3++;
_95:
        level[maxerr] = levcur;
        level[last] = levcur;
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
        if (error2 > error1) goto _100;
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
        goto _110;
_100:
        alist[maxerr] = a2;
        alist[last] = a1;
        blist[last] = b1;
        rlist[maxerr] = area2;
        rlist[last] = area1;
        elist[maxerr] = error2;
        elist[last] = error1;

/* Call dqsort to maintain the descending ordering in the list of error
    estimates and select the subinterval with nrmax-th largest
    error estimate (to be bisected next). */
_110:
        dqsort(limit,last,&maxerr,&errmax,elist,iord,&nrmax);
        if (errsum <= errbnd) goto _190;
        if (*ier != 0) goto _170;
        if (noext) goto _160;
        erlarg -= erlast;
        if (levcur+1 <= levmax)
            erlarg += erro12;
        if (extrap) goto _120;

/* Test whether the interval to be bisected next is the smallest interval. */
        if ((level[maxerr]+1) <= levmax)
            goto _160;
        extrap = TRUE;
        nrmax = 1;        /* nrmax = 2 */
_120:
        if ((ierro == 3) || (erlarg <= ertest)) goto _140;

/* The smallest interval has the largest error. Before bisecting, decrease
    the sum of the errors over the larger intervals (erlarg) and
        perform extrapolation.) */
        id = nrmax;
        jupbnd = last;
        if (last > (2 + limit/2))
            jupbnd = limit + 3 - last;
        for (k = id;k <= jupbnd; k++) {
            maxerr = iord[nrmax];
            errmax = elist[maxerr];
            if (level[maxerr]+1 <= levmax)
                goto _160;
            nrmax++;
        }

/* Perform extrapolation. */
_140:
        numrl2++;
        rlist2[numrl2] = area;
        if (numrl2 <= 1) goto _155;
        reseps=dqext(&numrl2,rlist2,&abseps,res3la,&nres);
        ktmin++;
        if ((ktmin > 5) && (*abserr < 1.0e-3 * errsum)) *ier = 5;
        if (abseps >= *abserr) goto _150;
        ktmin = 0;
        *abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs,epsrel * fabs(reseps));
        if (*abserr <= ertest) goto _170;

/* Prepare bisection of the smallest interval. */
_150:
        if (numrl2 == 0) noext = TRUE;
        if (*ier == 5) goto _170;
_155:
        maxerr = iord[0];
        errmax = elist[maxerr];
        nrmax = 0;
        extrap = FALSE;
        levmax += 1.0;
        erlarg = errsum;
_160:
        ;
    }
_170:
    if (*abserr == oflow) goto _190;
    if ((*ier + ierro) == 0) goto _180;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if ((result != 0.0) && (area != 0.0)) goto _175;
    if (*abserr > errsum) goto _190;
    if (area == 0.0) goto _210;
    goto _180;
_175:
    if (*abserr/fabs(result) > errsum/fabs(area)) goto _190;

/* Test on divergence. */
_180:
    if ((ksgn == -1) && (max(fabs(result),fabs(area)) <= defabs * .01))
        goto _210;
    if ((0.01 > result/area) || (result/area > 100.0) ||
        (errsum > fabs(area))) *ier = 6;
    goto _210;

/* Compute global integral. */
_190:
    result = 0.0;
    for (k = 0; k <= last; k++)
        result += rlist[k];
    *abserr = errsum;
_210:
    if (*ier > 2) (*ier)--;
    result = result * sign;
_999:
    return result;
}
