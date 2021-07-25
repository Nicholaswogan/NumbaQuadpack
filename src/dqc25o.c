#include "cquadpack.h"

#define NMAC 27

double dqc25o(dq_function_type f,double a,double b,double omega,int sincos,
    int nrmom,int maxp1,int ksave,double *abserr,int *neval,
    double *resabs,double *resasc,int *momcom,double **chebmo, void* user_data)
{
    static double x[11] = {
        0.99144486137381041114,
        0.96592582628906828675,
        0.92387953251128675613,
        0.86602540378443864676,
        0.79335334029123516458,
        0.70710678118654752440,
        0.60876142900872063942,
        0.50000000000000000000,
        0.38268343236508977173,
        0.25881904510252076235,
        0.13052619222005159155};
    double ac,an,an2,as,asap,ass,centr,conc,cons,cospar;
    double estc,ests,hlgth,parint,par2,par22;
    double resc12,resc24,ress12,ress24,result,sinpar;
    double cheb12[13],cheb24[25],d[28],d1[28],d2[28];
    double d3[28],fval[25],v[28];
    int unitialized_value = 0xCCCCCCCC;
    double p2 = unitialized_value, p3 = unitialized_value, p4 = unitialized_value;
    int i,isym,j,k,m,noequ,noeq1,mm1;

    centr = 0.5 * (b + a);
    hlgth = 0.5 * (b - a);
    parint = omega * hlgth;

/* Compute the integral using the 15-point Gauss-Kronrod formula
 * if the value of the parameter in the integrand is small or
 * is less than (bb-aa)/2^(maxp1-2), where (aa,bb) is the original
 * integration interval.
 */
    if (fabs(parint) > 2.0) goto _10;
     result = G_K15W(f,dqwgto,omega,p2,p3,p4,sincos,a,b,
              abserr,resabs,resasc, user_data);
     *neval = 15;
     goto _190;

 /* Compute the integral using the generalized Clenshaw-Curtis method. */
 _10:
    conc = hlgth * cos(centr * omega);
    cons = hlgth * sin(centr * omega);
     *resasc = oflow;
    *neval = 25;

/* Check whether the Chebyshev moments for this interval have
 * already been computed.
 */
     if ((nrmom < *momcom) || (ksave == 1)) goto _140;

/* Compute a new set of Chebyshev moments. */
    m = *momcom + 1;
/*** Add variable mm1 to ease transliteration from FORTRAN array
 *** indexing to C indexing.
 ***/
    mm1 = m - 1;
    par2 = parint * parint;
    par22 = par2 + 2.0;
    sinpar = sin(parint);
    cospar = cos(parint);

/* Compute the Chebyshev moments with respect to cosine. */
    v[0] = 2.0 * sinpar / parint;
    v[1] = (8.0 * cospar + (par2 + par2 - 8.0) * sinpar / parint) / par2;
    v[2] = (32.0 * (par2 - 12.0) * cospar + (2.0 * ((par2 - 80.0) *
        par2 + 192.0) * sinpar) / parint) / (par2 * par2);
    ac = 8.0 * cospar;
    as = 24.0 * parint * sinpar;
    if (fabs(parint) > 24.0) goto _70;

/* Compute the Chebyshev moments as the solutions of a boundary
 * value problem with 1 initial value (v[2]) and 1 end value
 * (computed using an asymptotic formula).
 */

     noequ = NMAC - 3;
     noeq1 = noequ - 1;
     an = 6.0;
    for (k = 0; k <= noeq1; k++) {
         an2 = an * an;
         d[k] = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
         d2[k] = (an - 1.0) * (an - 2.0) * par2;
         d1[k] = (an + 3.0) * (an + 4.0) * par2;
         v[k+3] = as - (an2 - 4.0) * ac;
         an += 2.0;
     }
     an2 = an * an;
     d[noequ] = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
     v[noequ+3] = as - (an2 - 4.0) * ac;
     v[3] -= (56.0 * par2 * v[2]);
     ass = parint * sinpar;
    asap = (((((210.0 * par2 -1.0) * cospar - (105.0 * par2 - 63.0) *
         ass) / an2 - (1.0 - 15.0 * par2) * cospar + 15.0 * ass) /
         an2 - cospar + 3.0 * ass) / an2 - cospar) / an2;
     v[noequ+3] -= (2.0 * asap * par2 * (an - 1.0) * (an - 2.0));
/* Solve the tridiagonal system by means of Gaussian elimination
 * with partial pivoting.
 */
     for (i = 0; i <= noequ; i++)
         d3[i] = 0.0;
    d2[noequ] = 0.0;
    for (i = 0; i <= noeq1; i++) {
        if (fabs(d1[i]) <= fabs(d[i])) goto _40;
        an = d1[i];
        d1[i] = d[i];
        d[i] = an;
        an = d2[i];
        d2[i] = d[i+1];
        d[i+1] = an;
        d3[i] = d2[i+1];
        d2[i+1] = 0.0;
        an = v[i+4];
        v[i+4] = v[i+3];
        v[i+3] = an;
_40:
        d[i+1] -= (d2[i] * d1[i] / d[i]);
        d2[i+1] -= (d3[i] * d1[i] / d[i]);
        v[i+4] -= (v[i+3] * d1[i] / d[i]);
    }
    v[noequ+3] /= d[noequ];
    v[noequ+2] = (v[noequ+2] - d2[noeq1] * v[noequ+3]) / d[noeq1];
    for (i = 1; i <= noeq1; i++) {
        k = noequ - i - 1;
        v[k+3] = (v[k+3] - d3[k] * v[k+5] - d2[k] * v[k+4]) / d[k];
    }
    goto _90;

/* Compute the Chebyshev moments by means of forward recursion. */
_70:
    an = 4.0;
    for (i = 3; i < 13; i++) {
        an2 = an * an;
        v[i] = ((an2 - 4.0) * (2.0 * (par22 - an2 - an2) *
            v[i-1] - ac) + as - par2 * (an + 1.0) *
            (an + 2.0) * v[i-2]) / (par2 * (an - 1.0) *
            (an - 2.0));
        an += 2.0;
    }
_90:
    for (j = 0; j < 13; j++)
        chebmo[mm1][2*j] = v[j];

/* Compute the Chebyshev moments with respect to sine. */
    v[0] = 2.0 * (sinpar - parint * cospar) / par2;
    v[1] = (18.0 - 48.0 / par2) * sinpar / par2 +(-2.0 + 48.0 / par2) *
        cospar / parint;
    ac = -24.0 * parint * cospar;
    as = -8.0 * sinpar;
    chebmo[mm1][1] = v[0];
    chebmo[mm1][3] = v[1];
    if (fabs(parint) > 24.0) goto _120;
    for (k = 2; k < 12; k++) {
        an = k+1;
        chebmo[mm1][2*k+1] = - sinpar / (an * (2.0 * an - 2.0)) -
            0.25 * parint * (v[k+1] / an - v[k] / (an - 1.0));
    }
    goto _140;

/* Compute the Chebyshev moments by means of forward recursion. */
_120:
    an = 3.0;
    for (i = 2; i < 12; i++) {
        an2 = an * an;
        v[i] = ((an2 - 4.0) * (2.0 * (par22 - an2 - an2) *
            v[i-1] + as) + ac - par2 * (an + 1.0) *
            (an + 2.0) * v[i-2]) / (par2 * (an - 1.0) *
            (an - 2.0));
        an += 2.0;
        chebmo[mm1][2*i+1] = v[i];
    }
_140:
    if (nrmom < *momcom) {
        m = nrmom + 1;
        mm1 = m - 1;
    }
    if ((*momcom < (maxp1 - 1)) && (nrmom >= (*momcom)))
        (*momcom)++;

/* Compute the coefficients of the Chebyshev expansions of degrees
 * 12 and 24 of the function f.
 */
     fval[0] = 0.5 * f(centr+hlgth, user_data);
     fval[12] = f(centr, user_data);
     fval[24] = 0.5 * f(centr-hlgth, user_data);
     for (i = 1; i < 12; i++) {
         isym = 24 - i;
         fval[i] = f(hlgth*x[i-1]+centr, user_data);
         fval[isym] = f(centr-hlgth*x[i-1], user_data);
     }

     dqcheb(x,fval,cheb12,cheb24);

/* Compute the integral and error estimates. */
    resc12 = cheb12[12] * chebmo[mm1][12];
    ress12 = 0.0;
    estc = fabs(cheb24[24] * chebmo[mm1][24]) + fabs((cheb12[12] -
        cheb24[12]) * chebmo[mm1][12]);
    ests = 0.0;
    k = 10;
    for (j = 0; j < 6; j++) {
        resc12 += (cheb12[k] * chebmo[mm1][k]);
        ress12 += (cheb12[k+1] * chebmo[mm1][k+1]);
        estc += fabs((cheb12[k] - cheb24[k]) * chebmo[mm1][k]);
        ests += fabs((cheb12[k+1] - cheb24[k+1]) * chebmo[mm1][k+1]);
        k -= 2;
    }
    resc24 = cheb24[24] * chebmo[mm1][24];
    ress24 = 0.0;
    *resabs = fabs(cheb24[24]);
    k = 22;
    for (j = 0; j < 12; j++) {
        resc24 += (cheb24[k] * chebmo[mm1][k]);
        ress24 += (cheb24[k+1] * chebmo[mm1][k+1]);
        *resabs += (fabs(cheb24[k]) + fabs(cheb24[k+1]));
        if (j <= 4) {
            estc += (fabs(cheb24[k] * chebmo[mm1][k]));
            ests += (fabs(cheb24[k+1] * chebmo[mm1][k+1]));
        }
        k -= 2;
    }
    *resabs *= fabs(hlgth);
    if (sincos == SINE)
        goto _180;
    result = conc * resc24 - cons * ress24;
    *abserr = fabs(conc * estc) + fabs(cons * ests);
    goto _190;
_180:
    result = conc * ress24 + cons * resc24;
    *abserr = fabs(conc * ests) + fabs(cons * estc);
_190:
     return result;
}

