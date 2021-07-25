#include "cquadpack.h"

double dqc25s(dq_function_type f,double a,double b,double bl,double br,
    double alfa,double beta,double ri[],double rj[],double rg[],
    double rh[],double *abserr,double *resasc,int wgtfunc,int *nev, void* user_data)
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
    double centr,dc,factor,fix,hlgth,resabs,res12,res24,u,result;
    double cheb12[13],cheb24[25],fval[25];
    int i,isym;

    *nev = 25;
    if ((bl == a) && ((alfa != 0.0) || (wgtfunc == 2) || (wgtfunc == 4)))
        goto _10;
    if ((br == b) && ((beta != 0.0) || (wgtfunc == 3) || (wgtfunc == 4)))
        goto _140;

/*  If a>bl and b<br, apply the 15-point Gauss-Kronrod scheme. */
    result = G_K15W(f,dqwgts,a,b,alfa,beta,wgtfunc,bl,br,abserr,
                &resabs,resasc, user_data);
    *nev = 15;
    goto _270;

/*  This part is only executed if a = bl.
 *  Compute the Chebyshev series expansion of the following function:
 *  f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)^beta*
 *          f(0.5*(br-a)*x+0.5*(br+a))
 */
_10:
    hlgth = 0.5 * (br-bl);
    centr = 0.5 * (br+bl);
    fix = b-centr;
    fval[0]  = 0.5 * f(hlgth+centr, user_data)*pow(fix-hlgth,beta);
    fval[12] = f(centr, user_data) * pow(fix,beta);
    fval[24] = 0.5 * f(centr-hlgth, user_data)*pow(fix+hlgth,beta);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] = f(u+centr, user_data) * pow(fix-u,beta);
        fval[isym] = f(centr-u, user_data) * pow(fix+u,beta);
    }
    factor = pow(hlgth,alfa+1.0);
    result = 0.0;
    *abserr = 0.0;
    res12 = 0.0;
    res24 = 0.0;
    if (wgtfunc > 2) goto _70;
    dqcheb(x,fval,cheb12,cheb24);

/*  wgtfunc = 1  (or 2) */
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * ri[i]);
        res24 += (cheb24[i] * ri[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * ri[i]);
    }
    if (wgtfunc == 1) goto _130;

/*  wgtfunc = 2 */
    dc = log(br-bl);
    result = res24 * dc;
    res12 = 0.0;
    res24 = 0.0;
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rg[i]);
        res24 += (cheb24[i] * rg[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rg[i]);
    }
    goto _130;

/*  Compute the Chebyshev series expansion of the following function:
 *      f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
 */
_70:
    fval[0] *= log(fix-hlgth);
    fval[12] *= log(fix);
    fval[24] *= log(fix+hlgth);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] *= log(fix-u);
        fval[isym] *= log(fix+u);
    }
    dqcheb(x,fval,cheb12,cheb24);

/*  wgtfunc = 3  (or 4) */
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * ri[i]);
        res24 += (cheb24[i] * ri[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * ri[i]);
    }
    if (wgtfunc == 3) goto _130;

/*  wgtfunc = 4 */
    dc = log(br-bl);
    result = res24 * dc;
    *abserr = fabs((res24-res12)*dc);
    res12 = 0.0;
    res24 = 0.0;
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rg[i]);
        res24 += (cheb24[i] * rg[i]);
    }
    for (i=0;i<13;i++) {
        res24 += (cheb24[i] * rg[i]);
    }
_130:
    result = (result + res24) * factor;
    *abserr = (*abserr + fabs(res24-res12)) * factor;
    goto _270;

/*  This part is executed only if b = br
 *
 *  Compute the Chebyshev series expansion of the following function:
 *
 *  f2 = (0.5 *(b+bl-a-a)+0.5*(b-bl)*x)^alfa *
 *      f(0.5*(b-bl)*x+0.5*(b+bl))
 */
_140:
    hlgth = 0.5 * (b-bl);
    centr = 0.5 * (br+bl);
    fix = centr-a;
    fval[0]  = 0.5 * f(hlgth+centr, user_data) * pow(fix+hlgth,alfa);
    fval[12] = f(centr, user_data) * pow(fix,alfa);
    fval[24] = 0.5 * f(centr-hlgth, user_data) * pow(fix-hlgth,alfa);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] = f(u+centr, user_data) * pow(fix+u,alfa);
        fval[isym] = f(centr-u, user_data) * pow(fix-u,alfa);
    }
    factor = pow(hlgth,beta+1.0);
    result = 0.0;
    *abserr = 0.0;
    res12 = 0.0;
    res24 = 0.0;
    if ((wgtfunc == 2) || (wgtfunc == 4)) goto _200;

/*  wgtfunc = 1  (or 3)  */
    dqcheb(x,fval,cheb12,cheb24);
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rj[i]);
        res24 += (cheb24[i] * rj[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rj[i]);
    }
    if (wgtfunc == 1) goto _260;

/*  wgtfunc = 3  */
    dc = log(br-bl);
    result = res24 * dc;
    *abserr = fabs((res24 - res12) * dc);
    res12 = 0.0;
    res24 = 0.0;
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rh[i]);
        res24 += (cheb24[i] * rh[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rh[i]);
    }
_190:
    goto _260;

/*  Compute the Chebyshev series expansion of the following function:
 *
 *      f3 = f2 * log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
 */
_200:
    fval[0] *= log(hlgth+fix);
    fval[12] *= log(fix);
    fval[24] *= log(fix-hlgth);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] *= log(u+fix);
        fval[isym] *= log(fix-u);
    }
    dqcheb(x,fval,cheb12,cheb24);

/*  wgtfunc = 2  (or 4)  */
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rj[i]);
        res24 += (cheb24[i] * rj[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rj[i]);
    }
    if (wgtfunc == 2) goto _260;
    dc = log(br-bl);
    result = res24 * dc;
    *abserr = fabs((res24-res12) * dc);
    res12 = 0.0;
    res24 = 0.0;

/*  wgtfunc == 4  */
    for (i=0;i<13;i++) {
        res12 += (cheb12[i] * rh[i]);
        res24 += (cheb24[i] * rh[i]);
    }
    for (i=13;i<25;i++) {
        res24 += (cheb24[i] * rh[i]);
    }
_260:
    result = (result + res24)* factor;
    *abserr = (*abserr + fabs(res24-res12))*factor;
_270:
     return result;
}

