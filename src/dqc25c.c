#include "cquadpack.h"

/*  DQC25C - Integration rules for the computation of Cauchy
 *          principal value integrals.
 *
 *  PARAMETERS:
 *
 *      f() - double precision function defining the integrand
 *
 *      a   - left end point of the integration interval
 *
 *      b   - right end point of the integration interval
 *
 *      c   - parameter in the weight function
 *
 *      abserr  - estimate of the modulus of the absolute error
 *
 *      krul    - key which is decreased by 1 if the 15-point
 *                  Gauss-Kronrod scheme is used
 *
 *      neval   - number of function evaluations
 */
double dqc25c(dq_function_type f,double a,double b,double c,double *abserr,
        int *krul, int *neval, void* user_data)
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
    double ak22,amom0,amom1,amom2,cc,centr;
    double cheb12[13],cheb24[25],fval[25];
    double hlgth,resabs,resasc,res12,res24,u,result;
    int i,isym,k;
    int unitialized_value = 0xCCCCCCCC;
    int kp = unitialized_value;
    double p2 = unitialized_value, p3 = unitialized_value, p4 = unitialized_value;

    cc = (2.0 * c - b - a) / (b - a);
    if (fabs(cc) < 1.1) goto _10;

/*  Apply the 15-point Gauss-Kronrod scheme.    */
    (*krul)--;
    result = G_K15W(f,dqwgtc,c,p2,p3,p4,kp,a,b,abserr,&resabs,&resasc, user_data);
    *neval = 15;
    if (resasc == *abserr) (*krul)++;
    goto _50;

/*  Use the generalized Clenshaw-Curtis method. */
_10:
    hlgth = 0.5 * (b - a);
    centr = 0.5 * (b + a);
    *neval = 25;
    fval[0] = 0.5 * f(hlgth+centr, user_data);
    fval[12] = f(centr, user_data);
    fval[24] = 0.5 * f(centr-hlgth, user_data);
    for (i=1;i<12;i++) {
        u = hlgth * x[i-1];
        isym = 24 - i;
        fval[i] = f(u+centr, user_data);
        fval[isym] = f(centr-u, user_data);
    }

/*  Compute the Chebyshev series expansion. */
    dqcheb(x,fval,cheb12,cheb24);

/*  The modified Chebyshev moments are computed by forward
 *  recursion, using amom0 and amom1 as starting values.
 */
    amom0 = log(fabs((1.0-cc)/(1.0+cc)));
    amom1 = 2.0 + cc * amom0;
    res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
    res24 = cheb24[0] * amom0 + cheb24[1] * amom1;
    for (k = 2;k < 13;k++) {
        amom2 = 2.0 * cc * amom1 - amom0;
        ak22 = (k-1) * (k-1);
        if ((k/2)*2 != k) amom2 -= (4.0 / (ak22 - 1.0));
        res12 += (cheb12[k] * amom2);
        res24 += (cheb24[k] * amom2);
        amom0 = amom1;
        amom1 = amom2;
    }
    for (k = 13;k < 25;k++) {
        amom2 = 2.0 * cc * amom1 - amom0;
        ak22 = (k-1) * (k-1);
        if ((k/2)*2 != k) amom2 -= (4.0 /(ak22 - 1.0));
        res24 += (cheb24[k] * amom2);
        amom0 = amom1;
        amom1 = amom2;
    }
    result = res24;
    *abserr = fabs(res24-res12);
_50:
    return result;
}
