#include "cquadpack.h"


double dqext(int *n,double epstab[],double *abserr,
    double res3la[],int *nres)
{
    double delta1,delta2,delta3,epsinf;
    double error,err1,err2,err3,e0,e1,e1abs,e2,e3;
    double res,result,ss,tol1,tol2,tol3;
    int NN,i,ib,ib2,ie,indx,k1,k2,k3,limexp,newelm,num;

    (*nres)++;
    NN = *n;
    NN++;   /* make NN a FORTRAN array index */
    *abserr = oflow;
    result = epstab[*n];
    if (NN < 3) goto _100;        /* N < 3 */
    limexp = 50;            /* limexp = 50 */
    epstab[*n+2] = epstab[*n];
    newelm = (*n)/2;      /* (n-1)/2 */
    epstab[*n] = oflow;
    num = NN;
    k1 = NN;
    for (i = 1; i <= newelm; i++) {
        k2 = k1 - 1;
        k3 = k1 - 2;
        res = epstab[k1+1];
        e0 = epstab[k3-1];
        e1 = epstab[k2-1];
        e2 = res;
        e1abs = fabs(e1);
        delta2 = e2 - e1;
        err2 = fabs(delta2);
        tol2 = max(fabs(e2),e1abs) * epmach;
        delta3 = e1 - e0;
        err3 = fabs(delta3);
        tol3 = max(e1abs,fabs(e0)) * epmach;
        if ((err2 > tol2) || (err3 > tol3)) goto _10;
        result = res;
        *abserr = err2 + err3;
        goto _100;
_10:
        e3 = epstab[k1-1];
        epstab[k1-1] = e1;
        delta1 = e1 - e3;
        err1 = fabs(delta1);
        tol1 = max(e1abs,fabs(e3)) * epmach;
        if ((err1 <= tol1) || (err2 <= tol2) || (err3 <= tol3)) goto _20;
        ss = 1.0/delta1 + 1.0/delta2 - 1.0/delta3;
        epsinf = fabs(ss*e1);
        if (epsinf > 1.0e-4) goto _30;
_20:
        NN = i + i - 1;
        goto _50;
_30:
        res = e1 + 1.0 / ss;
        epstab[k1-1] = res;
        k1 -= 2;
        error = err2 + fabs(res - e2) + err3;
        if (error > (*abserr)) goto _40;
        *abserr = error;
        result = res;
_40:
        ;
    }
_50:
    if (NN == limexp) NN = 2 * (limexp/2) - 1;
    ib = 1;                        /* ib = 1 */
    if (((num/2) * 2 ) == num) ib = 2;        /* ib = 2 */
    ie = newelm + 1;
    for (i = 1;i <= ie; i++) {
        ib2 = ib + 2;
        epstab[ib-1] = epstab[ib2-1];
        ib = ib2;
    }
    if (num == NN) goto _80;
    indx = num - NN + 1;
    for (i = 1;i <= NN; i++) {
        epstab[i-1] = epstab[indx-1];
        indx++;
    }
_80:
    if (*nres > 3) goto _90;       /* nres >= 4 */
    res3la[(*nres)-1] = result;
    *abserr = oflow;
    goto _100;
_90:
    *abserr = fabs(result - res3la[2]) + fabs(result - res3la[1]) +
        fabs(result - res3la[0]);
    res3la[0] = res3la[1];
    res3la[1] = res3la[2];
    res3la[2] = result;
_100:
    *abserr = max(*abserr,5.0 * epmach * fabs(result));
    *n = NN - 1;
    return result;
    }


