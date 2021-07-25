#include <stdio.h>
#include <math.h>
#include "cquadpack.h"

void dqmomo(double alpha,double beta,double ri[],double rj[],
    double rg[], double rh[],int wgtfunc)
{
    double alfp1,alfp2,an,anm1,betp1,betp2,ralf,rbet;
    int i,im1;

    alfp1 = alpha + 1.0;
    betp1 = beta + 1.0;
    alfp2 = alpha + 2.0;
    betp2 = beta + 2.0;
    ralf = pow(2.0,alfp1);
    rbet = pow(2.0,betp1);

/* Compute ri, rj using a forward recurrence relation. */
    ri[0] = ralf / alfp1;
    rj[0] = rbet / betp1;
    ri[1] = ri[0] * alpha / alfp2;
    rj[1] = rj[0] * beta / betp2;
    an = 2.0;
    anm1 = 1.0;
    for(i = 2;i < 25; i++) {
        ri[i] = -(ralf + an * (an - alfp2) * ri[i - 1]) /
            (anm1 * (an + alfp1));
        rj[i] = -(rbet + an * (an - betp2) * rj[i - 1]) /
            (anm1 * (an + betp1));
        anm1 = an;
        an += 1.0;
    }
    if (wgtfunc == 1)
        goto _70;
    if (wgtfunc == 3)
        goto _40;

/* Compute rg using a forward recurrence formula. */
    rg[0] = -ri[0] / alfp1;
    rg[1] = -(ralf + ralf) / (alfp2 * alfp2) - rg[0];
    an = 2.0;
    im1 = 1;    /* FORTRAN uses im1 = 2 */
    for (i = 2; i < 25; i++) {
        rg[i] = -(an * (an - alfp2) * rg[im1] - an * ri[im1] +
            anm1 * ri[i]) / (anm1 * (an + alfp1));
        anm1 = an;
        an += 1.0;
        im1 = i;
    }
    if (wgtfunc == 2)
        goto _70;

/* Compute rh using a forward recurrence relation. */
_40:
    rh[0] = -rj[0] / betp1;
    rh[1] = -(rbet + rbet) / (betp2 * betp2) - rh[0];
    an = 2.0;
    anm1 = 1.0;
    im1 = 1;    /* FORTRAN uses im1 = 2 */
    for (i = 2; i < 25; i++) {
        rj[i] = -(an * (an - betp2) * rh[im1] - an * rj    [im1] +
            anm1 * rj[i]) / ( anm1 * (an + betp1));
        anm1 = an;
        an += 1.0;
        im1 = i;
    }
    for (i = 1; i < 25; i += 2) {
        rh[i] = -rh[i];
    }
_70:
    for (i = 1; i < 25; i += 2) {
        rj[i] = -rj[i];
    }
}
