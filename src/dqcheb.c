#include "cquadpack.h"

void dqcheb(double x[],double fval[],double cheb12[],double cheb24[])
{
    double alam,alam1,alam2,part1,part2,part3;
    double v[12];
    int i,j;

/* Dimensions of input vectors are:
 *        x[11], fval[25], cheb12[13], cheb24[25]
 */
     for (i = 0; i < 12; i++) {
         j = 24 - i;
         v[i] = fval[i] - fval[j];
         fval[i] += fval[j];
     }
     alam1 = v[0] - v[8];
     alam2 = x[5] * (v[2] - v[6] - v[10]);
     cheb12[3] = alam1 + alam2;
     cheb12[9] = alam1 - alam2;
     alam1 = v[1] - v[7] - v[9];
     alam2 = v[3] - v[5] - v[11];
     alam = x[2] * alam1 + x[8] * alam2;
     cheb24[3] = cheb12[3] + alam;
     cheb24[21] = cheb12[3] - alam;
     alam = x[8] * alam1 - x[2] * alam2;
     cheb24[9] = cheb12[9] + alam;
     cheb24[15] = cheb12[9] - alam;
     part1 = x[3] * v[4];
     part2 = x[7] * v[8];
     part3 = x[5] * v[6];
     alam1 = v[0] + part1 + part2;
     alam2 = x[1] * v[2] + part3 + x[9] * v[10];
     cheb12[1] = alam1 + alam2;
     cheb12[11] = alam1 - alam2;
     alam = x[0] * v[1] + x[2] * v[3] + x[4] * v[5] +
         x[6] * v[7] + x[8] * v[9] + x[10] * v[11];
     cheb24[1] = cheb12[1] + alam;
     cheb24[23] = cheb12[1] - alam;
     alam = x[10] * v[1] - x[8] * v[3] + x[6] * v[5] -
         x[4] * v[7] + x[2] * v[9] - x[0] * v[11];
     cheb24[11] = cheb12[11] + alam;
     cheb24[13] = cheb12[11] - alam;
     alam1 = v[0] - part1 + part2;
     alam2 = x[9] * v[2] - part3 + x[1] * v[10];
     cheb12[5] = alam1 + alam2;
     cheb12[7] = alam1 - alam2;
     alam = x[4] * v[1] - x[8] * v[3] - x[0] * v[5] -
         x[10] * v[7] + x[2] * v[9] + x[6] * v[11];
     cheb24[5] = cheb12[5] + alam;
     cheb24[19] = cheb12[5] - alam;
     alam = x[6] * v[1] - x[2] * v[3] - x[10] * v[5] +
         x[0] * v[7] - x[8] * v[9] - x[4] * v[11];
     cheb24[7] = cheb12[7] + alam;
     cheb24[17] = cheb12[7] - alam;
     for (i = 0; i < 6; i++) {
         j = 12 - i;
         v[i] = fval[i] - fval[j];
         fval[i] += fval[j];
     }
     alam1 = v[0] + x[7] * v[4];
     alam2 = x[3] * v[2];
     cheb12[2] = alam1 + alam2;
     cheb12[10] = alam1 - alam2;
     cheb12[6] = v[0] - v[4];
     alam = x[1] * v[1] + x[5] * v[3] + x[9] * v[5];
     cheb24[2] = cheb12[2] + alam;
     cheb24[22] = cheb12[2] - alam;
     alam = x[5] * (v[1] - v[3] - v[5]);
     cheb24[6] = cheb12[6] + alam;
     cheb24[18] = cheb12[6] - alam;
     alam = x[9] * v[1] - x[5] * v[3] + x[1] * v[5];
     cheb24[10] = cheb12[10] + alam;
     cheb24[14] = cheb12[10] - alam;
     for (i = 0; i < 3; i++) {
         j = 6 - i;
         v[i] = fval[i] -fval[j];
         fval[i] += fval[j];
     }
     cheb12[4] = v[0] + x[7] * v[2];
     cheb12[8] = fval[0] - x[7] * fval[2];
     alam = x[3] * v[1];
     cheb24[4] = cheb12[4] + alam;
     cheb24[20] = cheb12[4] - alam;
     alam = x[7] * fval[1] - fval[3];
     cheb24[8] = cheb12[8] + alam;
     cheb24[16] = cheb12[8] - alam;
     cheb12[0] = fval[0] + fval[2];
     alam = fval[1] + fval[3];
     cheb24[0] = cheb12[0] + alam;
     cheb24[24] = cheb12[0] - alam;
     cheb12[12] = v[0] - v[2];
     cheb24[12] = cheb12[12];
      alam = 1.0 / 6.0;
      for (i = 1; i < 12; i++)
          cheb12[i] *= alam;
      alam *= 0.5;
      cheb12[0] *= alam;
      cheb12[12] *= alam;
      for (i = 1; i < 24; i ++)
          cheb24[i] *= alam;
      cheb24[0] *= (0.5 * alam);
      cheb24[24] *= (0.5 * alam);
}

