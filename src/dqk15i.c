#include "cquadpack.h"

double G_K15I(dq_function_type f, double boun, int inf, double a, double b,
    double *abserr,double *resabs,double *resasc, void* user_data)
{
    static double XGK15[8] = {
        0.99145537112081263921,
        0.94910791234275852453,
        0.86486442335976907279,
        0.74153118559939443986,
        0.58608723546769113029,
        0.40584515137739716691,
        0.20778495500789846760,
        0.00000000000000000000};
    static double WGK15[8] = {
        0.02293532201052922496,
        0.06309209262997855329,
        0.10479001032225018384,
        0.14065325971552591875,
        0.16900472663926790283,
        0.19035057806478540991,
        0.20443294007529889241,
        0.20948214108472782801};
    static double WG7[4] = {
        0.12948496616886969327,
        0.27970539148927666790,
        0.38183005050511894495,
        0.41795918367346938776};
    double fv1[8],fv2[8];
    double absc,absc1,absc2,centr,dinf;
    double fc,fsum,fval1,fval2,hlgth,resg,resk;
    double reskh,result,tabsc1,tabsc2;
    int j;

    dinf = min((double)(1.0),(double)inf);
    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    tabsc1 = boun + dinf * (1.0 - centr)/centr;
    fval1 = (*f)(tabsc1, user_data);
    if (inf == 2)
        fval1 += (*f)(-tabsc1, user_data);
    fc=(fval1/centr)/centr;
    resg = fc * WG7[3];
    resk = fc * WGK15[7];
    *resabs = fabs(resk);
    for (j = 0; j < 7; j++) {
        absc = hlgth * XGK15[j];
        absc1 = centr - absc;
        absc2 = centr + absc;
        tabsc1 = boun + dinf * (1.0 - absc1)/absc1;
        tabsc2 = boun + dinf * (1.0 - absc2)/absc2;
        fval1 = (*f)(tabsc1, user_data);
        fval2 = (*f)(tabsc2, user_data);
        if (inf == 2) {
            fval1 += (*f)(-tabsc1, user_data);
            fval2 += (*f)(-tabsc2, user_data);
        }
        fval1 = (fval1/absc1)/absc1;
        fval2 = (fval2/absc2)/absc2;
        fv1[j] = fval1;
        fv2[j] = fval2;
        fsum = fval1 + fval2;
        if (j & 1) resg += WG7[j/2] * fsum; /* odd 'j's are truncated */
        resk += WGK15[j] * fsum;
        *resabs = (*resabs) + WGK15[j] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK15[7] * fabs(fc - reskh);
    for (j = 0; j < 7; j++ )
        *resasc = (*resasc) + WGK15[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * hlgth;
    *resasc = (*resasc) * hlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}
