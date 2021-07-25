#include <float.h>
#include <math.h>
#include "cquadpack.h"

double G_K51(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
/* Gauss-Kronrod abscissae and weights for 51 - 25 rules */
    static double XGK51[26] = {
        0.99926210499260983419,
        0.99555696979049809791,
        0.98803579453407724764,
        0.97666392145951751150,
        0.96161498642584251242,
        0.94297457122897433941,
        0.92074711528170156175,
        0.89499199787827536885,
        0.86584706529327559545,
        0.83344262876083400142,
        0.79787379799850005941,
        0.75925926303735763058,
        0.71776640681308438819,
        0.67356636847346836449,
        0.62681009901031741279,
        0.57766293024122296772,
        0.52632528433471918260,
        0.47300273144571496052,
        0.41788538219303774885,
        0.36117230580938783774,
        0.30308953893110783017,
        0.24386688372098843205,
        0.18371893942104889202,
        0.12286469261071039639,
        0.06154448300568507889,
        0.00000000000000000000};
    static double WGK51[26] = {
        0.00198738389233031593,
        0.00556193213535671376,
        0.00947397338617415161,
        0.01323622919557167481,
        0.01684781770912829823,
        0.02043537114588283546,
        0.02400994560695321622,
        0.02747531758785173780,
        0.03079230016738748889,
        0.03400213027432933784,
        0.03711627148341554356,
        0.04008382550403238207,
        0.04287284502017004948,
        0.04550291304992178891,
        0.04798253713883671391,
        0.05027767908071567196,
        0.05236288580640747586,
        0.05425112988854549014,
        0.05595081122041231731,
        0.05743711636156783285,
        0.05868968002239420796,
        0.05972034032417405998,
        0.06053945537604586295,
        0.06112850971705304831,
        0.06147118987142531666,
        0.06158081806783293508};
    static double WG25[13] = {
        0.01139379850102628795,
        0.02635498661503213726,
        0.04093915670130631266,
        0.05490469597583519193,
        0.06803833381235691721,
        0.08014070033500101801,
        0.09102826198296364981,
        0.10053594906705064420,
        0.10851962447426365312,
        0.11485825914571164834,
        0.11945576353578477223,
        0.12224244299031004169,
        0.12317605372671545120};

    double fv1[25],fv2[25];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc=(*f)(centr, user_data);
    resg = fc * WG25[12];
    resk = fc * WGK51[25];
    *resabs = fabs(resk);
    for (j = 0; j < 12; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK51[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG25[j] * fsum;
        resk += WGK51[jtw] * fsum;
        *resabs = *resabs + WGK51[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 13; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK51[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK51[jtwm1] * fsum;
        *resabs = (*resabs) + WGK51[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK51[25] * fabs(fc - reskh);
    for (j = 0; j < 25; j++ )
        *resasc = (*resasc) + WGK51[j] * (fabs(fv1[j] - reskh) +
            fabs(fv2[j] - reskh));
    result = resk * hlgth;
    *resabs = (*resabs) * dhlgth;
    *resasc = (*resasc) * dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if ((*resasc != 0.0) && (*abserr != 0.0))
        *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
    if (*resabs > uflow/(50.0 * epmach))
        *abserr = max(epmach * 50.0 * (*resabs),(*abserr));
    return result;
}
