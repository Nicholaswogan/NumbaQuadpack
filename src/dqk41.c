#include <float.h>
#include <math.h>
#include "cquadpack.h"

double G_K41(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
/* Gauss-Kronrod abscissae and weights for 41 - 20 rules */
    static double XGK41[21] = {
        0.99885903158827766384,
        0.99312859918509492479,
        0.98150787745025025919,
        0.96397192727791379127,
        0.94082263383175475352,
        0.91223442825132590587,
        0.87827681125228197608,
        0.83911697182221882339,
        0.79504142883755119835,
        0.74633190646015079261,
        0.69323765633475138481,
        0.63605368072651502545,
        0.57514044681971031534,
        0.51086700195082709800,
        0.44359317523872510320,
        0.37370608871541956067,
        0.30162786811491300432,
        0.22778585114164507808,
        0.15260546524092267551,
        0.07652652113349733375,
        0.00000000000000000000};
    static double WGK41[21] = {
        0.00307358371852053150,
        0.00860026985564294220,
        0.01462616925697125298,
        0.02038837346126652360,
        0.02588213360495115883,
        0.03128730677703279896,
        0.03660016975820079803,
        0.04166887332797368626,
        0.04643482186749767472,
        0.05094457392372869193,
        0.05519510534828599474,
        0.05911140088063957237,
        0.06265323755478116803,
        0.06583459713361842211,
        0.06864867292852161935,
        0.07105442355344406831,
        0.07303069033278666750,
        0.07458287540049918899,
        0.07570449768455667466,
        0.07637786767208073671,
        0.07660071191799965645};
    static double WG20[10] = {
        0.01761400713915211831,
        0.04060142980038694133,
        0.06267204833410906357,
        0.08327674157670474872,
        0.10193011981724043504,
        0.11819453196151841731,
        0.13168863844917662690,
        0.14209610931838205133,
        0.14917298647260374679,
        0.15275338713072585070};
    double fv1[20],fv2[20];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    resg = 0.0;
    fc = (*f)(centr, user_data);
    resk = fc * WGK41[20];
    *resabs = fabs(resk);
    for (j = 0; j < 10; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK41[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG20[j] * fsum;
        resk += WGK41[jtw] * fsum;
        *resabs = *resabs + WGK41[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 10; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK41[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK41[jtwm1] * fsum;
        *resabs = (*resabs) + WGK41[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK41[20] * fabs(fc - reskh);
    for (j = 0; j < 20; j++ )
        *resasc = (*resasc) + WGK41[j] * (fabs(fv1[j] - reskh) +
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
