#include "cquadpack.h"

double G_K31(dq_function_type f,double a,double b,double *abserr,
    double *resabs,double *resasc, void* user_data)
{
    static double XGK31[16] = {
        0.99800229869339706029,
        0.98799251802048542849,
        0.96773907567913913426,
        0.93727339240070590431,
        0.89726453234408190088,
        0.84820658341042721620,
        0.79041850144246593297,
        0.72441773136017004742,
        0.65099674129741697053,
        0.57097217260853884754,
        0.48508186364023968069,
        0.39415134707756336990,
        0.29918000715316881217,
        0.20119409399743452230,
        0.10114206691871749903,
        0.00000000000000000000};
    static double WGK31[16] = {
        0.00537747987292334899,
        0.01500794732931612254,
        0.02546084732671532019,
        0.03534636079137584622,
        0.04458975132476487661,
        0.05348152469092808727,
        0.06200956780067064029,
        0.06985412131872825871,
        0.07684968075772037889,
        0.08308050282313302104,
        0.08856444305621177065,
        0.09312659817082532123,
        0.09664272698362367851,
        0.09917359872179195933,
        0.10076984552387559504,
        0.10133000701479154902};
    static double WG15[8] = {
        0.03075324199611726835,
        0.07036604748810812471,
        0.10715922046717193501,
        0.13957067792615431445,
        0.16626920581699393355,
        0.18616100001556221103,
        0.19843148532711157646,
        0.20257824192556127288};

    double fv1[15],fv2[15];
    double absc,centr,dhlgth;
    double fc,fsum,fval1,fval2,hlgth;
    double resg,resk,reskh,result;
    int j,jtw,jtwm1;

    centr = 0.5 * (a + b);
    hlgth = 0.5 * (b - a);
    dhlgth = fabs(hlgth);

    fc=(*f)(centr, user_data);
    resg = fc * WG15[7];
    resk = fc * WGK31[15];
    *resabs = fabs(resk);
    for (j = 0; j < 7; j++) {
        jtw = 2 * j + 1;
        absc = hlgth * XGK31[jtw];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1 + fval2;
        resg += WG15[j] * fsum;
        resk += WGK31[jtw] * fsum;
        *resabs = *resabs + WGK31[jtw] * (fabs(fval1) + fabs(fval2));
    }
    for (j = 0; j < 8; j++) {
        jtwm1 = j * 2;
        absc = hlgth * XGK31[jtwm1];
        fval1 = (*f)(centr-absc, user_data);
        fval2 = (*f)(centr+absc, user_data);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1 + fval2;
        resk = resk + WGK31[jtwm1] * fsum;
        *resabs = (*resabs) + WGK31[jtwm1] * (fabs(fval1) + fabs(fval2));
    }
    reskh = resk * 0.5;
    *resasc = WGK31[15] * fabs(fc - reskh);
    for (j = 0; j < 15; j++ )
        *resasc = (*resasc) + WGK31[j] * (fabs(fv1[j] - reskh) +
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
