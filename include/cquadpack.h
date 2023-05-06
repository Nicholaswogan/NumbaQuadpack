#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>

#include "cquadpack_export.h"

#define uflow     DBL_MIN
#define oflow     DBL_MAX
#define epmach     DBL_EPSILON
#define LIMIT     500
#define MAXP1     21
#ifdef M_PI
#define Pi      M_PI
#else
#define Pi      3.14159265358979323846
#endif
#define COSINE     1
#define SINE    2

#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif
#ifndef min
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif
#ifndef max
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef double(*dq_function_type)(double, void*);

/* Log Sum / Diff functions */
CQUADPACK_EXPORT double logsumexp(double a, double b);
CQUADPACK_EXPORT double logsubexp(double a, double b);
#define LOGDIFF(x, y) ((x) < (y) ? logsubexp(y, x) : logsubexp(x, y))

/* Integration routines */
/* Gauss-Kronrod for integration over finite range. */
CQUADPACK_EXPORT double G_K15(dq_function_type f,double a,double b,double *abserr,
    double *resabs, double *resasc, void* user_data);
CQUADPACK_EXPORT double G_K21(dq_function_type f, double a, double b, double *abserr,
    double *resabs, double *resasc, void* user_data);
CQUADPACK_EXPORT double G_K31(dq_function_type f, double a, double b, double *abserr,
    double *resabs, double *resasc, void* user_data);
CQUADPACK_EXPORT double G_K41(dq_function_type f, double a, double b, double *abserr,
    double *resabs, double *resasc, void* user_data);
CQUADPACK_EXPORT double G_K51(dq_function_type f, double a, double b, double *abserr,
    double *resabs, double *resasc, void* user_data);
CQUADPACK_EXPORT double G_K61(dq_function_type f, double a, double b, double *abserr,
    double *resabs, double *resasc, void* user_data);

/* Gauss-Kronrod for integration over infinite range. */
CQUADPACK_EXPORT double G_K15I(dq_function_type f, double boun, int inf, double a, double b,
    double *abserr,double *resabs, double *resasc, void* user_data);

/* Gauss-Kronrod for integration of weighted function. */
CQUADPACK_EXPORT double G_K15W(dq_function_type f, double w(), double p1, double p2, double p3,
    double p4,int kp,double a,double b,double *abserr,
    double *resabs, double *resasc, void* user_data);
CQUADPACK_EXPORT double dqext(int *n, double epstab [], double *abserr,
    double res3la[],int *nres);
CQUADPACK_EXPORT void dqsort(int limit, int last, int *maxerr, double *ermax,
    double elist[],int iord[],int *nrmax);
CQUADPACK_EXPORT double dqagi(dq_function_type f, double bound, int inf, double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT double dqags(dq_function_type f, double a, double b, double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT double dqagp(dq_function_type f, double a, double b, int npts2, double *points,
    double epsabs,double epsrel,double *abserr,int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT double dqng(dq_function_type f, double a, double b, double epsabs, double epsrel,
    double *abserr,int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT double dqag(dq_function_type f, double a, double b, double epsabs, double epsrel,
    int irule,double *abserr,int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT double dqage(dq_function_type f, double a, double b, double epsabs, double epsrel,
    int irule,double *abserr,int *neval,int *ier,int *last, void* user_data);
CQUADPACK_EXPORT double dqwgtc(double x, double c, double p2, double p3, double p4,
    int kp);
CQUADPACK_EXPORT double dqwgto(double x, double omega, double p2, double p3, double p4,
    int integr);
CQUADPACK_EXPORT double dqwgts(double x, double a, double b, double alpha, double beta,
    int integr);
CQUADPACK_EXPORT void dqcheb(double *x, double *fval, double *cheb12, double *cheb24);
CQUADPACK_EXPORT double dqc25o(dq_function_type f, double a, double b, double omega, int integr,
    int nrmom,int maxp1,int ksave,double *abserr,int *neval,
    double *resabs,double *resasc,int *momcom,double **chebmo, void* user_data);
CQUADPACK_EXPORT double dqfour(dq_function_type f, double a, double b, double omega, int integr,
    double epsabs,double epsrel,int icall,int maxp1,
    double *abserr,int *neval,int *ier,int *momcom,
    double **chebmo, void* user_data);
CQUADPACK_EXPORT double dqawfe(dq_function_type f, double a, double omega, int integr, double epsabs,
    int limlst,int maxp1,double *abserr,int *neval,int *ier,
    double *rslst,double *erlist,int *ierlst,double **chebmo, void* user_data);
CQUADPACK_EXPORT double dqawf(dq_function_type f, double a, double omega, int integr, double epsabs,
    double *abserr,int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT double dqawo(dq_function_type f, double a, double b, double omega, int integr, double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT double dqaws(dq_function_type f, double a, double b, double alfa, double beta, int wgtfunc,
    double epsabs,double epsrel,double *abserr,int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT double dqawse(dq_function_type f, double a, double b, double alfa, double beta,
    int wgtfunc,double epsabs,double epsrel,double *abserr,
    int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT void dqmomo(double alfa, double beta, double ri [], double rj [], double rg [],
    double rh[],int wgtfunc);
CQUADPACK_EXPORT double dqc25s(dq_function_type f, double a, double b, double bl, double br, double alfa,
    double beta,double ri[],double rj[],double rg[],double rh[],
    double *abserr,double *resasc,int wgtfunc,int *nev, void* user_data);
CQUADPACK_EXPORT double dqc25c(dq_function_type f, double a, double b, double c, double *abserr,
    int *krul,int *neval, void* user_data);
CQUADPACK_EXPORT double dqawc(dq_function_type f, double a, double b, double c, double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data);
CQUADPACK_EXPORT double dqawce(dq_function_type f, double a, double b, double c, double epsabs,
    double epsrel,double *abserr,int *neval,int *ier, void* user_data);

CQUADPACK_EXPORT double G_B15(dq_function_type f, double a, double b, double *abserr,
    double *resabs, double *resasc, void* user_data);

#ifdef __cplusplus
}
#endif /* __cplusplus */
