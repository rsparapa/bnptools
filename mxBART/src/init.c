#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
/* extern SEXP cmbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); */
extern SEXP cpwbart(SEXP, SEXP, SEXP);
extern SEXP cgbmm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/*extern SEXP cspbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);*/
extern SEXP mc_cores_openmp();
extern SEXP crtnorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP crtgamma(SEXP, SEXP, SEXP, SEXP);
extern SEXP cdraw_lambda_i(SEXP, SEXP, SEXP, SEXP);
/*extern SEXP cdpmbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cdpmwbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cdpgbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);*/

static const R_CallMethodDef CallEntries[] = {
/*  {"cmbart",  (DL_FUNC) &cmbart,  29},*/
    {"cpwbart", (DL_FUNC) &cpwbart,  3},
    {"cgbmm",   (DL_FUNC) &cgbmm,   39}, 
/*  {"cspbart",  (DL_FUNC) &cspbart,  30}, */
    {"mc_cores_openmp",(DL_FUNC) &mc_cores_openmp,0},
    {"crtnorm", (DL_FUNC) &crtnorm,  4},
    {"crtgamma",(DL_FUNC) &crtgamma, 4},
    {"cdraw_lambda_i", (DL_FUNC) &cdraw_lambda_i, 4},
/*  {"cdpgbart",(DL_FUNC) &cdpgbart,35},
    {"cdpmbart",(DL_FUNC) &cdpmbart,25},
    {"cdpmwbart",(DL_FUNC) &cdpmwbart,27}, */
    {NULL, NULL, 0}
};

void R_init_BART(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
