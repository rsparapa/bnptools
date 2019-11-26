#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP call_mutau(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_NoGa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_rcat(SEXP);
extern SEXP call_rgamma(SEXP, SEXP, SEXP);
extern SEXP call_rlgamma(SEXP, SEXP);
extern SEXP call_rtnorm(SEXP, SEXP, SEXP, SEXP);
extern SEXP call_quantile(SEXP, SEXP);
extern SEXP call_floor(SEXP);
extern SEXP call_ceiling(SEXP);
extern SEXP call_sort(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"call_mutau",       (DL_FUNC) &call_mutau,       7},
  {"call_NoGa",        (DL_FUNC) &call_NoGa,        7},
  {"call_rcat",        (DL_FUNC) &call_rcat,        1},
  {"call_rgamma",      (DL_FUNC) &call_rgamma,      3},
  {"call_rlgamma",     (DL_FUNC) &call_rlgamma,     2},
  {"call_rtnorm",      (DL_FUNC) &call_rtnorm,      4},
  {"call_quantile",    (DL_FUNC) &call_quantile,    2},
  {"call_floor",       (DL_FUNC) &call_floor,       1},
  {"call_ceiling",     (DL_FUNC) &call_ceiling,     1},
  {"call_sort",        (DL_FUNC) &call_sort,        1},
  {NULL, NULL, 0}
};

void R_init_DPM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
