#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cpsambrt_predict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cprnft(SEXP, SEXP, SEXP);
extern SEXP cnft(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP cpsambrt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP cpsambrt_vartivity(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP cpsambrt_save(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP cpsambrt_Rexport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP cpsambrt_Rimport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP crtnorm(SEXP, SEXP, SEXP, SEXP);
//extern SEXP chft(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP CallDPMLIOdensity(SEXP, SEXP, SEXP);
//extern SEXP CallDPMLIOmutau(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP CallDPMlambda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cpsambrt_predict", (DL_FUNC) &cpsambrt_predict,  7},
    {"cprnft",           (DL_FUNC) &cprnft,            3},
    {"cnft",             (DL_FUNC) &cnft,             35},
//    {"cpsambrt",         (DL_FUNC) &cpsambrt,         32},
//   {"cpsambrt_vartivity",(DL_FUNC) &cpsambrt_vartivity,6},
//    {"cpsambrt_save",    (DL_FUNC) &cpsambrt_save,     7},
//    {"cpsambrt_Rexport", (DL_FUNC) &cpsambrt_Rexport,  6},
//    {"cpsambrt_Rimport", (DL_FUNC) &cpsambrt_Rimport,  6},
//    {"crtnorm",          (DL_FUNC) &crtnorm,           4},
//    {"chft",             (DL_FUNC) &chft,             35},
//  {"CallDPMLIOdensity",(DL_FUNC) &CallDPMLIOdensity, 3},
//    {"CallDPMLIOmutau",  (DL_FUNC) &CallDPMLIOmutau,   7},
//    {"CallDPMlambda",    (DL_FUNC) &CallDPMlambda,     7},
    {NULL, NULL, 0}
};

void R_init_nftbart(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
