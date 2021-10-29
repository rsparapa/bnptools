#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cdotree(SEXP, SEXP, SEXP, SEXP);
extern SEXP cmbrt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpsambrt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpsambrt_predict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpsambrt_vartivity(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpsambrt_save(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpsambrt_Rexport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpsambrt_Rimport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cphbart(SEXP, SEXP, SEXP);
/* extern SEXP cnft(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); */

static const R_CallMethodDef CallEntries[] = {
    {"cdotree",          (DL_FUNC) &cdotree,           4},
    {"cmbrt",            (DL_FUNC) &cmbrt,            17},
    {"cpsambrt",         (DL_FUNC) &cpsambrt,         32},
    {"cpsambrt_predict", (DL_FUNC) &cpsambrt_predict,  8},
   {"cpsambrt_vartivity",(DL_FUNC) &cpsambrt_vartivity,6},
    {"cpsambrt_save",    (DL_FUNC) &cpsambrt_save,     7},
    {"cpsambrt_Rexport", (DL_FUNC) &cpsambrt_Rexport,  6},
    {"cpsambrt_Rimport", (DL_FUNC) &cpsambrt_Rimport,  6},
    {"cphbart",          (DL_FUNC) &cphbart,           3},
    //{"cnft",             (DL_FUNC) &cnft,             36},
    {NULL, NULL, 0}
};

void R_init_hbart(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
