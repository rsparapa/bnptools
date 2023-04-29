#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
//extern SEXP cpsambrt_predict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpsambrt_predict2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP cprnft(SEXP, SEXP, SEXP, SEXP);
extern SEXP cprnft2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//extern SEXP cnft(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cnft2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
//    {"cpsambrt_predict", (DL_FUNC) &cpsambrt_predict,  7},
    {"cpsambrt_predict2",(DL_FUNC) &cpsambrt_predict2, 9},
//    {"cprnft",           (DL_FUNC) &cprnft,            4},
    {"cprnft2",          (DL_FUNC) &cprnft2,           6},
//    {"cnft",             (DL_FUNC) &cnft,             32},
    {"cnft2",            (DL_FUNC) &cnft2,            33},
    {NULL, NULL, 0}
};

void R_init_nftbart(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
