#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP sparseHessianFD_get_colors(SEXP, SEXP, SEXP);
extern SEXP sparseHessianFD_subst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"sparseHessianFD_get_colors", (DL_FUNC) &sparseHessianFD_get_colors, 3},
    {"sparseHessianFD_subst",      (DL_FUNC) &sparseHessianFD_subst,      7},
    {NULL, NULL, 0}
};

void R_init_sparseHessianFD(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
