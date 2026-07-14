#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP C_genomicsem_gls_batch(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_genomicsem_gls_columns(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_genomicsem_gls_results_columns(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"genomicsem_gls_batch", (DL_FUNC) &C_genomicsem_gls_batch, 6},
    {"genomicsem_gls_columns", (DL_FUNC) &C_genomicsem_gls_columns, 8},
    {"genomicsem_gls_results_columns", (DL_FUNC) &C_genomicsem_gls_results_columns, 9},
    {NULL, NULL, 0}
};

void attribute_visible R_init_GenomicSEM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
