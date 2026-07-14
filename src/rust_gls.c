#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <stddef.h>

extern int genomicsem_gls_batch(
    const double *betas,
    const double *ses,
    size_t n,
    size_t traits,
    const double *loadings,
    size_t factors,
    const double *corr,
    const double *q_corr,
    size_t requested_threads,
    double *beta_out,
    double *se_out,
    double *q_out,
    int *status_out);

static void matrix_dimensions(SEXP matrix, const char *name, int *rows, int *cols) {
    if (!isReal(matrix) || !isMatrix(matrix)) {
        error("'%s' must be a double matrix", name);
    }
    SEXP dimensions = getAttrib(matrix, R_DimSymbol);
    *rows = INTEGER(dimensions)[0];
    *cols = INTEGER(dimensions)[1];
}

SEXP C_genomicsem_gls_batch(
    SEXP betas,
    SEXP ses,
    SEXP loadings,
    SEXP corr,
    SEXP q_corr,
    SEXP threads) {
    int n, traits, se_n, se_traits, loading_traits, factors, corr_rows, corr_cols;
    int q_corr_rows, q_corr_cols;
    matrix_dimensions(betas, "betas", &n, &traits);
    matrix_dimensions(ses, "ses", &se_n, &se_traits);
    matrix_dimensions(loadings, "loadings", &loading_traits, &factors);
    matrix_dimensions(corr, "corr", &corr_rows, &corr_cols);
    matrix_dimensions(q_corr, "q_corr", &q_corr_rows, &q_corr_cols);

    if (n < 1 || traits < 1 || factors < 1) {
        error("betas, ses, and loadings must have non-zero dimensions");
    }
    if (se_n != n || se_traits != traits) {
        error("'ses' must have the same dimensions as 'betas'");
    }
    if (loading_traits != traits) {
        error("nrow(loadings) must equal ncol(betas)");
    }
    if (corr_rows != traits || corr_cols != traits) {
        error("'corr' must be a square matrix matching the number of traits");
    }
    if (q_corr_rows != traits || q_corr_cols != traits) {
        error("'q_corr' must be a square matrix matching the number of traits");
    }
    if ((TYPEOF(threads) != INTSXP && TYPEOF(threads) != REALSXP) || XLENGTH(threads) != 1) {
        error("'threads' must be a numeric scalar");
    }

    int thread_count = asInteger(threads);
    if (thread_count == NA_INTEGER || thread_count < 1) {
        error("'threads' must be at least 1");
    }

    SEXP beta_out = PROTECT(allocMatrix(REALSXP, n, factors));
    SEXP se_out = PROTECT(allocMatrix(REALSXP, n, factors));
    SEXP q_out = PROTECT(allocVector(REALSXP, n));
    SEXP status_out = PROTECT(allocVector(INTSXP, n));

    R_CheckUserInterrupt();
    int result = genomicsem_gls_batch(
        REAL(betas),
        REAL(ses),
        (size_t)n,
        (size_t)traits,
        REAL(loadings),
        (size_t)factors,
        REAL(corr),
        REAL(q_corr),
        (size_t)thread_count,
        REAL(beta_out),
        REAL(se_out),
        REAL(q_out),
        INTEGER(status_out));
    R_CheckUserInterrupt();

    if (result != 0) {
        UNPROTECT(4);
        if (result == 1) {
            error("the intercept correlation matrix contains non-finite values");
        } else if (result == 2) {
            error("the intercept correlation matrix is not symmetric");
        } else if (result == 3) {
            error("the intercept correlation matrix used for Q is singular");
        } else {
            error("the Rust analytic kernel failed unexpectedly");
        }
    }

    SEXP output = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(output, 0, beta_out);
    SET_VECTOR_ELT(output, 1, se_out);
    SET_VECTOR_ELT(output, 2, q_out);
    SET_VECTOR_ELT(output, 3, status_out);

    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("beta"));
    SET_STRING_ELT(names, 1, mkChar("se"));
    SET_STRING_ELT(names, 2, mkChar("q"));
    SET_STRING_ELT(names, 3, mkChar("status"));
    setAttrib(output, R_NamesSymbol, names);

    UNPROTECT(6);
    return output;
}
