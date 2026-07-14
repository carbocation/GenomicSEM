#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>

typedef double (*normal_p_value_fn)(double);
typedef double (*chi_square_p_value_fn)(double, double);

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

extern int genomicsem_gls_batch_columns(
    const double *const *betas,
    const double *const *ses,
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

extern int genomicsem_gls_batch_columns_results(
    const double *const *betas,
    const double *const *ses,
    size_t n,
    size_t traits,
    const double *loadings,
    size_t factors,
    const double *corr,
    const double *q_corr,
    size_t requested_threads,
    double q_df,
    normal_p_value_fn normal_p_value,
    chi_square_p_value_fn chi_square_p_value,
    double *results_out,
    int *status_out);

static void matrix_dimensions(SEXP matrix, const char *name, int *rows, int *cols) {
    if (!isReal(matrix) || !isMatrix(matrix)) {
        error("'%s' must be a double matrix", name);
    }
    SEXP dimensions = getAttrib(matrix, R_DimSymbol);
    *rows = INTEGER(dimensions)[0];
    *cols = INTEGER(dimensions)[1];
}

typedef struct {
    int traits;
    int factors;
    int threads;
    size_t row_count;
    const double **betas;
    const double **ses;
} column_inputs;

static column_inputs prepare_column_inputs(
    SEXP betas,
    SEXP ses,
    SEXP loadings,
    SEXP corr,
    SEXP q_corr,
    SEXP start,
    SEXP count,
    SEXP threads) {
    if (TYPEOF(betas) != VECSXP || TYPEOF(ses) != VECSXP) {
        error("'betas' and 'ses' must be lists of double columns");
    }
    R_xlen_t traits_x = XLENGTH(betas);
    if (traits_x < 1 || traits_x > INT_MAX || XLENGTH(ses) != traits_x) {
        error("'betas' and 'ses' must contain the same non-zero number of columns");
    }
    int traits = (int) traits_x;
    int loading_traits, factors, corr_rows, corr_cols, q_corr_rows, q_corr_cols;
    matrix_dimensions(loadings, "loadings", &loading_traits, &factors);
    matrix_dimensions(corr, "corr", &corr_rows, &corr_cols);
    matrix_dimensions(q_corr, "q_corr", &q_corr_rows, &q_corr_cols);
    if (loading_traits != traits || factors < 1) {
        error("nrow(loadings) must match the number of trait columns");
    }
    if (corr_rows != traits || corr_cols != traits) {
        error("'corr' must be a square matrix matching the number of traits");
    }
    if (q_corr_rows != traits || q_corr_cols != traits) {
        error("'q_corr' must be a square matrix matching the number of traits");
    }

    if ((TYPEOF(start) != INTSXP && TYPEOF(start) != REALSXP) || XLENGTH(start) != 1) {
        error("'start' must be a numeric scalar");
    }
    if ((TYPEOF(count) != INTSXP && TYPEOF(count) != REALSXP) || XLENGTH(count) != 1) {
        error("'count' must be a numeric scalar");
    }
    double start_value = asReal(start);
    double count_value = asReal(count);
    if (!R_FINITE(start_value) || start_value < 0 || start_value > INT_MAX ||
        start_value != (size_t) start_value) {
        error("'start' must be a non-negative whole number");
    }
    if (!R_FINITE(count_value) || count_value < 1 || count_value > INT_MAX ||
        count_value != (size_t) count_value) {
        error("'count' must be a positive whole number no larger than INT_MAX");
    }
    size_t start_index = (size_t) start_value;
    size_t row_count = (size_t) count_value;

    if ((TYPEOF(threads) != INTSXP && TYPEOF(threads) != REALSXP) || XLENGTH(threads) != 1) {
        error("'threads' must be a numeric scalar");
    }
    int thread_count = asInteger(threads);
    if (thread_count == NA_INTEGER || thread_count < 1) {
        error("'threads' must be at least 1");
    }

    const double **beta_pointers =
        (const double **) R_alloc((size_t) traits, sizeof(double *));
    const double **se_pointers =
        (const double **) R_alloc((size_t) traits, sizeof(double *));
    for (int trait = 0; trait < traits; ++trait) {
        SEXP beta_column = VECTOR_ELT(betas, trait);
        SEXP se_column = VECTOR_ELT(ses, trait);
        if (TYPEOF(beta_column) != REALSXP || TYPEOF(se_column) != REALSXP) {
            error("all beta and SE columns must be double vectors");
        }
        R_xlen_t beta_length = XLENGTH(beta_column);
        R_xlen_t se_length = XLENGTH(se_column);
        if (start_index > (size_t) beta_length ||
            row_count > (size_t) beta_length - start_index ||
            start_index > (size_t) se_length ||
            row_count > (size_t) se_length - start_index) {
            error("the requested row range exceeds a beta or SE column");
        }
        beta_pointers[trait] = REAL(beta_column) + start_index;
        se_pointers[trait] = REAL(se_column) + start_index;
    }

    column_inputs inputs = {
        traits,
        factors,
        thread_count,
        row_count,
        beta_pointers,
        se_pointers
    };
    return inputs;
}

static void check_kernel_result(int result) {
    if (result == 0) {
        return;
    }
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

/* These Rmath entry points are pure for the validated parameters used here and
   may therefore be invoked by Rust's worker threads without touching the R API. */
static double two_sided_normal_p(double z) {
    return 2.0 * pnorm5(-fabs(z), 0.0, 1.0, 1, 0);
}

static double upper_chi_square_p(double q, double df) {
    return pchisq(q, df, 0, 0);
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

SEXP C_genomicsem_gls_columns(
    SEXP betas,
    SEXP ses,
    SEXP loadings,
    SEXP corr,
    SEXP q_corr,
    SEXP start,
    SEXP count,
    SEXP threads) {
    column_inputs inputs = prepare_column_inputs(
        betas, ses, loadings, corr, q_corr, start, count, threads);
    int n = (int) inputs.row_count;
    SEXP beta_out = PROTECT(allocMatrix(REALSXP, n, inputs.factors));
    SEXP se_out = PROTECT(allocMatrix(REALSXP, n, inputs.factors));
    SEXP q_out = PROTECT(allocVector(REALSXP, n));
    SEXP status_out = PROTECT(allocVector(INTSXP, n));

    R_CheckUserInterrupt();
    int result = genomicsem_gls_batch_columns(
        inputs.betas,
        inputs.ses,
        inputs.row_count,
        (size_t) inputs.traits,
        REAL(loadings),
        (size_t) inputs.factors,
        REAL(corr),
        REAL(q_corr),
        (size_t) inputs.threads,
        REAL(beta_out),
        REAL(se_out),
        REAL(q_out),
        INTEGER(status_out));
    R_CheckUserInterrupt();

    check_kernel_result(result);

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

SEXP C_genomicsem_gls_results_columns(
    SEXP betas,
    SEXP ses,
    SEXP loadings,
    SEXP corr,
    SEXP q_corr,
    SEXP start,
    SEXP count,
    SEXP threads,
    SEXP q_df) {
    column_inputs inputs = prepare_column_inputs(
        betas, ses, loadings, corr, q_corr, start, count, threads);
    if ((TYPEOF(q_df) != INTSXP && TYPEOF(q_df) != REALSXP) || XLENGTH(q_df) != 1) {
        error("'q_df' must be a numeric scalar");
    }
    double q_degrees_freedom = asReal(q_df);
    if (!R_FINITE(q_degrees_freedom) || q_degrees_freedom < 0) {
        error("'q_df' must be a non-negative finite number");
    }
    if (inputs.factors > (INT_MAX - 3) / 4) {
        error("too many factors for an R result matrix");
    }

    int n = (int) inputs.row_count;
    int result_columns = 4 * inputs.factors + 3;
    SEXP results_out = PROTECT(allocMatrix(REALSXP, n, result_columns));
    SEXP status_out = PROTECT(allocVector(INTSXP, n));

    R_CheckUserInterrupt();
    int result = genomicsem_gls_batch_columns_results(
        inputs.betas,
        inputs.ses,
        inputs.row_count,
        (size_t) inputs.traits,
        REAL(loadings),
        (size_t) inputs.factors,
        REAL(corr),
        REAL(q_corr),
        (size_t) inputs.threads,
        q_degrees_freedom,
        two_sided_normal_p,
        upper_chi_square_p,
        REAL(results_out),
        INTEGER(status_out));
    R_CheckUserInterrupt();
    check_kernel_result(result);

    SEXP output = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(output, 0, results_out);
    SET_VECTOR_ELT(output, 1, status_out);

    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("results"));
    SET_STRING_ELT(names, 1, mkChar("status"));
    setAttrib(output, R_NamesSymbol, names);

    UNPROTECT(4);
    return output;
}
