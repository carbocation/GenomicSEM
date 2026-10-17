#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <stddef.h>

extern int genomicsem_ldsc_regressions(
    const double *l2,
    const double *wld,
    const double *n,
    const double *z,
    const int *a1,
    const int *id,
    const int *complete,
    const int *offsets,
    size_t n_traits,
    size_t n_universe,
    const int *pairs_j,
    const int *pairs_k,
    size_t n_pairs,
    double m_tot,
    size_t n_blocks,
    size_t requested_threads,
    double *reg_tot_out,
    double *intercept_out,
    double *n_bar_out,
    double *n_snps_out,
    double *mean_y_out,
    double *intercept_se_out,
    double *tot_se_out,
    double *pseudo_out,
    int *status_out);

static void check_type(SEXP x, int type, const char *name) {
    if (TYPEOF(x) != type) {
        error("'%s' has the wrong storage type", name);
    }
}

/* Every heritability and genetic-covariance regression of ldsc(). Trait
   columns are concatenated across traits; offsets (length n_traits + 1)
   delimit them. Pairs are 0-based trait indices; j == k is a heritability. */
SEXP C_genomicsem_ldsc(
    SEXP l2, SEXP wld, SEXP n, SEXP z, SEXP a1, SEXP id, SEXP complete,
    SEXP offsets, SEXP n_universe, SEXP pairs_j, SEXP pairs_k, SEXP m_tot,
    SEXP n_blocks, SEXP threads) {
    check_type(l2, REALSXP, "l2");
    check_type(wld, REALSXP, "wld");
    check_type(n, REALSXP, "n");
    check_type(z, REALSXP, "z");
    check_type(a1, INTSXP, "a1");
    check_type(id, INTSXP, "id");
    check_type(complete, INTSXP, "complete");
    check_type(offsets, INTSXP, "offsets");
    check_type(pairs_j, INTSXP, "pairs_j");
    check_type(pairs_k, INTSXP, "pairs_k");

    R_xlen_t total = XLENGTH(l2);
    if (XLENGTH(wld) != total || XLENGTH(n) != total || XLENGTH(z) != total ||
        XLENGTH(a1) != total || XLENGTH(id) != total || XLENGTH(complete) != total) {
        error("trait columns must have equal lengths");
    }
    R_xlen_t n_traits = XLENGTH(offsets) - 1;
    if (n_traits < 1 || INTEGER(offsets)[0] != 0 || INTEGER(offsets)[n_traits] != total) {
        error("'offsets' must start at 0 and end at the column length");
    }
    for (R_xlen_t t = 0; t < n_traits; ++t) {
        if (INTEGER(offsets)[t + 1] < INTEGER(offsets)[t]) {
            error("'offsets' must be non-decreasing");
        }
    }
    int universe = asInteger(n_universe);
    if (universe == NA_INTEGER || universe < 1) {
        error("'n_universe' must be a positive integer");
    }
    for (R_xlen_t i = 0; i < total; ++i) {
        int v = INTEGER(id)[i];
        if (v == NA_INTEGER || v < 0 || v >= universe) {
            error("'id' values must lie in 0..n_universe - 1");
        }
    }
    R_xlen_t n_pairs = XLENGTH(pairs_j);
    if (XLENGTH(pairs_k) != n_pairs || n_pairs < 1) {
        error("'pairs_j' and 'pairs_k' must have the same positive length");
    }
    for (R_xlen_t p = 0; p < n_pairs; ++p) {
        int j = INTEGER(pairs_j)[p], k = INTEGER(pairs_k)[p];
        if (j == NA_INTEGER || k == NA_INTEGER || j < 0 || k < 0 || j >= n_traits || k >= n_traits) {
            error("pair indices must lie in 0..n_traits - 1");
        }
    }
    int blocks = asInteger(n_blocks);
    if (blocks == NA_INTEGER || blocks < 2) {
        error("'n_blocks' must be at least 2");
    }
    int thread_count = asInteger(threads);
    if (thread_count == NA_INTEGER || thread_count < 1) {
        error("'threads' must be at least 1");
    }

    const char *names[] = {"reg_tot", "intercept", "n_bar", "n_snps", "mean_y",
                           "intercept_se", "tot_se", "pseudo", "status", ""};
    SEXP out = PROTECT(mkNamed(VECSXP, names));
    double *vec_out[7];
    for (int i = 0; i < 7; ++i) {
        SEXP v = allocVector(REALSXP, n_pairs);
        SET_VECTOR_ELT(out, i, v);
        vec_out[i] = REAL(v);
    }
    SEXP pseudo = allocMatrix(REALSXP, blocks, (int) n_pairs);
    SET_VECTOR_ELT(out, 7, pseudo);
    SEXP status = allocVector(INTSXP, n_pairs);
    SET_VECTOR_ELT(out, 8, status);

    int result = genomicsem_ldsc_regressions(
        REAL(l2), REAL(wld), REAL(n), REAL(z), INTEGER(a1), INTEGER(id),
        INTEGER(complete), INTEGER(offsets), (size_t) n_traits, (size_t) universe,
        INTEGER(pairs_j), INTEGER(pairs_k), (size_t) n_pairs, asReal(m_tot),
        (size_t) blocks, (size_t) thread_count,
        vec_out[0], vec_out[1], vec_out[2], vec_out[3], vec_out[4], vec_out[5],
        vec_out[6], REAL(pseudo), INTEGER(status));
    if (result != 0) {
        UNPROTECT(1);
        error("the Rust LD-score regression kernel failed unexpectedly");
    }
    UNPROTECT(1);
    return out;
}
