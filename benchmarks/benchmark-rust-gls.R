library(GenomicSEM)

args <- commandArgs(trailingOnly = TRUE)
n_snps <- if (length(args)) as.integer(args[[1L]]) else 50000L
n_traits <- 12L
n_factors <- if (length(args) >= 3L) as.integer(args[[3L]]) else 2L
if (is.na(n_factors) || n_factors < 1L || n_factors > n_traits) {
  stop("The factor count must be between 1 and ", n_traits, ".")
}
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) detected_cores <- 1L
threads <- max(1L, min(4L, detected_cores))
if (length(args) >= 2L) threads <- max(1L, as.integer(args[[2L]]))

set.seed(20260714)
betas <- matrix(rnorm(n_snps * n_traits, sd = 0.05), n_snps, n_traits)
ses <- matrix(runif(n_snps * n_traits, 0.02, 0.09), n_snps, n_traits)
loadings <- qr.Q(qr(matrix(
  runif(n_traits * n_factors, 0.3, 1.2), n_traits, n_factors
)))[, seq_len(n_factors), drop = FALSE]
q_corr <- matrix(0.06, n_traits, n_traits)
diag(q_corr) <- seq(0.9, 1.1, length.out = n_traits)
sampling_corr <- q_corr
diag(sampling_corr)[diag(sampling_corr) < 1] <- 1

timed <- function(label, expression) {
  gc()
  timing <- system.time(value <- force(expression))
  cat(sprintf("%-18s %.3f seconds\n", label, timing[["elapsed"]]))
  list(value = value, elapsed = unname(timing[["elapsed"]]))
}

cat(sprintf("Workload: %d SNPs, %d traits, %d factors\n",
            n_snps, n_traits, n_factors))
reference <- timed("R reference", GenomicSEM:::.analytic_gls_batch_r(
  betas, ses, loadings, sampling_corr, q_corr
))
rust_serial <- timed("Rust, 1 thread", GenomicSEM:::.analytic_gls_batch_rust(
  betas, ses, loadings, sampling_corr, q_corr, threads = 1L
))
rust_parallel <- timed(sprintf("Rust, %d threads", threads),
  GenomicSEM:::.analytic_gls_batch_rust(
    betas, ses, loadings, sampling_corr, q_corr, threads = threads
  ))

stopifnot(
  isTRUE(all.equal(reference$value$beta, rust_serial$value$beta, tolerance = 1e-11)),
  isTRUE(all.equal(reference$value$se, rust_serial$value$se, tolerance = 1e-11)),
  isTRUE(all.equal(reference$value$q, rust_serial$value$q, tolerance = 1e-10)),
  identical(rust_serial$value, rust_parallel$value)
)
cat(sprintf("Rust, 1 thread speedup: %.1fx\n",
            reference$elapsed / rust_serial$elapsed))
cat(sprintf("Rust, %d threads speedup: %.1fx\n",
            threads, reference$elapsed / rust_parallel$elapsed))
