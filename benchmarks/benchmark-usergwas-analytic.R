library(GenomicSEM)

args <- commandArgs(trailingOnly = TRUE)
n_snps <- if (length(args) >= 1L) as.integer(args[[1L]]) else 500000L
threads <- if (length(args) >= 2L) as.integer(args[[2L]]) else 4L
n_factors <- if (length(args) >= 3L) as.integer(args[[3L]]) else 2L
batch_size <- if (length(args) >= 4L) as.integer(args[[4L]]) else 100000L
repetitions <- if (length(args) >= 5L) as.integer(args[[5L]]) else 5L
n_traits <- 12L

set.seed(20260714)
traits <- paste0("T", seq_len(n_traits))
factors <- paste0("F", seq_len(n_factors))
betas <- matrix(rnorm(n_snps * n_traits, sd = 0.05), n_snps, n_traits)
ses <- matrix(runif(n_snps * n_traits, 0.02, 0.09), n_snps, n_traits)
loadings <- matrix(runif(n_traits * n_factors, 0.3, 1.2), n_traits, n_factors)
if (n_factors > 1L) {
  loadings <- qr.Q(qr(loadings))[, seq_len(n_factors), drop = FALSE]
}
q_corr <- matrix(0.06, n_traits, n_traits)
diag(q_corr) <- seq(0.9, 1.1, length.out = n_traits)
sampling_corr <- q_corr
diag(sampling_corr)[diag(sampling_corr) < 1] <- 1

sumstats <- data.frame(
  SNP = paste0("rs", seq_len(n_snps)),
  CHR = rep.int(1L, n_snps),
  BP = seq_len(n_snps),
  MAF = rep.int(0.2, n_snps),
  A1 = rep.int("A", n_snps),
  A2 = rep.int("G", n_snps),
  stringsAsFactors = FALSE
)
for (trait_index in seq_along(traits)) {
  sumstats[[paste0("beta.", traits[[trait_index]])]] <- betas[, trait_index]
}
for (trait_index in seq_along(traits)) {
  sumstats[[paste0("se.", traits[[trait_index]])]] <- ses[, trait_index]
}
rm(betas, ses)

ldsc <- list(S = diag(n_traits), I = q_corr)
colnames(ldsc$S) <- rownames(ldsc$S) <- traits
usermod <- do.call(rbind, lapply(seq_along(factors), function(factor_index) {
  data.frame(
    lhs = factors[[factor_index]], op = "=~", rhs = traits,
    Unstand_Est = loadings[, factor_index]
  )
}))
model <- paste(
  vapply(seq_along(factors), function(index) {
    paste0(factors[[index]], " =~ ", paste(traits, collapse = " + "))
  }, character(1)),
  vapply(factors, function(factor) paste0(factor, " ~ SNP"), character(1)),
  collapse = "\n"
)

now <- function() unname(proc.time()[["elapsed"]])

run_public <- function() {
  invisible(capture.output(
    result <- suppressWarnings(userGWAS(
      covstruc = ldsc,
      SNPs = sumstats,
      model = model,
      cores = threads,
      analytic = TRUE,
      batch_size = batch_size,
      usermod = usermod
    ))
  ))
  result
}

run_copying_input <- function() {
  phase <- c(initialize = 0, extract = 0, kernel = 0, output = 0, assemble = 0)
  start <- now()

  beta_columns <- match(paste0("beta.", traits), names(sumstats))
  se_columns <- match(paste0("se.", traits), names(sumstats))
  factor_beta_columns <- seq.int(1L, by = 4L, length.out = n_factors)
  factor_se_columns <- factor_beta_columns + 1L
  factor_z_columns <- factor_beta_columns + 2L
  factor_p_columns <- factor_beta_columns + 3L
  q_column <- 4L * n_factors + 1L
  q_df_column <- q_column + 1L
  q_p_column <- q_column + 2L
  numeric_names <- c(
    unlist(lapply(factors, function(factor) {
      c(
        paste0("beta_", factor), paste0("SE_", factor),
        paste0("Z_beta_", factor), paste0("p_val_", factor)
      )
    })),
    "Q_omnibus", "Q_omnibus_df", "Q_omnibus_pval"
  )
  numeric_result <- matrix(
    NA_real_, nrow = n_snps, ncol = length(numeric_names),
    dimnames = list(NULL, numeric_names)
  )
  numeric_result[, q_df_column] <- n_traits - n_factors
  phase[["initialize"]] <- now() - start

  for (batch_start in seq.int(1L, n_snps, by = batch_size)) {
    batch_end <- min(batch_start + batch_size - 1L, n_snps)
    rows <- batch_start:batch_end

    phase_start <- now()
    batch_betas <- as.matrix(sumstats[rows, beta_columns, drop = FALSE])
    batch_ses <- as.matrix(sumstats[rows, se_columns, drop = FALSE])
    phase[["extract"]] <- phase[["extract"]] + now() - phase_start

    phase_start <- now()
    native <- GenomicSEM:::.analytic_gls_batch_rust(
      batch_betas, batch_ses, loadings, sampling_corr, q_corr,
      threads = threads
    )
    phase[["kernel"]] <- phase[["kernel"]] + now() - phase_start

    phase_start <- now()
    z_values <- native$beta / native$se
    p_values <- matrix(
      2 * stats::pnorm(-abs(z_values)),
      nrow = length(rows), ncol = n_factors
    )
    numeric_result[rows, factor_beta_columns] <- native$beta
    numeric_result[rows, factor_se_columns] <- native$se
    numeric_result[rows, factor_z_columns] <- z_values
    numeric_result[rows, factor_p_columns] <- p_values
    numeric_result[rows, q_column] <- native$q
    numeric_result[rows, q_p_column] <- stats::pchisq(
      native$q, df = n_traits - n_factors, lower.tail = FALSE
    )
    phase[["output"]] <- phase[["output"]] + now() - phase_start
  }

  phase_start <- now()
  result <- data.frame(
    sumstats[, 1:6, drop = FALSE],
    as.data.frame(numeric_result, optional = TRUE),
    check.names = FALSE
  )
  phase[["assemble"]] <- now() - phase_start
  list(result = result, phase = phase, elapsed = now() - start)
}

# Warm native code and R method dispatch before collecting alternating runs.
public_result <- run_public()
copying_run <- run_copying_input()
stopifnot(
  identical(names(public_result), names(copying_run$result)),
  isTRUE(all.equal(public_result, copying_run$result, tolerance = 0))
)

public_times <- numeric(repetitions)
copying_times <- numeric(repetitions)
copying_phases <- matrix(
  0, nrow = repetitions, ncol = length(copying_run$phase),
  dimnames = list(NULL, names(copying_run$phase))
)

for (iteration in seq_len(repetitions)) {
  if (iteration %% 2L == 1L) {
    invisible(gc(FALSE))
    start <- now()
    public_result <- run_public()
    public_times[[iteration]] <- now() - start

    invisible(gc(FALSE))
    copying_run <- run_copying_input()
    copying_times[[iteration]] <- copying_run$elapsed
  } else {
    invisible(gc(FALSE))
    copying_run <- run_copying_input()
    copying_times[[iteration]] <- copying_run$elapsed

    invisible(gc(FALSE))
    start <- now()
    public_result <- run_public()
    public_times[[iteration]] <- now() - start
  }
  copying_phases[iteration, ] <- copying_run$phase
}

stopifnot(isTRUE(all.equal(public_result, copying_run$result, tolerance = 0)))

median_public <- median(public_times)
median_copying <- median(copying_times)
median_phases <- apply(copying_phases, 2L, median)

cat(sprintf(
  "Workload: %d SNPs, %d traits, %d factors, batch size %d, %d threads\n",
  n_snps, n_traits, n_factors, batch_size, threads
))
cat(sprintf("Input size: %.1f MiB\n", as.numeric(object.size(sumstats)) / 1024^2))
cat(sprintf("Repetitions: %d alternating runs after warmup\n", repetitions))
cat(sprintf("Current public median:  %.3f seconds\n", median_public))
cat(sprintf("Copying-input median:   %.3f seconds\n", median_copying))
cat(sprintf("  initialize:           %.3f seconds\n", median_phases[["initialize"]]))
cat(sprintf("  extract matrices:     %.3f seconds\n", median_phases[["extract"]]))
cat(sprintf("  Rust kernel:          %.3f seconds\n", median_phases[["kernel"]]))
cat(sprintf("  statistics + output:  %.3f seconds\n", median_phases[["output"]]))
cat(sprintf("  final assembly:       %.3f seconds\n", median_phases[["assemble"]]))
cat(sprintf("Zero-copy speedup:      %.2fx\n", median_copying / median_public))
