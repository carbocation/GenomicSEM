.analytic_gls_batch_r <- function(betas, ses, loadings, corr, q_corr = corr) {
  betas <- as.matrix(betas)
  ses <- as.matrix(ses)
  loadings <- as.matrix(loadings)
  corr <- as.matrix(corr)
  q_corr <- as.matrix(q_corr)

  n_snps <- nrow(betas)
  n_factors <- ncol(loadings)
  beta_out <- matrix(NA_real_, nrow = n_snps, ncol = n_factors)
  se_out <- matrix(NA_real_, nrow = n_snps, ncol = n_factors)
  q_out <- rep(NA_real_, n_snps)
  corr_inv <- solve(q_corr)

  for (row in seq_len(n_snps)) {
    trait_se <- ses[row, ]
    trait_beta <- betas[row, ]
    sampling_cov <- tcrossprod(trait_se) * corr
    weights <- 1 / diag(sampling_cov)
    normal <- crossprod(loadings, loadings * weights)
    bread <- solve(normal)
    factor_beta <- bread %*% crossprod(loadings, weights * trait_beta)
    weighted_loadings <- loadings * weights
    meat <- crossprod(weighted_loadings, sampling_cov %*% weighted_loadings)
    sandwich <- bread %*% meat %*% bread
    residual <- (trait_beta - drop(loadings %*% factor_beta)) / trait_se

    beta_out[row, ] <- factor_beta
    se_out[row, ] <- sqrt(diag(sandwich))
    q_out[row] <- drop(crossprod(residual, corr_inv %*% residual))
  }

  list(beta = beta_out, se = se_out, q = q_out, status = integer(n_snps))
}

.check_analytic_gls_output <- function(output) {
  bad <- which(output$status != 0L)
  if (length(bad) > 0L) {
    status_messages <- c(
      "non-finite beta, SE, or loading, or a non-positive SE/intercept diagonal",
      "singular factor normal matrix",
      "non-finite or negative sampling variance"
    )
    first <- bad[[1L]]
    code <- output$status[[first]]
    stop(
      "The analytic estimator failed for batch row ", first, ": ",
      status_messages[[code]],
      call. = FALSE
    )
  }
  output
}

.analytic_gls_batch_rust <- function(betas, ses, loadings, corr, q_corr = corr,
                                     threads = 1L) {
  betas <- as.matrix(betas)
  ses <- as.matrix(ses)
  loadings <- as.matrix(loadings)
  corr <- as.matrix(corr)
  q_corr <- as.matrix(q_corr)
  storage.mode(betas) <- "double"
  storage.mode(ses) <- "double"
  storage.mode(loadings) <- "double"
  storage.mode(corr) <- "double"
  storage.mode(q_corr) <- "double"

  output <- .Call(
    .genomicsem_gls_batch,
    betas,
    ses,
    loadings,
    corr,
    q_corr,
    as.integer(threads)
  )

  .check_analytic_gls_output(output)
}

.analytic_gls_columns_rust <- function(betas, ses, loadings, corr, q_corr = corr,
                                       start = 0L, count, threads = 1L) {
  if (!is.list(betas) || !is.list(ses)) {
    stop("'betas' and 'ses' must be lists of numeric columns.", call. = FALSE)
  }
  if (!all(vapply(betas, is.double, logical(1))) ||
      !all(vapply(ses, is.double, logical(1)))) {
    stop("All beta and SE columns must use double storage.", call. = FALSE)
  }
  loadings <- as.matrix(loadings)
  corr <- as.matrix(corr)
  q_corr <- as.matrix(q_corr)
  storage.mode(loadings) <- "double"
  storage.mode(corr) <- "double"
  storage.mode(q_corr) <- "double"

  output <- .Call(
    .genomicsem_gls_columns,
    betas,
    ses,
    loadings,
    corr,
    q_corr,
    start,
    count,
    as.integer(threads)
  )

  .check_analytic_gls_output(output)
}

.analytic_gls_results_columns_rust <- function(
    betas, ses, loadings, corr, q_corr = corr,
    start = 0L, count, threads = 1L, q_df) {
  if (!is.list(betas) || !is.list(ses)) {
    stop("'betas' and 'ses' must be lists of numeric columns.", call. = FALSE)
  }
  if (!all(vapply(betas, is.double, logical(1))) ||
      !all(vapply(ses, is.double, logical(1)))) {
    stop("All beta and SE columns must use double storage.", call. = FALSE)
  }
  loadings <- as.matrix(loadings)
  corr <- as.matrix(corr)
  q_corr <- as.matrix(q_corr)
  storage.mode(loadings) <- "double"
  storage.mode(corr) <- "double"
  storage.mode(q_corr) <- "double"

  output <- .Call(
    .genomicsem_gls_results_columns,
    betas,
    ses,
    loadings,
    corr,
    q_corr,
    start,
    count,
    as.integer(threads),
    q_df
  )

  .check_analytic_gls_output(output)$results
}

.analytic_gls_batch <- function(betas, ses, loadings, corr, q_corr = corr, threads = 1L,
                                backend = getOption("GenomicSEM.analytic_backend", "rust")) {
  backend <- match.arg(backend, c("rust", "R"))
  if (backend == "rust") {
    return(.analytic_gls_batch_rust(betas, ses, loadings, corr, q_corr, threads))
  }
  .analytic_gls_batch_r(betas, ses, loadings, corr, q_corr)
}
