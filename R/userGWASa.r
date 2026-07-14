#' userGWASa: Ultra-fast multivariate GWAS with flexible analytic estimation
#'
#' Runs a multivariate GWAS across a set of
#' GWAS summary statistics and a user-specified factor model. Factor-specific
#' SNP effects (betas, SEs, Z-statistics, p-values) and an omnibus
#' heterogeneity statistic (Q_omnibus) are computed.
#'
#' @param sumstats A \code{data.frame} of merged GWAS summary statistics,
#'   as produced by the \code{sumstats()} function in GenomicSEM. Must contain
#'   columns \code{SNP}, \code{A1}, \code{A2}, \code{MAF}, \code{N}, and
#'   trait-specific \code{beta.*} and \code{se.*} columns.
#' @param LDSCoutput A list object returned by the \code{ldsc()} function.
#' @param model A character string specifying the factor model in
#'   \code{lavaan}-style syntax. Ignored if
#'   \code{usermod} is provided.
#' @param usermod Optional. A pre-fitted no-SNP model results data frame
#'   (the \code{$results} element from a \code{usermodel()} call). When
#'   supplied, the function skips fitting the no-SNP model internally and uses
#'   these parameter estimates directly to extract lambda coefficients.
#'   Default is \code{NULL}.
#' @param batch_size Integer. Number of SNPs to process per batch. Larger
#'   values increase memory use but reduce overhead. Default is \code{100000}.
#' @param cores Optional positive integer controlling the number of Rust worker
#'   threads. The default, \code{NULL}, uses one thread.
#' @param backend Internal diagnostic backend selector. Either \code{"rust"}
#'   (the default) or \code{"R"} for the reference implementation.
#'
#' @return A \code{data.frame} with one row per SNP and the following columns:
#'   \itemize{
#'     \item The first 6 columns from \code{sumstats} (SNP identifiers and
#'           allele information).
#'     \item \code{beta_<factor>}: GLS-estimated SNP effect on each factor.
#'     \item \code{SE_<factor>}: Sandwich-corrected standard error of the
#'           factor beta.
#'     \item \code{Z_beta_<factor>}: Z-statistic for the factor beta.
#'     \item \code{p_val_<factor>}: Two-sided p-value for the factor beta.
#'     \item \code{Q_omnibus}, \code{Q_omnibus_df}, \code{Q_omnibus_pval}:
#'           Omnibus Q_SNP statistic across all traits, its degrees of freedom,
#'           and p-value.
#'   }
#'
#' @details
#' The function implements a two-stage approach. First, a no-SNP factor model
#' is fit using \code{\link[GenomicSEM]{usermodel}} with DWLS estimation to
#' obtain factor loading estimates (lambdas). Second, for each batch of SNPs,
#' SNP-to-factor betas are estimated via GLS using the diagonal of the
#' SNP-specific sampling covariance matrix as weights, with a sandwich
#' variance estimator for the standard errors.
#'
#' The Q_omnibus statistic tests whether the observed SNP-trait association
#' vector is consistent with the implied factor model.
#'
#' @seealso \code{\link[GenomicSEM]{usermodel}}, \code{\link[GenomicSEM]{userGWAS}}
#' @noRd
#'
#' @examples
#' \dontrun{
#' load("LDSC_PSYCH.RData")
#' sumstats <- data.table::fread("Psych_sumstats_4GLS.txt", data.table = FALSE)
#'
#' model <- '
#'   Psych =~ a*SCZ + a*BIP
#'   Neuro =~ ADHD + MDD + ASD
#'   Psych ~~ Neuro
#'   Psych ~ SNP
#'   Neuro ~ SNP
#' '
#'
#' results <- userGWASa(
#'   sumstats   = sumstats,
#'   LDSCoutput = LDSC_P,
#'   model      = model,
#'   batch_size = 50000
#' )
#' }
#'
#' @keywords internal
.userGWASa <- function(sumstats, LDSCoutput, model, usermod = NULL,
                       batch_size = 100000, cores = NULL,
                       backend = getOption("GenomicSEM.analytic_backend", "rust")) {

  start_time <- Sys.time()
  cat("userGWASa started at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

  # Coerce to plain data.frame to ensure consistent column subsetting
  # regardless of whether input is a data.table, tibble, or other tabular class
  sumstats <- as.data.frame(sumstats)
  if (nrow(sumstats) == 0L) {
    stop("'sumstats' must contain at least one SNP.", call. = FALSE)
  }
  if (length(batch_size) != 1L || !is.numeric(batch_size) ||
      !is.finite(batch_size) || batch_size < 1 || batch_size != floor(batch_size)) {
    stop("'batch_size' must be a positive whole number.", call. = FALSE)
  }
  if (batch_size > .Machine$integer.max) {
    stop("'batch_size' is too large for an R matrix.", call. = FALSE)
  }
  batch_size <- as.integer(batch_size)
  if (is.null(cores)) {
    threads <- 1L
  } else {
    if (length(cores) != 1L || !is.numeric(cores) || !is.finite(cores) ||
        cores < 1 || cores != floor(cores)) {
      stop("'cores' must be a positive whole number.", call. = FALSE)
    }
    if (cores > .Machine$integer.max) {
      stop("'cores' is too large.", call. = FALSE)
    }
    threads <- as.integer(cores)
  }
  backend <- match.arg(backend, c("rust", "R"))

  # ── No-SNP model ──────────────────────────────────────────────────────────────
  if (is.character(model) & is.null(usermod)) {
    model_lines  <- strsplit(model, "\n")[[1]]
    nosnp_model  <- paste(grep("~.*\\bSNP\\b", model_lines, value = TRUE, invert = TRUE),
                          collapse = "\n")

    captured_output <- capture.output(
    suppressWarnings(suppressMessages({
        nosnpmod <- usermodel(
        LDSCoutput, estimation = "DWLS", model = nosnp_model,
        CFIcalc = FALSE, std.lv = FALSE, imp_cov = FALSE
        )
    }))
    )
    nosnpmod <- nosnpmod$results

    # Re-emit smoothing warning if it occurred
    if (any(grepl("smoothed", captured_output))) {
    warning("The S matrix was smoothed prior to model estimation. ")
    }

  } else {
    nosnpmod <- usermod
  }

  # ── Extract lambda coefficients ───────────────────────────────────────────────
  factors    <- unique(nosnpmod$lhs[nosnpmod$op == "=~"])
  traits     <- colnames(LDSCoutput$S)
  num_traits  <- ncol(LDSCoutput$S)
  num_factors <- length(factors)
  if (num_factors == 0L) {
    stop("No factor loadings were found in the no-SNP model.", call. = FALSE)
  }

  combinations <- expand.grid(traits = traits, factors = factors)
  column_names <- paste0("lambda.", combinations$traits, "_", combinations$factors)

  extract_lambdas <- function(df, factors, traits, num_traits, num_factors) {
    lambdas <- rep(0, num_traits * num_factors)
    for (factor_idx in seq_along(factors)) {
      factor <- factors[factor_idx]
      for (trait_idx in seq_along(traits)) {
        trait <- traits[trait_idx]
        row   <- df[df$lhs == factor & df$rhs == trait & df$op == "=~", ]
        lambda_value <- if (nrow(row) > 0) row$Unstand_Est else 0
        lambdas[(factor_idx - 1) * num_traits + trait_idx] <- lambda_value
      }
    }
    lambdas_df <- data.frame(matrix(lambdas, nrow = 1, byrow = TRUE))
    colnames(lambdas_df) <- column_names
    return(lambdas_df)
  }

  lambdas <- extract_lambdas(nosnpmod, factors, traits, num_traits, num_factors)
  loadings <- matrix(as.numeric(lambdas), nrow = num_traits, ncol = num_factors)
  colnames(loadings) <- factors

  sampling_corr <- as.matrix(LDSCoutput$I)
  q_corr <- sampling_corr
  diag(sampling_corr)[diag(sampling_corr) < 1] <- 1

  beta_names <- paste0("beta.", traits)
  se_names <- paste0("se.", traits)
  missing_columns <- setdiff(c(beta_names, se_names), colnames(sumstats))
  if (length(missing_columns) > 0L) {
    stop(
      "'sumstats' is missing trait columns required by 'LDSCoutput': ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  beta_columns <- match(beta_names, colnames(sumstats))
  se_columns <- match(se_names, colnames(sumstats))
  if (backend == "rust") {
    beta_data <- lapply(sumstats[beta_columns], function(column) {
      if (is.double(column)) column else as.double(column)
    })
    se_data <- lapply(sumstats[se_columns], function(column) {
      if (is.double(column)) column else as.double(column)
    })
  }

  # ── Preallocate numeric output ────────────────────────────────────────────────
  # Filling a matrix avoids repeated [<-.data.frame dispatch and column copies
  # for every batch. Metadata is joined once after all numeric results are ready.
  factor_beta_columns <- seq.int(1L, by = 4L, length.out = num_factors)
  factor_se_columns <- factor_beta_columns + 1L
  factor_z_columns <- factor_beta_columns + 2L
  factor_p_columns <- factor_beta_columns + 3L
  q_column <- 4L * num_factors + 1L
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
  numeric_results <- matrix(
    NA_real_,
    nrow = nrow(sumstats),
    ncol = length(numeric_names),
    dimnames = list(NULL, numeric_names)
  )
  numeric_results[, q_df_column] <- num_traits - num_factors

  # ── Batch loop ────────────────────────────────────────────────────────────────
  total_batches <- ceiling(nrow(sumstats) / batch_size)
  pb <- txtProgressBar(min = 0, max = total_batches, style = 3)

  for (batch_num in seq_len(total_batches)) {

    i             <- (batch_num - 1) * batch_size + 1
    batch_end     <- min(i + batch_size - 1, nrow(sumstats))
    batch_indices <- i:batch_end
    if (backend == "rust") {
      batch_results <- .analytic_gls_results_columns_rust(
        betas = beta_data,
        ses = se_data,
        loadings = loadings,
        corr = sampling_corr,
        q_corr = q_corr,
        start = i - 1L,
        count = length(batch_indices),
        threads = threads,
        q_df = num_traits - num_factors
      )
      numeric_results[batch_indices, ] <- batch_results
    } else {
      betas <- as.matrix(sumstats[batch_indices, beta_columns, drop = FALSE])
      ses <- as.matrix(sumstats[batch_indices, se_columns, drop = FALSE])
      kernel_results <- .analytic_gls_batch(
        betas = betas,
        ses = ses,
        loadings = loadings,
        corr = sampling_corr,
        q_corr = q_corr,
        threads = threads,
        backend = backend
      )
      z_values <- kernel_results$beta / kernel_results$se

      # ── Write batch results ────────────────────────────────────────────────────
      p_values <- matrix(
        2 * stats::pnorm(-abs(z_values)),
        nrow = length(batch_indices),
        ncol = num_factors
      )
      numeric_results[batch_indices, factor_beta_columns] <- kernel_results$beta
      numeric_results[batch_indices, factor_se_columns] <- kernel_results$se
      numeric_results[batch_indices, factor_z_columns] <- z_values
      numeric_results[batch_indices, factor_p_columns] <- p_values
      numeric_results[batch_indices, q_column] <- kernel_results$q
      numeric_results[batch_indices, q_p_column] <- stats::pchisq(
        kernel_results$q,
        df = num_traits - num_factors,
        lower.tail = FALSE
      )
    }

    setTxtProgressBar(pb, batch_num)
  }

  close(pb)

  GLS_mGWAS_results <- data.frame(
    sumstats[, 1:6, drop = FALSE],
    as.data.frame(numeric_results, optional = TRUE),
    check.names = FALSE
  )

  end_time     <- Sys.time()
  elapsed_time <- end_time - start_time
  cat("Finished at:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total time elapsed:", round(elapsed_time, 2), attr(elapsed_time, "units"), "\n")

  return(GLS_mGWAS_results)
}
