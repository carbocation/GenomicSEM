make_gls_inputs <- function(n_snps = 200L, n_traits = 8L, n_factors = 2L,
                            seed = 123L) {
  set.seed(seed)
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

  list(
    betas = betas,
    ses = ses,
    loadings = loadings,
    sampling_corr = sampling_corr,
    q_corr = q_corr
  )
}

test_that("Rust analytic estimates agree with the R reference", {
  for (n_factors in c(1L, 2L, 3L, 4L)) {
    inputs <- make_gls_inputs(n_factors = n_factors, seed = 100L + n_factors)
    reference <- GenomicSEM:::.analytic_gls_batch_r(
      inputs$betas,
      inputs$ses,
      inputs$loadings,
      inputs$sampling_corr,
      inputs$q_corr
    )
    rust <- GenomicSEM:::.analytic_gls_batch_rust(
      inputs$betas,
      inputs$ses,
      inputs$loadings,
      inputs$sampling_corr,
      inputs$q_corr,
      threads = 1L
    )

    expect_equal(rust$beta, reference$beta, tolerance = 1e-11)
    expect_equal(rust$se, reference$se, tolerance = 1e-11)
    expect_equal(rust$q, reference$q, tolerance = 1e-10)
    expect_identical(rust$status, integer(nrow(inputs$betas)))
  }
})

test_that("Rust thread counts produce deterministic results", {
  for (n_factors in 1:3) {
    inputs <- make_gls_inputs(n_snps = 503L, n_factors = n_factors)
    serial <- GenomicSEM:::.analytic_gls_batch_rust(
      inputs$betas, inputs$ses, inputs$loadings,
      inputs$sampling_corr, inputs$q_corr, threads = 1L
    )
    parallel <- GenomicSEM:::.analytic_gls_batch_rust(
      inputs$betas, inputs$ses, inputs$loadings,
      inputs$sampling_corr, inputs$q_corr, threads = 4L
    )
    expect_identical(parallel, serial)
  }
})

test_that("two-factor specialization remains accurate for correlated loadings", {
  inputs <- make_gls_inputs(n_snps = 257L, n_traits = 10L, n_factors = 2L)
  inputs$loadings[, 2] <- inputs$loadings[, 1] +
    seq(-0.025, 0.025, length.out = nrow(inputs$loadings))
  reference <- GenomicSEM:::.analytic_gls_batch_r(
    inputs$betas, inputs$ses, inputs$loadings,
    inputs$sampling_corr, inputs$q_corr
  )
  rust <- GenomicSEM:::.analytic_gls_batch_rust(
    inputs$betas, inputs$ses, inputs$loadings,
    inputs$sampling_corr, inputs$q_corr, threads = 2L
  )

  expect_equal(rust$beta, reference$beta, tolerance = 1e-9)
  expect_equal(rust$se, reference$se, tolerance = 1e-9)
  expect_equal(rust$q, reference$q, tolerance = 1e-9)
})

test_that("precomputed Q inverse preserves LU pivoting semantics", {
  inputs <- make_gls_inputs(n_snps = 31L, n_traits = 4L, n_factors = 1L)
  inputs$q_corr <- diag(4L)
  inputs$q_corr[1:2, 1:2] <- matrix(c(0, 0.4, 0.4, 0), 2L)
  inputs$sampling_corr <- inputs$q_corr
  diag(inputs$sampling_corr) <- 1
  reference <- GenomicSEM:::.analytic_gls_batch_r(
    inputs$betas, inputs$ses, inputs$loadings,
    inputs$sampling_corr, inputs$q_corr
  )
  rust <- GenomicSEM:::.analytic_gls_batch_rust(
    inputs$betas, inputs$ses, inputs$loadings,
    inputs$sampling_corr, inputs$q_corr
  )

  expect_equal(rust$beta, reference$beta, tolerance = 1e-11)
  expect_equal(rust$se, reference$se, tolerance = 1e-11)
  expect_equal(rust$q, reference$q, tolerance = 1e-10)
})

test_that("singular factor models report the failing batch row", {
  inputs <- make_gls_inputs(n_snps = 5L, n_factors = 2L)
  inputs$loadings[, 2] <- inputs$loadings[, 1]
  expect_error(
    GenomicSEM:::.analytic_gls_batch_rust(
      inputs$betas, inputs$ses, inputs$loadings,
      inputs$sampling_corr, inputs$q_corr
    ),
    "singular factor normal matrix"
  )
})

test_that("invalid sampling inputs fail with an informative row", {
  inputs <- make_gls_inputs(n_snps = 5L, n_factors = 2L)
  inputs$ses[3, 2] <- 0
  expect_error(
    GenomicSEM:::.analytic_gls_batch_rust(
      inputs$betas, inputs$ses, inputs$loadings,
      inputs$sampling_corr, inputs$q_corr
    ),
    "batch row 3.*non-positive SE"
  )
})

test_that("analytic batching requires whole-number workload controls", {
  expect_error(
    suppressMessages(GenomicSEM:::.userGWASa(
      data.frame(x = 1), NULL, model = "", batch_size = 1.5
    )),
    "positive whole number"
  )
  expect_error(
    suppressMessages(GenomicSEM:::.userGWASa(
      data.frame(x = 1), NULL, model = "", batch_size = 1L, cores = 1.5
    )),
    "positive whole number"
  )
})

test_that("userGWAS analytic batching preserves output", {
  inputs <- make_gls_inputs(n_snps = 211L, n_traits = 8L, n_factors = 2L)
  traits <- paste0("T", seq_len(8L))
  factors <- c("F1", "F2")
  sumstats <- data.frame(
    SNP = paste0("rs", seq_len(211L)), CHR = 1L, BP = seq_len(211L),
    MAF = 0.2, A1 = "A", A2 = "G"
  )
  for (trait_index in seq_along(traits)) {
    sumstats[[paste0("beta.", traits[[trait_index]])]] <-
      inputs$betas[, trait_index]
  }
  for (trait_index in seq_along(traits)) {
    sumstats[[paste0("se.", traits[[trait_index]])]] <-
      inputs$ses[, trait_index]
  }
  ldsc <- list(S = diag(8L), I = inputs$q_corr)
  colnames(ldsc$S) <- rownames(ldsc$S) <- traits
  usermod <- do.call(rbind, lapply(seq_along(factors), function(factor_index) {
    data.frame(
      lhs = factors[[factor_index]], op = "=~", rhs = traits,
      Unstand_Est = inputs$loadings[, factor_index]
    )
  }))

  reference <- suppressWarnings(GenomicSEM:::.userGWASa(
    sumstats, ldsc, model = "", usermod = usermod,
    batch_size = 37L, cores = 1L, backend = "R"
  ))
  rust <- suppressWarnings(GenomicSEM:::.userGWASa(
    sumstats, ldsc, model = "", usermod = usermod,
    batch_size = 29L, cores = 4L, backend = "rust"
  ))

  expect_identical(names(rust), names(reference))
  expect_identical(rust[, 1:6], reference[, 1:6])
  expect_equal(rust[, 7:ncol(rust)], reference[, 7:ncol(reference)],
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("userGWAS aligns trait columns by LDSC names", {
  inputs <- make_gls_inputs(n_snps = 13L, n_traits = 4L, n_factors = 1L)
  traits <- paste0("T", seq_len(4L))
  sumstats <- data.frame(
    SNP = paste0("rs", seq_len(13L)), CHR = 1L, BP = seq_len(13L),
    MAF = 0.2, A1 = "A", A2 = "G"
  )
  for (trait_index in rev(seq_along(traits))) {
    sumstats[[paste0("se.", traits[[trait_index]])]] <- inputs$ses[, trait_index]
    sumstats[[paste0("beta.", traits[[trait_index]])]] <- inputs$betas[, trait_index]
  }
  ldsc <- list(S = diag(4L), I = inputs$q_corr)
  colnames(ldsc$S) <- rownames(ldsc$S) <- traits
  usermod <- data.frame(
    lhs = "F1", op = "=~", rhs = traits,
    Unstand_Est = inputs$loadings[, 1]
  )

  model <- paste0("F1 =~ ", paste(traits, collapse = " + "), "\nF1 ~ SNP")
  result <- suppressWarnings(userGWAS(
    covstruc = ldsc, SNPs = sumstats, model = model,
    analytic = TRUE, usermod = usermod, batch_size = 7L, cores = 2L
  ))
  reference <- GenomicSEM:::.analytic_gls_batch_rust(
    inputs$betas, inputs$ses, inputs$loadings,
    pmax(inputs$q_corr, diag(4L)), inputs$q_corr
  )

  expect_equal(result$beta_F1, reference$beta[, 1], tolerance = 1e-11)
  expect_equal(result$SE_F1, reference$se[, 1], tolerance = 1e-11)
  expect_equal(result$Q_omnibus, reference$q, tolerance = 1e-10)
})
