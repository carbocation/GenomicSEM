make_ldsc_fixture <- function(dir, n_per_chr = 1500L, n_traits = 4L, seed = 7L) {
  set.seed(seed)
  ld_dir <- file.path(dir, "ld")
  dir.create(ld_dir, recursive = TRUE, showWarnings = FALSE)
  snps <- list()
  for (chr in 1:22) {
    n <- n_per_chr + chr * 10L
    bp <- sort(sample.int(5e7, n))
    l2 <- rgamma(n, shape = 2, scale = 30)
    d <- data.frame(CHR = chr, SNP = sprintf("rs%d_%d", chr, seq_len(n)), BP = bp,
                    CM = 0, MAF = runif(n, 0.05, 0.5), L2 = round(l2, 3))
    con <- gzfile(file.path(ld_dir, paste0(chr, ".l2.ldscore.gz")), "w")
    write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE)
    close(con)
    writeLines(as.character(n * 3L), file.path(ld_dir, paste0(chr, ".l2.M_5_50")))
    snps[[chr]] <- d
  }
  ld <- do.call(rbind, snps)
  m <- nrow(ld) * 3
  # Correlated per-SNP effects: Z ~ N(0, 1 + N h2 L2 / M) with shared signal.
  shared <- rnorm(nrow(ld))
  files <- character(n_traits)
  for (t in seq_len(n_traits)) {
    n_t <- c(50000, 80000, 120000, 30000)[t]
    h2 <- c(0.2, 0.3, 0.1, 0.4)[t]
    signal <- sqrt(n_t * h2 * ld$L2 / m) * (0.6 * shared + 0.8 * rnorm(nrow(ld)))
    z <- signal + rnorm(nrow(ld), sd = sqrt(1.02))
    a1 <- sample(c("A", "C", "G", "T"), nrow(ld), replace = TRUE)
    keep <- runif(nrow(ld)) > 0.08                 # SNPs missing per trait
    keep[nrow(ld)] <- FALSE                        # an LD-score SNP no trait has
    d <- data.frame(SNP = ld$SNP, N = round(n_t * runif(nrow(ld), 0.9, 1)), Z = round(z, 5),
                    A1 = a1, A2 = "N")[keep, ]
    if (t == 2L) d$Z[5] <- 40                      # dropped by the chi-square filter
    if (t > 1L) {                                  # the other traits are re-coded to A2
      flip <- runif(nrow(d)) < 0.3
      d$A1[flip] <- paste0(d$A1[flip], "x")
      d$Z[flip] <- -d$Z[flip]
    }
    d <- d[sample.int(nrow(d)), ]                  # file order differs from LD order
    files[t] <- file.path(dir, paste0("trait", t, ".sumstats.gz"))
    con <- gzfile(files[t], "w")
    write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE)
    close(con)
  }
  list(ld = ld_dir, traits = files)
}

run_ldsc_backend <- function(fixture, backend, dir, ..., threads = 2L) {
  old <- options(GenomicSEM.ldsc_backend = backend, GenomicSEM.ldsc_threads = threads)
  on.exit(options(old))
  log <- file.path(dir, paste0("ldsc-", backend, "-", threads))
  out <- suppressMessages(capture.output(res <- ldsc(
    traits = fixture$traits, sample.prev = c(NA, 0.3, NA, NA),
    population.prev = c(NA, 0.05, NA, NA), ld = fixture$ld, wld = fixture$ld,
    trait.names = paste0("T", seq_along(fixture$traits)), ldsc.log = log, stand = TRUE, ...)))
  lines <- readLines(paste0(log, "_ldsc.log"))
  list(res = res, log = lines[!grepl("began at|finished running|took", lines)])
}

test_that("Rust ldsc() agrees with the reference R implementation", {
  dir <- tempfile("ldsc-fixture")
  dir.create(dir)
  fixture <- make_ldsc_fixture(dir)
  reference <- run_ldsc_backend(fixture, "R", dir)
  rust <- run_ldsc_backend(fixture, "rust", dir)

  for (name in c("S", "V", "I", "N", "m", "S_Stand", "V_Stand")) {
    expect_equal(rust$res[[name]], reference$res[[name]], tolerance = 1e-10, info = name)
  }
  expect_identical(rust$log, reference$log)

  serial <- run_ldsc_backend(fixture, "rust", dir, threads = 1L)
  expect_identical(serial$res, rust$res)
})

test_that("Rust ldsc() agrees with the reference with chisq.max and selected chromosomes", {
  dir <- tempfile("ldsc-fixture")
  dir.create(dir)
  fixture <- make_ldsc_fixture(dir, n_traits = 3L, seed = 11L)
  fixture$traits <- fixture$traits[1:3]
  reference <- run_ldsc_backend(fixture, "R", dir, select = "ODD", chisq.max = 30)
  rust <- run_ldsc_backend(fixture, "rust", dir, select = "ODD", chisq.max = 30)
  for (name in c("S", "V", "I", "N", "m", "S_Stand", "V_Stand")) {
    expect_equal(rust$res[[name]], reference$res[[name]], tolerance = 1e-10, info = name)
  }
  expect_identical(rust$log, reference$log)
})
