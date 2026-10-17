ldsc <- function(traits, sample.prev, population.prev, ld, wld,
                trait.names = NULL, sep_weights = FALSE, chr = 22,
                n.blocks = 200, ldsc.log = NULL, stand = FALSE,select=FALSE,chisq.max = NA) {

  # The regressions run in the Rust kernel (src/rust/src/ldsc.rs) unless
  # options(GenomicSEM.ldsc_backend = "R") selects the reference R
  # implementation (.ldsc_legacy, R/ldsc_legacy.R). Inputs the fast path does
  # not reproduce exactly (duplicated SNP IDs, extra LD-score columns, a block
  # with no SNPs, a singular regression) also use the reference
  # implementation. Worker threads: options(GenomicSEM.ldsc_threads).
  backend <- match.arg(getOption("GenomicSEM.ldsc_backend", "rust"), c("rust", "R"))
  args <- list(traits = traits, sample.prev = sample.prev, population.prev = population.prev,
               ld = ld, wld = wld, trait.names = trait.names, sep_weights = sep_weights, chr = chr,
               n.blocks = n.blocks, ldsc.log = ldsc.log, stand = stand, select = select,
               chisq.max = chisq.max)
  if (backend == "R") {
    return(do.call(.ldsc_legacy, args))
  }
  out <- do.call(.ldsc_fast, args)
  if (is.character(out)) {
    message("ldsc: using the R implementation (", out, ")")
    return(do.call(.ldsc_legacy, args))
  }
  out
}

.ldsc_workers <- function() {
  threads <- as.integer(getOption("GenomicSEM.ldsc_threads", parallel::detectCores()))
  if (.Platform$OS.type == "windows" || is.na(threads) || threads < 1L) 1L else threads
}

# Forked workers (one item each) on Unix; serial on Windows.
.ldsc_lapply <- function(X, FUN) {
  workers <- min(.ldsc_workers(), length(X))
  if (workers <= 1L) return(lapply(X, FUN))
  # Each forked worker reads one file; fread's own threads would oversubscribe.
  out <- parallel::mclapply(X, function(item) {
    options(GenomicSEM.ldsc_fread_threads = 1L)
    FUN(item)
  }, mc.cores = workers, mc.preschedule = FALSE)
  failed <- vapply(out, inherits, logical(1), "try-error")
  if (any(failed)) stop(out[[which(failed)[1]]], call. = FALSE)
  out
}

.ldsc_read_ld <- function(dir, chroms, suffix) {
  out <- do.call(rbind, .ldsc_lapply(chroms, function(i) {
    data.table::fread(file.path(dir, paste0(i, suffix)), showProgress = FALSE, data.table = FALSE,
                      nThread = getOption("GenomicSEM.ldsc_fread_threads", data.table::getDTthreads()))
  }))
  out$CM <- NULL
  out$MAF <- NULL
  out
}

# One trait: read, drop incomplete rows, join to the LD scores as
# merge(y1, w) then merge(., x) do with sort = FALSE (y1's rows present in
# both, in y1's order), order by CHR, BP (stable) and apply the chi-square
# filter. Returns the kernel's columns and the counts the log reports.
.ldsc_prepare_trait <- function(chi1, x, w, chisq.max) {
  y1 <- data.table::fread(chi1, sep = "\t", na.strings = c("", "NA"), showProgress = FALSE,
                          colClasses = list(character = c("SNP", "A1")), data.table = FALSE,
                          nThread = getOption("GenomicSEM.ldsc_fread_threads", data.table::getDTthreads()))
  y1 <- na.omit(y1)
  if (anyDuplicated(y1$SNP)) return(list(duplicated = TRUE))
  iw <- match(y1$SNP, w$SNP)
  ix <- match(y1$SNP, x$SNP)
  keep <- which(!is.na(iw) & !is.na(ix))
  ord <- keep[order(x$CHR[ix[keep]], x$BP[ix[keep]])]
  iw <- iw[ord]
  ix <- ix[ord]
  N <- as.double(y1$N[ord])
  Z <- as.double(y1$Z[ord])
  trait.chisq.max <- if (is.na(chisq.max)) max(0.001 * max(N), 80) else chisq.max
  ok <- !(Z^2 > trait.chisq.max)
  cols <- list(N = N[ok], Z = Z[ok], A1 = y1$A1[ord][ok], wLD = w$L2[iw][ok], L2 = x$L2[ix][ok],
               id = ix[ok] - 1L)
  cols$complete <- as.integer(complete.cases(as.data.frame(cols), x$CHR[ix][ok], x$BP[ix][ok]))
  c(cols, list(duplicated = FALSE, n.read = nrow(y1), n.merged = length(ord), n.removed = sum(!ok),
               chisq.max = trait.chisq.max, median.chi = median(cols$Z^2)))
}

# Returns the ldsc() result, or a string giving why the reference R
# implementation must be used instead.
.ldsc_fast <- function(traits, sample.prev, population.prev, ld, wld,
                       trait.names = NULL, sep_weights = FALSE, chr = 22,
                       n.blocks = 200, ldsc.log = NULL, stand = FALSE, select = FALSE, chisq.max = NA) {

  begin.time <- Sys.time()

  if(is.null(ldsc.log)){
    logtraits<-gsub(".*/","",traits)
    log2<-paste(logtraits,collapse="_")
    if(object.size(log2) > 200){
      log2<-substr(log2,1,100)
    }
    log.file <- file(paste0(log2, "_ldsc.log"),open="wt")
  }else{log.file<-file(paste0(ldsc.log, "_ldsc.log"),open="wt")}
  fallback <- function(reason) { close(log.file); reason }

  if(!isFALSE(select) && !identical(select, "ODD") && !identical(select, "EVEN") && !is.numeric(select)){
    close(log.file)
    stop("select must be one of the following values: FALSE, 'ODD', 'EVEN', or a numeric vector of chromosome numbers.")
  }

  .LOG("Multivariate ld-score regression of ", length(traits), " traits ", "(", paste(traits, collapse = " "), ")", " began at: ", begin.time, file=log.file)

  chroms <- if (isFALSE(select)) 1:chr else if (identical(select, "ODD")) seq(1, chr, 2) else
    if (identical(select, "EVEN")) seq(2, chr, 2) else select

  # Dimensions
  n.traits <- length(traits)
  n.V <- n.traits * (n.traits + 1) / 2

  if(n.traits > 18){
    n.blocks<-(((n.traits+1)*(n.traits+2))/2)+1
    .LOG("     ", file=log.file, print = FALSE)
    .LOG("Setting the number of blocks used to perform the block jacknife used to estimate the sampling covariance matrix (V) to ", n.blocks, file=log.file)
    .LOG("This reflects the need to estimate V using at least one more block then their are nonredundant elements in the genetic covariance matrix that includes individual SNPs.", file=log.file)
    .LOG("If the n.blocks is > 1000 you should carefully inspect output for any strange results, such as extremely significant Q_SNP estimates.", file=log.file)
    .LOG("     ", file=log.file, print = FALSE)
    if(n.blocks > 1000){
      warning("The number of blocks needed to estimate V is > 1000, which may result in sampling dependencies across the blocks used to estimate standard errors and can bias results.")
    }
  }

  if(!(is.null(trait.names))){
    check_names<-str_detect(trait.names, "[+\\-*/^]")
    if(any(check_names))
      warning("Your trait names specified include mathematical arguments (e.g., +, -, *, /) that will be misread by lavaan. Please rename the traits using the trait.names argument.")
  }

  if(length(traits)==1)
    warning("Our version of ldsc requires 2 or more traits. Please include an additional trait.")

  # Storage:
  cov <- matrix(NA,nrow=n.traits,ncol=n.traits)
  V.hold <- matrix(NA,nrow=n.blocks,ncol=n.V)
  N.vec <- matrix(NA,nrow=1,ncol=n.V)
  Liab.S <- rep(1, n.traits)
  I <- matrix(NA,nrow=n.traits,ncol=n.traits)

  #########  READ LD SCORES:
  .LOG("Reading in LD scores", file=log.file)
  x <- .ldsc_read_ld(ld, chroms, ".l2.ldscore.gz")
  w <- if (sep_weights) .ldsc_read_ld(wld, chroms, ".l2.ldscore.gz") else x
  ld.cols <- c("CHR", "SNP", "BP", "L2")
  if (!setequal(names(x), ld.cols) || !setequal(names(w), ld.cols) || names(w)[ncol(w)] != "L2") {
    return(fallback("LD-score files with columns other than CHR, SNP, BP, CM, MAF, L2"))
  }
  if (anyDuplicated(x$SNP) || anyDuplicated(w$SNP)) {
    return(fallback("duplicated SNP IDs in the LD-score files"))
  }

  ### READ M
  m <- sum(unlist(lapply(chroms, function(i) scan(file.path(ld, paste0(i, ".l2.M_5_50")), quiet = TRUE))))
  M.tot <- m

  ### READ ALL CHI2 + MERGE WITH LDSC FILES
  all_y <- .ldsc_lapply(traits, function(chi1) .ldsc_prepare_trait(chi1, x, w, chisq.max))
  for (s in seq_len(n.traits)) {
    chi1 <- traits[s]
    y <- all_y[[s]]
    .LOG("Read in summary statistics [", s, "/", n.traits, "] from: ", chi1, file=log.file)
    if (y$duplicated) {
      return(fallback(paste("duplicated SNP IDs in", chi1)))
    }
    .LOG("Out of ", y$n.read, " SNPs, ", y$n.merged, " remain after merging with LD-score files", file=log.file)
    .LOG("Removing ", y$n.removed, " SNPs with Chi^2 > ", y$chisq.max, "; ", length(y$Z), " remain", file=log.file)
  }

  # Pairs in the order of the reference implementation's loops.
  pairs <- do.call(rbind, lapply(seq_len(n.traits), function(j) cbind(j, j:n.traits)))
  alleles <- unique(unlist(lapply(all_y, `[[`, "A1")))
  rows <- vapply(all_y, function(y) length(y$Z), integer(1))
  col <- function(name) unlist(lapply(all_y, `[[`, name), use.names = FALSE)
  res <- .Call(.genomicsem_ldsc,
               as.double(col("L2")), as.double(col("wLD")), col("N"), col("Z"),
               match(col("A1"), alleles), as.integer(col("id")),
               col("complete"),
               as.integer(c(0L, cumsum(rows))), as.integer(nrow(x)),
               as.integer(pairs[, 1] - 1L), as.integer(pairs[, 2] - 1L),
               as.double(M.tot), as.integer(n.blocks),
               .ldsc_workers())
  if (any(res$status != 0L)) {
    return(fallback(if (any(res$status == 2L)) "fewer SNPs than jackknife blocks" else "a singular regression"))
  }

  for (s in seq_len(n.V)) {
    j <- pairs[s, 1]; k <- pairs[s, 2]
    chi1 <- traits[j]
    reg.tot <- res$reg_tot[s]
    intercept <- res$intercept[s]
    intercept.se <- res$intercept_se[s]
    tot.se <- res$tot_se[s]
    V.hold[, s] <- res$pseudo[, s]
    N.vec[1, s] <- res$n_bar[s]

    if (j == k) {
      .LOG("     ", "     ", file=log.file, print = FALSE)
      .LOG("Estimating heritability [", s, "/", n.V, "] for: ", chi1, file=log.file)
      samp.prev <- sample.prev[j]
      pop.prev <- population.prev[j]
      if(is.na(pop.prev)==F & is.na(samp.prev)==F){
        conversion.factor <- (pop.prev^2*(1-pop.prev)^2)/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
        Liab.S[j] <- conversion.factor
        .LOG("     ", file=log.file, print = FALSE)
        .LOG("Please note that the results initially printed to the screen and log file reflect the NON-liability h2 and cov_g. However, a liability conversion is being used for trait ",
            chi1, " when creating the genetic covariance matrix used as input for Genomic SEM and liability scale results are printed at the end of the log file.", file=log.file)
        .LOG("     ", file=log.file, print = FALSE)
      }
      cov[j,j] <- reg.tot
      I[j,j] <- intercept

      lambda.gc <- all_y[[j]]$median.chi / qchisq(0.5, df = 1)
      mean.Chi <- res$mean_y[s]
      ratio <- (intercept - 1) / (mean.Chi - 1)
      ratio.se <- intercept.se / (mean.Chi - 1)

      .LOG("Heritability Results for trait: ", chi1, file=log.file)
      .LOG("Mean Chi^2 across remaining SNPs: ", round(mean.Chi, 4), file=log.file)
      .LOG("Lambda GC: ", round(lambda.gc, 4), file=log.file)
      .LOG("Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")", file=log.file)
      .LOG("Ratio: ", round(ratio, 4), " (", round(ratio.se, 4), ")", file=log.file)
      .LOG("Total Observed Scale h2: ", round(reg.tot, 4), " (", round(tot.se, 4), ")", file=log.file)
      .LOG("h2 Z: ", format(reg.tot / tot.se, digits = 3), file=log.file)
    } else {
      .LOG("     ", file=log.file, print = FALSE)
      chi2 <- traits[k]
      .LOG("Calculating genetic covariance [", s, "/", n.V, "] for traits: ", chi1, " and ", chi2, file=log.file)
      .LOG(res$n_snps[s], " SNPs remain after merging ", chi1, " and ", chi2, " summary statistics", file=log.file)
      cov[k, j] <- cov[j, k] <- reg.tot
      I[k, j] <- I[j, k] <- intercept
      .LOG("Results for genetic covariance between: ", chi1, " and ", chi2, file=log.file)
      .LOG("Mean Z*Z: ", round(res$mean_y[s], 4), file=log.file)
      .LOG("Cross trait Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")", file=log.file)
      .LOG("Total Observed Scale Genetic Covariance (g_cov): ", round(reg.tot, 4), " (", round(tot.se, 4), ")", file=log.file)
      .LOG("g_cov Z: ", format(reg.tot / tot.se, digits = 3), file=log.file)
      .LOG("g_cov P-value: ", format(2 * pnorm(abs(reg.tot / tot.se), lower.tail = FALSE), digits = 5), file=log.file)
    }
  }
  ## Scale V to N per study (assume m constant)
  # /!\ crossprod instead of tcrossprod because N.vec is a one-row matrix
  v.out <- cov(V.hold) / crossprod(N.vec * (sqrt(n.blocks) / m))

  ### Scale S and V to liability:
  ratio <- tcrossprod(sqrt(Liab.S))
  S <- cov * ratio

  #calculate the ratio of the rescaled and original S matrices
  scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)

  #rescale the sampling correlation matrix by the appropriate diagonals
  V <- v.out * tcrossprod(scaleO)


  #name traits according to trait.names argument
  #use general format of V1-VX if no names provided
  colnames(S) <- if (is.null(trait.names)) paste0("V", 1:ncol(S)) else trait.names

  if(mean(Liab.S)!=1){
    r<-nrow(S)
    SE<-matrix(0, r, r)
    SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(V))

    .LOG(c("     ", "     "), file=log.file, print = FALSE)
    .LOG("Liability Scale Results", file=log.file)

    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j == k){
          .LOG("     ", file=log.file, print = FALSE)
          .LOG("Liability scale results for: ", chi1, file=log.file)
          .LOG("Total Liability Scale h2: ", round(S[j, j], 4), " (", round(SE[j, j], 4), ")", file=log.file)
        }

        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          .LOG("Total Liability Scale Genetic Covariance between ", chi1, " and ",
              chi2, ": ", round(S[k, j], 4), " (", round(SE[k, j], 4), ")", file=log.file)
          .LOG("     ", file=log.file, print = FALSE)
        }
      }
    }
  }


  if(all(diag(S) > 0)){

    ##calculate standardized results to print genetic correlations to log and screen
    ratio <- tcrossprod(1 / sqrt(diag(S)))
    S_Stand <- S * ratio

    #calculate the ratio of the rescaled and original S matrices
    scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)

    ## MAke sure that if ratio in NaN (devision by zero) we put the zero back in
    # -> not possible because of 'all(diag(S) > 0)'
    # scaleO[is.nan(scaleO)] <- 0

    #rescale the sampling correlation matrix by the appropriate diagonals
    V_Stand <- V * tcrossprod(scaleO)

    #enter SEs from diagonal of standardized V
    r<-nrow(S)
    SE_Stand<-matrix(0, r, r)
    SE_Stand[lower.tri(SE_Stand,diag=TRUE)] <-sqrt(diag(V_Stand))


    .LOG(c("     ", "     "), file=log.file, print = FALSE)
    .LOG("Genetic Correlation Results", file=log.file)

    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          .LOG("Genetic Correlation between ", chi1, " and ", chi2, ": ",
              round(S_Stand[k, j], 4), " (", round(SE_Stand[k, j], 4), ")", file=log.file)
          .LOG("     ", file=log.file, print = FALSE)
        }
      }
    }
  }else{
    warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
    .LOG("Your genetic covariance matrix includes traits estimated to have a negative heritability.", file=log.file, print = FALSE)
    .LOG("Genetic correlation results could not be computed due to negative heritability estimates.", file=log.file)
  }

  end.time <- Sys.time()

  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- floor(total.time-mins*60)

  .LOG("     ", file=log.file, print = FALSE)
  .LOG("LDSC finished running at ", end.time, file=log.file)
  .LOG("Running LDSC for all files took ", mins, " minutes and ", secs, " seconds", file=log.file)
  .LOG("     ", file=log.file, print = FALSE)

  flush(log.file)
  close(log.file)

  if(stand){
    list(V=V,S=S,I=I,N=N.vec,m=m,V_Stand=V_Stand,S_Stand=S_Stand)
  } else {
    list(V=V,S=S,I=I,N=N.vec,m=m)
  }
}

