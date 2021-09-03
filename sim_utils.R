# Extract subject indices -----------------------------------------------------

#' Extract subject indices from tree sequence
#'
#' This is an internal function used to extract subject indices from a tskit
#' tree sequence.

get_subjects <- function(ts){

  if (class(ts)[1] != "tskit.trees.TreeSequence")
    stop("input must be a tskit TreeSequence")

  #TODO: extend to handle different numbers of samples per population
  nPop <- ts$num_populations

  # Haploid indices
  indHap <- lapply(seq(0,nPop-1), ts$samples)

  # For subjects, take even indices only, with one-indexing
  lapply(indHap, function(x) as.integer(1 + x[!(x %% 2)]/2))
}

# Filter tree sequence on subjects --------------------------------------------

#' Filter tree sequence by subjects
#'
#' This is an internal function used to filter a tskit treesequence by subjects.

filter_sub <- function(ts, indSub){

  if (class(ts)[1] != "tskit.trees.TreeSequence")
    stop("input must be a tskit TreeSequence")

  nPop <- ts$num_populations
  if (is.list(indSub)){
    if (length(indSub) > nPop) {
      stop("indSub must be a single vector of subject indices
           or a list of indices for each population")
    }
    indSub <- unlist(indSub)
  }

  # Get haploid indices for each subject
  indHap <- as.integer(sort(c(2*(indSub-1), 2*(indSub-1) + 1)))

  ts$simplify(indHap, filter_sites=FALSE)
}


# Compute proportion of variance explained ------------------------------------

#' Compute r2
#'
#' This is a utility function for computing the r2 (proportion of variance
#' explained) of a prediction of a continuous phenotype.

compute_r2 <- function(y, pred){

  ssTot <- sum((y - mean(y))^2)
  ssRes <- sum((y - pred)^2)

  1 - (ssRes/ssTot)
}


## Simulate true effects ------------------------------------------------------

#' Simulate effect sizes
#'
#' simul_effects() generates effect sizes for each of the variants along a
#' tree sequence. Causal variants (i.e. non-zero effects) are assumed to be
#' evenly spaced along the genome. Non-zero effect sizes are independent and
#' normally distributed with mean zero and total variance 1.

simul_effects <- function(p,
                          n_causal = p,
                          beta_cov = matrix(1),
                          ind_causal = NULL){

  if (n_causal > p | n_causal < 0)
    stop("n_causal must be an integer between 0 and nSNP")

  if (is.null(ind_causal))
    ind_causal <- seq(sample(seq(1, p/n_causal), 1), p, length.out = n_causal)

  n_pop <- nrow(beta_cov)
  true_beta <- matrix(0, p, n_pop)
  rownames(true_beta) <- seq(1, p)
  true_beta[ind_causal, ] <- mvtnorm::rmvnorm(n_causal, sigma = beta_cov)

  true_beta
}


# Calculate genetic risk ------------------------------------------------------

#' Calculate genetic risk
#'
#' This function computes the genetic risk for individuals in the tree sequence
#' by multiplying their genotype with the true variant effect sizes.

get_snp_risk <- function(ts, true_beta, ind_snp = ts$num_sites, ind_sub = NULL){

  # Filter subjects
  if (!is.null(ind_sub))
    ts <- filter_sub(ts, ind_sub)

  if (length(true_beta) != length(ind_snp))
    stop("Length of true_beta must be equal to length of ind_snp")

  # Remove non-causal SNPs
  ind_causal <- ind_snp[true_beta != 0]
  ind_delete <- as.integer(setdiff(1:ts$num_sites, ind_causal) - 1)
  g <- ts$delete_sites(ind_delete)$genotype_matrix()

  # Get diploid samples from haploids
  n <- ncol(g)
  g <- g[ ,seq(1, n-1, 2)] + g[ ,seq(2, n, 2)]

  drop(true_beta[true_beta != 0] %*% g)
}


# Simulate disease ------------------------------------------------------------

#' Simulate phenotypes
#'
#' simul_pheno() generates continuous or binary phenotypes for each individual
#' based on their genetic risk.

simul_pheno <- function(snp_risk,
                        env_var = var(snp_risk) * ((1 / h2) - 1),
                        h2 = 1 - env_var / (var(snp_risk) + env_var),
                        center = FALSE,
                        response = "cts",
                        prevalence = 0.1) {

  if (h2 < 0 | h2 > 1)
    stop("h2 should be a number between 0 and 1")

  response <- match.arg(response)
  n <- length(snp_risk)

  env_risk <- rnorm(n, sd = sqrt(env_var))

  total_risk <- scale(snp_risk + env_risk, center, scale = FALSE)

  if (response == "bin") {
    func <- function(y) {
      (mean(1 / (1 + exp(-snp_risk - y))) - prevalence)^2
    }
    intercept <- optimise(func, c(-100, 100))$minimum
    eta <- 1 / (1 + exp(-snp_risk - intercept))

    pheno <- as.integer(runif(n) < eta)
  } else
    pheno <- total_risk

  pheno
}


# Extract genotype matrix -----------------------------------------------------

#' Extract genotype matrix from tree sequence
#'
#' get_genotype() extracts a genotype matrix from a tskit tree sequence after
#' applying any individual or variant (i.e.) SNP filters. The function assumes
#' that individuals are stored as haploid pairs and combines them into diploid
#' samples so that each entry is 0, 1, or 2.

get_genotype <- function(ts,
                         ind_sub = seq(1, ts$num_samples),
                         ind_snp = NULL,
                         sparse = TRUE){

  # Filter SNPs
  p <- ts$num_sites
  ind_delete <- as.integer(setdiff(1:p, ind_snp) - 1)
  ts <- ts$delete_sites(ind_delete)

  # Filter subjects
  ts <- filter_sub(ts, ind_sub)

  # Extract (diploid) genotype matrix from tree sequence
  G <- ts$genotype_matrix()
  n <- ncol(G)
  G <- t(G[ ,seq(1, n-1, 2)] + G[ ,seq(2, n, 2)])

  rownames(G) <- unlist(ind_sub)
  colnames(G) <- ind_snp

  if (sparse)
    G <- Matrix::Matrix(G)

  G
}