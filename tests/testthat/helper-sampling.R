library(checkmate)
library(assertthat)
library(Formula)
library(abind)
library(rstan)

set_sampling_default <- function(iter, warmup, chains, cores = getOption("mc.cores", 1), save_warmup = FALSE, backend = "rstan", control=list()) {
  options(
    OncoBayes2.MC.iter = iter,
    OncoBayes2.MC.warmup = warmup,
    OncoBayes2.MC.chains = chains,
    mc.cores = cores,
    OncoBayes2.MC.save_warmup = save_warmup,
    OncoBayes2.MC.backend = backend,
    OncoBayes2.MC.control = control
  )
}

very_fast_sampling <- function() {
  message("Tests running with very fast sampling")
  ## note: 250 warmups are needed to get Stans NUTS adaptation to work
  set_sampling_default(500, 250, 1, 1, control=list(adapt_delta=0.85))
}

fake_sampling <- function() {
  message("Tests running with fake sampling")
  set_sampling_default(10, 2, 1, 1, control=list(adapt_delta=0.85))
}

default_sampling <- function() {
  set_sampling_default(2000, 1000, 4, 1)
}

fast_sampling <- function() {
  message("Tests running with fast sampling")
  set_sampling_default(1000, 500, 1, 1, control=list(adapt_delta=0.85))
}

run_example <- function(example) {
  env <- new.env()
  suppressWarnings(example_model(example, env, silent = TRUE))
  invisible(env)
}



## set up slim sampling in case we are on CRAN
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  very_fast_sampling()
} else {
  fake_sampling()
}


## takes a flat array of a mixture mvn and returns an RBesT mixture
## mvn
array2mix <- function(a, p) {
  Nc <- dim(a)[1]
  l <- dim(a)[2]
  mix <- list()
  for (i in seq_len(Nc)) {
    w <- a[i, 1]
    m <- a[i, 2:(2 + p - 1)]
    s <- a[i, (2 + p):(2 + p + p - 1)]
    if (p > 1) {
      rho <- a[i, (2 + p + p):l]
      R <- diag(1, p, p)
      R[lower.tri(R)] <- rho
      R[upper.tri(R)] <- t(R)[upper.tri(R)]
    } else {
      R <- diag(1, 1, 1)
    }
    mix[[i]] <- unname(c(w, m, diag(s, nrow = p) %*% R %*% diag(s, nrow = p)))
  }
  do.call(RBesT::mixmvnorm, mix)
}

## takes a fitted blrmfit object and returns the blrmfit object with a
## posterior containing one draw with all slopes and intercept set to
## the prior mean (exactly).
sample_prior_mean <- function(blrmfit) {
  ps <- prior_summary(blrmfit)

  num_groups <- ps$num_groups
  num_strata <- ps$num_strata
  num_comp <- dim(ps$EX_mu_log_beta)[2]
  has_inter <- ps$has_inter
  num_inter <- 0
  if (has_inter) {
    num_inter <- length(grep("^m", dimnames(ps$EX_mu_eta)$prior))
  }

  draw <- list(
    log_beta_raw = array(0, c(2 * num_groups, num_comp, 2)),
    eta_raw = array(0, c(2 * num_groups, num_inter)),
    mu_log_beta = array(0, c(num_comp, 2)),
    tau_log_beta_raw = array(1, c(num_strata, num_comp, 2)),
    L_corr_log_beta = abind(replicate(num_comp, diag(2), FALSE), along = -1),
    mu_eta = array(0, c(num_inter)),
    tau_eta_raw = array(1, c(num_strata, num_inter)),
    L_corr_eta = diag(num_inter)
  )
  
  EX_m_intercept <- apply(ps$EX_mu_log_beta["weight", , , drop=FALSE] *
                            ps$EX_mu_log_beta["m_intercept", , , drop=FALSE],
                          c(1,2), sum)
  EX_m_log_slope <- apply(ps$EX_mu_log_beta["weight", , , drop=FALSE] *
                            ps$EX_mu_log_beta["m_log_slope", , , drop=FALSE],
                          c(1,2), sum)

  NEX_m_intercept <- apply(ps$NEX_mu_log_beta["weight", , , drop=FALSE] *
                            ps$NEX_mu_log_beta["m_intercept", , , drop=FALSE],
                          c(1,2), sum)
  NEX_m_log_slope <- apply(ps$NEX_mu_log_beta["weight", , , drop=FALSE] *
                            ps$NEX_mu_log_beta["m_log_slope", , , drop=FALSE],
                          c(1,2), sum)

  draw$log_beta_raw <- abind(c(replicate(num_groups,
                                         abind(EX_m_intercept,
                                               EX_m_log_slope, along=3),
                                       simplify=FALSE),
                               replicate(num_groups,
                                         abind(NEX_m_intercept,
                                               NEX_m_log_slope, along=3),
                                         simplify=FALSE))
                               , along=1)

  if (has_inter) {
    EX_mu_eta <- array(summary(array2mix(ps$EX_mu_eta, num_inter))$mean, num_inter)
    NEX_mu_eta <- array(summary(array2mix(ps$NEX_mu_eta, num_inter))$mean, num_inter)
    draw$eta_raw <- abind(c(
      replicate(num_groups, EX_mu_eta, simplify=FALSE),
      replicate(num_groups, NEX_mu_eta, simplify=FALSE)),
      along=-1)
  }
  
  msg <- capture.output(blrmfit$stanfit <- sampling(OncoBayes2:::stanmodels$blrm_exnex, data = blrmfit$standata, chains = 1, iter = 1, warmup = 0, seed = 23542, init = list(draw), algorithm = "Fixed_param", open_progress = FALSE))
  blrmfit
}
