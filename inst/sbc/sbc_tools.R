#'
#' Utilities for SBC validation
#'

load_OB2_dev <- function(ob2_lib_dir) {
  if (ob2_lib_dir != .libPaths()[1]) {
    if ("OncoBayes2" %in% .packages()) {
      warning("Unloading non-development OncoBayes2")
    }
    unloadNamespace("OncoBayes2")
    cat("Setting libPaths to dev OncoBayes2 install.\n")
    .libPaths(c(ob2_lib_dir, .libPaths()))
  }

  if (!("OncoBayes2" %in% .packages())) {
    library("OncoBayes2", lib.loc = ob2_lib_dir, character.only = TRUE)
  }
}

setup_lecuyer_seeds <- function(lecuyer_seed, num) {
  ## note: seed have the format from L'Ecuyer. Just set
  ## RNGkind("L'Ecuyer-CMRG")
  ## and then use .Random.seed
  job_seeds <- list()
  job_seeds[[1]] <- parallel::nextRNGStream(lecuyer_seed)
  i <- 2
  while (i < num + 1) {
    job_seeds[[i]] <- parallel::nextRNGStream(job_seeds[[i - 1]])
    i <- i + 1
  }
  job_seeds
}

stan_adaptation_phases <- function(warmup = 1000, init = 75, window = 25, term = 50) {
  iter <- 0
  iter <- iter + init
  total_learn_cov <- warmup - iter - term
  phases <- list(init = init, mass_matrix = c(), term = term)
  i <- 1
  cur_window <- window
  while (total_learn_cov - cur_window > 0) {
    phases$mass_matrix <- c(phases$mass_matrix, cur_window)
    total_learn_cov <- total_learn_cov - cur_window
    cur_window <- 2 * cur_window
    i <- i + 1
  }
  if (total_learn_cov > 0) {
    phases$mass_matrix <- c(phases$mass_matrix, total_learn_cov)
  }
  ## phases$mass_matrix[i-1] <- phases$mass_matrix[i-1] - sum(c(init, term, phases$mass_matrix[- (i-1)]))
  phases
}

## Sample prior

## First sample EX hyperparameters
## prior_EX_mu_mean_comp + prior_EX_mu_sd_comp => EX_mu_comp (group means EX)
## prior_EX_tau_mean_comp + prior_EX_tau_sd_comp => EX_tau_comp (group taus EX)
## prior_EX_corr_eta_comp => rho (only 0 for now)

## prior_EX_mu_mean_inter + prior_EX_mu_sd_inter => EX_mu_inter (group means EX)
## prior_EX_tau_mean_inter + prior_EX_tau_sd_inter => EX_tau_comp (group taus EX)
## prior_EX_corr_eta_inter => rho (only 0 for now)
## => EX_eta (group specific means for inter (mvn normal))

## Then sample group parameters, EX
## => log_beta (group specific means for comp (mvn normal))
## => eta (group specific means for comp (mvn normal))

## and NEX parameters are in the same data structure with index starting at num_groups+1

## prior_NEX_mu_mean_comp + prior_NEX_mu_sd_comp => NEX_beta = NEX_comp (group means NEX iid groupwise)
## prior_NEX_mu_mean_inter + prior_NEX_mu_sd_inter => NEX_eta = NEX_inter (group means NEX iid groupwise)

## prior_is_EXNEX_comp + prior_EX_prob_comp => pick which one
## prior_is_EXNEX_inter + prior_EX_prob_inter => pick which one

sample_prior <- function(model) {
  prior <- prior_summary(model$base_fit)
  num_strata <- prior$num_strata
  num_groups <- prior$num_groups
  num_comp <- dim(prior$EX_mu_log_beta)[2]
  has_inter <- prior$has_inter
  num_inter <- 0
  if (has_inter) {
    num_inter <- dim(prior$EX_prob_inter)[2]
  } else {
    num_inter <- 0
  }
  blrm_args <- model$blrm_args
  standata <- model$base_fit$standata
  group_stratum <- standata$group_stratum_cid

  ## group specific parameters: EX, then NEX
  log_beta <- array(NA, dim = c(2 * num_groups, num_comp, 2))
  eta <- array(NA, dim = c(2 * num_groups, num_inter))

  ## sample EX hyperparameters

  ## mu
  EX_mu_comp <- array(NA, dim = c(num_comp, 2))
  EX_mu_inter <- array(NA, dim = c(num_inter))

  comp2mix <- function(prior_comp, comp) {
    max_Nc <- dim(prior_comp)[3]
    mix <- list()
    for (i in seq_len(max_Nc)) {
      if (prior_comp["weight", comp, i] != 0) {
        s1 <- prior_comp[4, comp, i]
        s2 <- prior_comp[5, comp, i]
        rho <- prior_comp["rho", comp, i]
        mix[[i]] <- unname(c(prior_comp[1, comp, i], prior_comp[2:3, comp, i], matrix(c(s1^2, s1 * s2 * rho, s1 * s2 * rho, s2^2), 2, 2)))
      }
    }
    do.call(RBesT::mixmvnorm, mix)
  }
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

  for (j in seq_len(num_comp)) {
    EX_mu_comp[j, 1:2] <- RBesT::rmix(comp2mix(prior$EX_mu_log_beta, j), 1)
  }

  if (has_inter) {
    EX_mu_inter[1:num_inter] <- RBesT::rmix(array2mix(prior$EX_mu_eta, num_inter), 1)
  }

  ## tau
  EX_tau_comp <- array(NA, dim = c(num_strata, num_comp, 2))
  EX_tau_inter <- array(NA, dim = c(num_strata, num_inter))
  ## correlation matrix
  EX_corr_comp <- array(NA, dim = c(num_strata, num_comp, 2, 2))
  EX_corr_inter <- array(NA, dim = c(num_strata, num_inter, num_inter))
  EX_Sigma_comp <- array(NA, dim = c(num_strata, num_comp, 2, 2))
  EX_Sigma_inter <- array(NA, dim = c(num_strata, num_inter, num_inter))

  sample_tau_prior <- function(dist, mix) {
    if (dist == 0) {
      m <- mix[2:3, , drop = FALSE]
      w <- mix[1, ]
      wm <- sweep(m, 2, w, "*")
      return(rowSums(wm))
    }
    if (dist == 1) {
      return(exp(RBesT::rmix(mix, 1)))
    }
    if (dist == 2) {
      return(abs(RBesT::rmix(mix, 1)))
    }
    stop("Unsupported tau prior density.")
  }

  for (s in seq_len(num_strata)) {
    for (j in seq_len(num_comp)) {
      mix_tau_comp <- comp2mix(adrop(prior$EX_tau_log_beta[s, , , , drop = FALSE], 1), j)
      EX_tau_comp[s, j, 1:2] <- sample_tau_prior(prior$tau_dist, mix_tau_comp)
      EX_corr_comp[s, j, , ] <- rcorvine(2, standata$prior_EX_corr_eta_comp[j], FALSE)
      EX_Sigma_comp[s, j, , ] <- diag(as.vector(EX_tau_comp[s, j, ]), 2, 2) %*% matrix(EX_corr_comp[s, j, , , drop = FALSE], 2, 2) %*% diag(as.vector(EX_tau_comp[s, j, ]), 2, 2)
    }

    if (has_inter) {
      mix_tau_inter <- array2mix(adrop(prior$EX_tau_eta[s, , , drop = FALSE], 1), num_inter)
      EX_tau_inter[s, 1:num_inter] <- sample_tau_prior(prior$tau_dist, mix_tau_inter)
    }
    EX_corr_inter[s, , ] <- diag(num_inter)
    if (num_inter > 1) {
      EX_corr_inter[s, , ] <- rcorvine(num_inter, standata$prior_EX_corr_eta_inter, FALSE)
    }
    EX_Sigma_inter[s, , ] <- diag(as.vector(EX_tau_inter[s, , drop = FALSE]), num_inter, num_inter) %*% matrix(EX_corr_inter[s, , , drop = FALSE], num_inter, num_inter) %*% diag(as.vector(EX_tau_inter[s, , drop = FALSE]), num_inter, num_inter)
  }

  ## EX - group-specific parameters
  for (g in seq_len(num_groups)) {
    s <- group_stratum[g]
    for (j in seq_len(num_comp)) {
      log_beta[g, j, 1:2] <- rmvnorm(1, EX_mu_comp[j, ], EX_Sigma_comp[s, j, , ])
      ## log_beta[g,j,1] <- rnorm(1, EX_mu_comp[j,1], EX_tau_comp[s,j,1] )
      ## log_beta[g,j,2] <- rnorm(1, EX_mu_comp[j,2], EX_tau_comp[s,j,2] )
      ## assert_that(standata$prior_EX_corr_eta_comp[j] == 1, msg="LKJ correlation == 1 is only supported.")
    }
    if (num_inter > 0) {
      if (num_inter > 1) {
        eta[g, ] <- rmvnorm(1, EX_mu_inter, EX_Sigma_inter[s, , ])
      } else {
        eta[g, 1] <- rnorm(1, EX_mu_inter, EX_tau_inter[s, 1])
      }
      ## for (j in 1:num_inter) {
      ##  eta[g,j] <- rnorm(1, EX_mu_inter[j], EX_tau_inter[s,j] )
      ## }
      ## assert_that(standata$prior_EX_corr_eta_inter == 1, msg="LKJ correlation == 1 is only supported.")
    }
  }
  ## NEX - group-specific parameters
  for (g in seq_len(num_groups)) {
    for (j in seq_len(num_comp)) {
      log_beta[num_groups + g, j, 1:2] <- RBesT::rmix(comp2mix(prior$NEX_mu_log_beta, j), 1)
    }

    if (has_inter) {
      eta[num_groups + g, 1:num_inter] <- adrop(RBesT::rmix(array2mix(prior$NEX_mu_eta, num_inter), 1), 1)
    }
  }

  ## convert slope to natural scale (enforced positivity)
  beta <- log_beta
  for (g in seq_len(2 * num_groups)) {
    for (j in seq_len(num_comp)) {
      beta[g, j, 2] <- exp(beta[g, j, 2])
    }
  }

  ## sample EX / NEX membership
  is_EX_comp <- array(NA, dim = c(num_groups, num_comp))
  is_EX_inter <- array(NA, dim = c(num_groups, num_inter))
  draw_beta <- array(NA, dim = c(num_groups, num_comp, 2))
  draw_eta <- array(NA, dim = c(num_groups, num_inter))
  for (g in seq_len(num_groups)) {
    for (j in seq_len(num_comp)) {
      if (standata$prior_is_EXNEX_comp[j] == 1) {
        is_EX_comp[g, j] <- rbinom(1, 1, standata$prior_EX_prob_comp[g, j])
      } else {
        is_EX_comp[g, j] <- 1
      }
      gidx <- ifelse(is_EX_comp[g, j] == 1, g, num_groups + g)
      draw_beta[g, j, 1] <- beta[gidx, j, 1]
      draw_beta[g, j, 2] <- beta[gidx, j, 2]
    }

    for (j in seq_len(num_inter)) {
      if (standata$prior_is_EXNEX_inter[j] == 1) {
        is_EX_inter[g, j] <- rbinom(1, 1, standata$prior_EX_prob_inter[g, j])
      } else {
        is_EX_inter[g, j] <- 1
      }
      gidx <- ifelse(is_EX_inter[g, j] == 1, g, num_groups + g)
      draw_eta[g, j] <- eta[gidx, j]
    }
  }

  ## name the array indices accordingly using the prior_summary
  ## structures
  ps <- prior_summary(model$base_fit)
  dimnames(EX_mu_comp) <- c(dimnames(ps$EX_mu_log_beta)[2], list(c("intercept", "log_slope")))
  dimnames(EX_tau_comp) <- c(dimnames(ps$EX_tau_log_beta)[c(1, 3)], list(c("tau_intercept", "tau_log_slope")))

  if (has_inter) {
    dimnames(EX_mu_inter) <- list(dimnames(ps$EX_mu_eta)[[2]][2:(2 + num_inter - 1)])
    dimnames(EX_tau_inter) <- c(dimnames(ps$EX_tau_eta)[c(1)], list(dimnames(ps$EX_tau_eta)[[c(3)]][2:(1 + num_inter)]))
  }

  dimnames(is_EX_comp) <- dimnames(ps$EX_prob_comp)
  dimnames(draw_beta) <- c(dimnames(ps$EX_prob_comp), list(coefficient = c("intercept", "log_slope")))

  dimnames(is_EX_inter) <- dimnames(ps$EX_prob_inter)
  dimnames(draw_eta) <- c(dimnames(ps$EX_prob_inter))

  list(
    draw_beta = draw_beta,
    draw_eta = draw_eta,
    EX_mu_comp = EX_mu_comp,
    EX_mu_inter = EX_mu_inter,
    EX_tau_comp = EX_tau_comp,
    EX_tau_inter = EX_tau_inter,
    EX_corr_comp = EX_corr_comp,
    EX_corr_inter = EX_corr_inter,
    log_beta = log_beta,
    eta = eta,
    is_EX_comp = is_EX_comp,
    is_EX_inter = is_EX_inter
  )
}

#'
#' Simulates a draw from the prior and fake data for it. This will be
#' the data generating step in the simulation. The function recieves
#' the problem data, job specifics and a blrmfit object which defines
#' the prior to sample and the design matrix.
#'
simulate_fake <- function(scenario) {
  prior_draw <- sample_prior(scenario)

  beta_group_rv <- as_rvar(prior_draw$draw_beta)
  eta_group_rv <- as_rvar(prior_draw$draw_eta)

  standata <- scenario$base_fit$standata

  ## logit by data-row
  ## draw_mu <- with(standata, blrm_logit_grouped_vec(group, stratum, X_comp, X_inter, prior_draw$draw_beta, prior_draw$draw_eta))

  ## TODO: can we avoid to use the internal function here?
  draw_mu <- with(standata, OncoBayes2:::blrm_logit_grouped_rv(group, X_comp, X_inter, beta_group_rv, eta_group_rv))

  num_trials <- standata$r + standata$nr

  yrep <- rbinom(length(num_trials), num_trials, inv_logit(draw_mu))

  list(yrep = yrep, draw = prior_draw)
}

restore_draw_dims <- function(standata, draw) {
  num_comp <- standata$num_comp
  num_inter <- standata$num_inter
  num_strata <- standata$num_strata
  num_groups <- standata$num_groups

  draw$mu_log_beta <- array(draw$mu_log_beta, c(num_comp, 2))
  draw$tau_log_beta_raw <- array(draw$tau_log_beta_raw, c(num_strata, num_comp, 2))
  draw$L_corr_log_beta <- array(draw$L_corr_log_beta, c(num_comp, 2, 2))
  draw$log_beta_raw <- array(draw$log_beta_raw, c(2 * num_groups, num_comp, 2))

  if (num_inter != 0) {
    draw$eta_raw <- array(draw$eta_raw, c(2 * num_groups, num_inter))
    draw$mu_eta <- array(draw$mu_eta, c(num_inter))
    draw$tau_eta_raw <- array(draw$tau_eta_raw, c(num_strata, num_inter))
    draw$L_corr_eta <- matrix(draw$L_corr_eta, num_inter, num_inter)
  } else {
    draw$eta_raw <- array(0, c(2 * num_groups, num_inter))
    draw$mu_eta <- array(0, c(num_inter))
    draw$tau_eta_raw <- array(0, c(num_strata, num_inter))
    draw$L_corr_eta <- matrix(1, num_inter, num_inter)
  }

  draw
}

#' extracts from a given fit the mass matrix, stepsize and a draw from
#' the typical set. The warmup info from multiple chains is being
#' averaged together to obtain less noisy estimates.
learn_warmup_info <- function(standata, stanfit) {
  gmean <- function(x) exp(mean(log(x)))
  have_inter <- standata$num_inter > 0
  sampled_params <- c("log_beta_raw", "mu_log_beta", "tau_log_beta_raw", "L_corr_log_beta")
  if (have_inter) {
    sampled_params <- c(sampled_params, "eta_raw", "mu_eta", "tau_eta_raw", "L_corr_eta")
  }
  posterior_draws <- merge_chains(as_draws_rvars(as.array(stanfit, pars = sampled_params)))
  s <- floor(seq.int(1, ndraws(posterior_draws), length = 10))
  draws <- list()
  for (i in s) {
    init <- subset_draws(posterior_draws, iteration = i)
    init <- lapply(lapply(init, draws_of), adrop, drop = 1, one.d.array = TRUE)
    draws <- c(draws, list(init))
  }
  warmup_info <- extract_adaptation_info_stanfit(stanfit)
  warmup_info$stepsize <- gmean(warmup_info$stepsize)
  warmup_info$inv_metric <- apply(warmup_info$inv_metric, 1, gmean)
  c(warmup_info, list(draws = lapply(draws, restore_draw_dims, standata = standata)))
}

#'
#'
#' Procedure to fit each fake data set using our fitting
#' procedure. This method obtains the problem data, job details and an
#' **instance** of the scenario as generated by `simulate_fake`.
#'

fit_exnex <- function(yrep, draw, scenario, ..., save_fit = FALSE) {
  ## yrep <- instance$yrep
  ## draw <- instance$draw
  group_draws <- list()
  group_draws$draw_beta <- draw$draw_beta
  group_draws$draw_eta <- draw$draw_eta

  ## pars <- job$pars$prob.pars

  dref <- scenario$dref
  sim_data <- scenario$base_fit$data
  sim_data$num_toxicities <- yrep

  blrm_args <- scenario$blrm_args

  have_warmup_info <- c("warmup_info") %in% names(scenario)

  if (have_warmup_info) {
    ## use a randomly selected warmup info from the ones provided
    fit_warmup_info <- sample(scenario$warmup_info, 1)[[1]]
    blrm_args <- scenario$blrm_args_with_warmup_info
    blrm_args$init <- sample(fit_warmup_info$draws, blrm_args$chains)
    blrm_args$control <- modifyList(
      blrm_args$control,
      list( ## adapt_inv_metric=fit_warmup_info$inv_metric,
        stepsize = 3 * fit_warmup_info$stepsize
      )
    )
  }

  ## NOTE: Previously the pattern
  ## "fit <- update(scenario$base_fit," had been used... BUT this
  ## causes major issues when running things in parallel (somehow R
  ## tried to update a vanilla OncoBayes2 object, which then was
  ## based on the installed OncoBayes2 package such that on remote
  ## nodes outdated code got used)
  fit <- do.call(
    blrm_exnex,
    modifyList(
      blrm_args,
      list(
        data = sim_data,
        verbose = FALSE,
        save_warmup = !have_warmup_info
      )
    )
  )

  np <- nuts_params(fit, inc_warmup = FALSE)

  n_divergent <- sum(subset(np, Parameter == "divergent__")$Value)
  accept_stat <- mean(subset(np, Parameter == "accept_stat__")$Value)

  params_comp <- c("mu_log_beta", "tau_log_beta", "beta_group")
  params_inter <- c()
  if (fit$has_inter) {
    params_inter <- c("mu_eta", "tau_eta", "eta_group")
  }
  params <- c(params_comp, params_inter)

  samp_diags_sum <- summarise_draws(as.array(fit$stanfit, pars = params), "rhat", "ess_bulk", "ess_tail") %>%
    summarize(max_rhat = max(rhat), min_ess_bulk = min(ess_bulk), min_ess_tail = min(ess_tail))

  samp_diags_lp_sum <- summarise_draws(as.array(fit$stanfit, pars = "lp__"), "rhat", "ess_bulk", "ess_tail")

  assert_that(nsamples(fit) > 1023)
  suppressMessages(post_thin <- subset_draws(as_draws_rvars(as.array(fit$stanfit, pars = params)), draw = seq(1, nsamples(fit), length = 1024 - 1)))

  dim(post_thin$tau_eta)
  lapply(post_thin, dim)

  if (fit$has_inter) {
    ## BUG in posterior??!!
    ## the last dimension gets dropped for tau_eta
    dim(post_thin$tau_eta) <- c(fit$standata$num_strata, fit$standata$num_inter)
  }

  calc_rank <- function(sample_rv, draw) {
    sample <- draws_of(sample_rv)
    sdims <- dim(sample)
    assert_that(all(sdims[-1] == dim(draw)))
    draw_margins <- 2:length(sdims)
    res <- array(apply(sweep(sample, draw_margins, draw) < 0, draw_margins, sum), dim = sdims[-1])
    dimnames(res) <- dimnames(draw)
    res
  }

  rank1 <- mapply(calc_rank,
    post_thin[params_comp],
    c(draw[c("EX_mu_comp", "EX_tau_comp")], list(beta_group = group_draws$draw_beta)),
    SIMPLIFY = FALSE
  )

  if (fit$has_inter) {
    rank1 <- c(
      rank1,
      mapply(calc_rank,
        post_thin[params_inter],
        c(draw[c("EX_mu_inter", "EX_tau_inter")], list(eta_group = group_draws$draw_eta)),
        SIMPLIFY = FALSE
      )
    )
  }

  flatten_array <- function(data, var) {
    if (missing(var)) {
      var <- deparse(substitute(data))
    }
    idx <- expand.grid(lapply(dim(data), seq))
    num_dim <- length(dim(data))
    size <- nrow(idx)
    values <- sapply(seq_len(size), function(r) {
      asub(data, idx[r, ], drop = FALSE)
    })
    idx[seq_len(num_dim)] <- lapply(seq_len(num_dim), function(col) dimnames(data)[[col]][idx[, col]])
    names(values) <- paste0(var, "[", do.call("paste", c(idx[seq_len(num_dim)], list(sep = ","))), "]")
    res <- matrix(values, nrow = 1, ncol = size)
    colnames(res) <- names(values)
    as.data.frame(res)
  }

  rank_wide <- bind_cols(mapply(flatten_array, rank1, names(rank1), SIMPLIFY = FALSE))

  res <- list(
    rank = rank_wide,
    min_Neff_bulk = pull(samp_diags_sum, "min_ess_bulk"),
    min_Neff_tail = pull(samp_diags_sum, "min_ess_tail"),
    n_divergent = n_divergent,
    max_Rhat = pull(samp_diags_sum, "max_rhat"),
    lp_ess_bulk = pull(samp_diags_lp_sum, "ess_bulk"),
    lp_ess_tail = pull(samp_diags_lp_sum, "ess_tail")
  )

  if (save_fit) {
    res$fit <- fit
  }

  if (!have_warmup_info) {
    res <- c(res, learn_warmup_info(fit$standata, fit$stanfit))
  }

  res$stepsize <- exp(mean(log(subset(np, Parameter == "stepsize__" & Iteration == 1)$Value)))
  res$accept_stat <- accept_stat

  return(res)
}


run_sbc_case <- function(job.id, repl, data_scenario, rng_seed, example_models_cmq) {
  source("sbc_job_startup.R")
  RNGkind("L'Ecuyer-CMRG")
  print(rng_seed)
  .Random.seed <<- rng_seed

  # if(is.null(example_models_cmq)) {
  #   example_models_cmq <- example_models
  # }

  runtime <- system.time({
    load_OB2_dev(ob2_lib_dir)
    scenario <- example_models_cmq[[data_scenario]]
    fake <- simulate_fake(scenario)
    fit <- fit_exnex(fake$yrep, fake$draw, scenario)
  })

  c(list(job.id = job.id, time.running = unname(runtime["elapsed"])), fit)
}



# AB: not currently using this function from rbest SBC...
# scale_ranks <- function(Nbins, scale=1) {
#   ## scale must evenly divide the total number of bins
#   assert_that(round(Nbins/scale) == Nbins/scale)
#   breaks <- (0:(Nbins/scale))
#   Nbreaks <- length(breaks)
#   function(scen) {
#     vars <- grep("^rank.", names(scen), value=TRUE)
#     res <- lapply(vars, function(v) hist(ceiling((scen[[v]]+1)/scale), breaks=breaks, plot=FALSE, include.lowest=FALSE)$counts)
#     names(res) <- gsub("^rank", "count", vars)
#     res$rank <- breaks[-Nbreaks]
#     res <- as.data.frame(do.call(cbind, res))
#     res
#   }
# }


#' extract for a stanfit object from rstan the adaptation information
#' @param fit cmdstanr fit
#' @keywords internal
extract_adaptation_info_stanfit <- function(fit) {
  info <- sapply(rstan::get_adaptation_info(fit), strsplit, "\n")
  ex_stepsize <- function(chain_info) {
    stepsize_line <- which(grepl("Step size", chain_info))
    as.numeric(strsplit(chain_info[stepsize_line], " = ")[[1]][2])
  }
  ex_mass <- function(chain_info) {
    metric_line <- which(grepl("inverse mass matrix", chain_info)) + 1
    as.numeric(strsplit(sub("^#", "", chain_info[metric_line]), ", ")[[1]])
  }
  stepsize <- sapply(info, ex_stepsize)
  inv_metric <- do.call(cbind, lapply(info, ex_mass))
  colnames(inv_metric) <- names(stepsize) <- paste0("chain_", seq_along(info))
  list(stepsize = stepsize, inv_metric = inv_metric)
}
