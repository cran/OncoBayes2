#' Summarise model prior
#'
#' @description
#' Extracts a summary of the prior in a structured data format.
#'
#' @param object \code{blrmfit} (\code{blrm_trial}) object as returned from \code{\link{blrm_exnex}} (\code{\link{blrm_trial}}) analysis
#' @param digits number of digits to show
#' @param ... ignored by the function
#'
#' @details The summary of the prior creates a structured
#'     representation of the specified prior from a
#'     \code{\link{blrm_exnex}} (\code{\link{blrm_trial}}) analysis.
#'
#' @return Returns an analysis specific list, which has it's own
#'     \code{print} function. The returned list contains arrays which
#'     represent the prior in a structured format.
#'
#' @template start-example
#' @examples
#' ## run combo2 analysis which defines blrmfit model object
#' example_model("combo2", silent = TRUE)
#'
#' prior_summary(blrmfit)
#'
#' prior_sum <- prior_summary(blrmfit)
#' names(prior_sum)
#'
#' ## the entries of the prior list are labelled arrays
#' dimnames(prior_sum$EX_mu_log_beta)
#'
#' @template stop-example
#'
#' @method prior_summary blrmfit
#' @aliases prior_summary
#' @export
prior_summary.blrmfit <- function(object, digits = 2, ...) {
  standata <- object$standata
  labels <- object$labels

  x <- list()
  x$has_inter <- object$has_inter

  x$is_EXNEX_comp <- .label_array(standata$prior_is_EXNEX_comp, component = labels$component)

  x$EX_prob_comp <- .label_array(standata$prior_EX_prob_comp, group = object$group_fct, component = labels$component)
  is_EX_comp_idx <- which(x$is_EXNEX_comp == 0)
  x$EX_prob_comp[, is_EX_comp_idx] <- 1

  x$EX_mu_log_beta <- with(standata, .parse_mu_log_beta_mix(prior_EX_mu_comp_Nc, prior_EX_mu_comp_w, prior_EX_mu_comp_m, prior_EX_mu_comp_sigma, labels))

  x$NEX_mu_log_beta <- with(standata, .parse_mu_log_beta_mix(prior_NEX_mu_comp_Nc, prior_NEX_mu_comp_w, prior_NEX_mu_comp_m, prior_NEX_mu_comp_sigma, labels))

  x$EX_tau_log_beta <- with(standata, .parse_tau_beta_mix(prior_EX_tau_comp_Nc, prior_EX_tau_comp_w, prior_EX_tau_comp_m, prior_EX_tau_comp_sigma, labels))
  x$EX_corr_eta_comp <- .label_array(standata$prior_EX_corr_eta_comp, component = labels$component)

  x$tau_dist <- standata$prior_tau_dist

  if (x$has_inter) {
    x$is_EXNEX_inter <- .label_array(standata$prior_is_EXNEX_inter, interaction = labels$param_eta)
    x$EX_prob_inter <- .label_array(standata$prior_EX_prob_inter, group = object$group_fct, interaction = labels$param_eta)
    is_EX_inter_idx <- which(x$is_EXNEX_inter == 0)
    x$EX_prob_inter[, is_EX_inter_idx] <- 1

    x$EX_mu_eta <- with(standata, .parse_mu_eta_mix(prior_EX_mu_inter_w, prior_EX_mu_inter_m, prior_EX_mu_inter_sigma, labels))
    x$NEX_mu_eta <- with(standata, .parse_mu_eta_mix(prior_NEX_mu_inter_w, prior_NEX_mu_inter_m, prior_NEX_mu_inter_sigma, labels))

    x$EX_tau_eta <- with(standata, .parse_tau_eta_mix(prior_EX_tau_inter_w, prior_EX_tau_inter_m, prior_EX_tau_inter_sigma, labels))

    x$EX_corr_eta_inter <- array(standata$prior_EX_corr_eta_inter)
    dimnames(x$EX_corr_eta_inter) <- list("interaction")
  }

  x$num_strata <- standata$num_strata
  x$num_groups <- standata$num_groups

  structure(x,
    class = "prior_summary.blrmfit",
    model_name = deparse(substitute(object)),
    print_digits = digits
  )
}

#' @export
#' @method print prior_summary.blrmfit
print.prior_summary.blrmfit <- function(x, digits, ...) {
  cat("Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability\n\n")

  if (missing(digits)) {
    digits <- attr(x, "print_digits")
  }

  tau_str <- if (x$tau_dist == 0) {
    "known"
  } else if (x$tau_dist == 1) {
    "log-normal"
  } else if (x$tau_dist == 2) {
    "half-normal"
  } else {
    stop("Unkown tau prior distribution.")
  }

  cat("Mixture configuration\n")
  cat("---------------------\n")
  cat("EXNEX components :", sum(x$is_EXNEX_comp), "\n")
  print(x$is_EXNEX_comp)
  cat("\n")
  if (x$has_inter) {
    cat("EXNEX interactions:", sum(x$is_EXNEX_inter), "\n")
    print(x$is_EXNEX_inter)
    cat("\n")
  }
  cat("Prior probability for exchangeability per group\n")
  print(x$EX_prob_comp)
  cat("\n")
  if (x$has_inter) {
    print(x$EX_prob_inter)
    cat("\n")
  }
  cat("EXchangable hyperparameter priors\n")
  cat("---------------------------------\n")
  cat("Component parameters\n")
  cat("Mean mu_log_beta\n")
  print(ftable(x$EX_mu_log_beta, row.vars = c(2, 3)), digits = digits)
  cat("\n")
  cat(paste0("Heterogeneity tau_log_beta (", tau_str, ")\n"))
  print(ftable(x$EX_tau_log_beta, row.vars = c(1, 3, 4)), digits = digits)
  cat("\nCorrelation LKJ\n")
  print(x$EX_corr_eta_comp, digits = digits)

  if (x$has_inter) {
    cat("\nInteraction parameters\n")
    cat("Mean mu_eta\n")
    print(ftable(x$EX_mu_eta, row.vars = 1), digits = digits)
    cat("\n")
    cat(paste0("Heterogeneity tau_eta (", tau_str, ")\n"))
    print(ftable(x$EX_tau_eta, row.vars = c(1, 2)), digits = digits)
    cat("\nCorrelation LKJ\n")
    print(x$EX_corr_eta_inter, digits = digits)
  } else {
    cat("\nModel has no interaction parameters.\n")
  }

  cat("\n")
  cat("NonEXchangable priors\n")
  cat("---------------------\n")
  cat("Component parameters\n")
  cat("Mean mu_log_beta\n")
  print(ftable(x$NEX_mu_log_beta, row.vars = c(2, 3)), digits = digits)

  if (x$has_inter) {
    cat("\nInteraction parameters\n")
    cat("Mean mu_eta\n")
    print(ftable(x$NEX_mu_eta, row.vars = 1), digits = digits)
  } else {
    cat("\nModel has no interaction parameters.\n")
  }

  invisible(x)
}

#' @method prior_summary blrm_trial
#' @export
prior_summary.blrm_trial <- function(object, ...) {
  .assert_is_blrm_trial_and_prior_is_set(object)

  x <- list()
  x$prior_summary.blrmfit <- prior_summary(object$blrmfit, ...)

  structure(x, class = "prior_summary.blrm_trial")
}

#' @method print prior_summary.blrm_trial
#' @export
print.prior_summary.blrm_trial <- function(x, ...) {
  print(x$prior_summary.blrmfit, ...)
}

## internal -----

.label_array <- function(data, ...) {
  labs <- list(...)
  assert_that(length(labs) == length(dim(data)), msg = "Number of labels must match dimension of input data.")
  assert_that(all(dim(data) == sapply(labs, nlevels)), msg = "Number of factor levels must match array dimensionality.")
  dimnames(data) <- lapply(labs, levels)
  data
}

.parse_mu_log_beta_mix <- function(mu_Nc_comp, mu_w_comp, mu_mean_comp, mu_sigma_comp, labels) {
  num_comp <- length(mu_Nc_comp)
  max_Nc <- max(mu_Nc_comp)
  mix_comp <- factor(paste0("comp_", 1:max_Nc))
  mu_sd_comp <- array(0.0, dim = dim(mu_mean_comp))
  mu_rho_comp <- array(0.0, dim = dim(mu_mean_comp)[c(1, 2)])
  for (i in 1:num_comp) {
    for (j in 1:max_Nc) {
      mu_sd_comp[i, j, 1] <- sqrt(mu_sigma_comp[i, j, 1, 1])
      mu_sd_comp[i, j, 2] <- sqrt(mu_sigma_comp[i, j, 2, 2])
      mu_rho_comp[i, j] <- mu_sigma_comp[i, j, 1, 2] / (mu_sd_comp[i, j, 1] * mu_sd_comp[i, j, 2])
    }
  }
  mu_log_beta_mean <- .label_array(mu_mean_comp, component = labels$component, mix = mix_comp, coefficient = labels$param_log_beta)
  mu_log_beta_sd <- .label_array(mu_sd_comp, component = labels$component, mix = mix_comp, coefficient = labels$param_log_beta)
  mu_log_beta_w <- .label_array(mu_w_comp, component = labels$component, mix = mix_comp)
  mu_log_beta_rho <- .label_array(mu_rho_comp, component = labels$component, mix = mix_comp)
  mu_log_beta_scalar <- abind(weight = mu_log_beta_w, correlation = mu_log_beta_rho, along = 0)
  mu_log_beta <- abind(
    w = mu_log_beta_scalar[1, , , drop = FALSE],
    m = aperm(mu_log_beta_mean, c(3, 1, 2)),
    s = aperm(mu_log_beta_sd, c(3, 1, 2)),
    rho = mu_log_beta_scalar[2, , , drop = FALSE],
    along = 1
  )
  names(dimnames(mu_log_beta)) <- c("prior", "component", "mix")
  dimnames(mu_log_beta)$prior <- c("weight", paste0("m_", labels$param_log_beta), paste0("s_", labels$param_log_beta), "rho")
  return(mu_log_beta)
}

.parse_tau_beta_mix <- function(tau_comp_Nc, tau_comp_w, tau_comp_m, tau_comp_sigma, labels) {
  num_strata <- dim(tau_comp_Nc)[1]
  tau_beta <- list()
  for (s in 1:num_strata) {
    tau_beta[[s]] <- .parse_mu_log_beta_mix(
      adrop(tau_comp_Nc[s, , drop = FALSE], 1),
      adrop(tau_comp_w[s, , , drop = FALSE], 1),
      adrop(tau_comp_m[s, , , , drop = FALSE], 1),
      adrop(tau_comp_sigma[s, , , , , drop = FALSE], 1),
      labels
    )
  }
  tau_beta <- abind(tau_beta, along = 0)
  dimnames(tau_beta)[[1]] <- paste0("stratum_", 1:num_strata)
  dimnames(tau_beta)[[2]] <- c("weight", "m_tau_intercept", "m_tau_log_slope", "s_tau_intercept", "s_tau_log_slope", "rho")
  names(dimnames(tau_beta)) <- c("stratum", "prior", "component", "mix")
  tau_beta
}

.parse_mu_eta_mix <- function(mu_w_inter, mu_mean_inter, mu_sigma_inter, labels) {
  Nc <- length(mu_w_inter)
  comps <- list()
  for (i in 1:Nc) {
    comps[[i]] <- c(mu_w_inter[i], mu_mean_inter[i, ], mu_sigma_inter[i, , ])
  }
  mv <- do.call(RBesT::mixmvnorm, comps)
  mvm <- matrix(mv, dim(mv))
  dimnames(mvm) <- dimnames(mv)
  names(dimnames(mvm)) <- c("prior", "mix")
  t(mvm)
}

.parse_tau_eta_mix <- function(tau_inter_w, tau_inter_m, tau_inter_sigma, labels) {
  num_strata <- dim(tau_inter_w)[1]
  tau_eta <- list()
  for (s in 1:num_strata) {
    tau_eta[[s]] <- .parse_mu_eta_mix(
      adrop(tau_inter_w[s, , drop = FALSE], 1),
      adrop(tau_inter_m[s, , , drop = FALSE], 1),
      adrop(tau_inter_sigma[s, , , , drop = FALSE], 1),
      labels
    )
  }
  tau_eta <- abind(tau_eta, along = 0)
  dimnames(tau_eta)[[1]] <- paste0("stratum_", 1:num_strata)
  names(dimnames(tau_eta)) <- c("stratum", "mix", "prior")
  tau_eta
}
