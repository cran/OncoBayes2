## ----SETTINGS-knitr, include=FALSE--------------------------------------------
## knitr settings used to build vignettes
library(OncoBayes2)
library(posterior)
library(RBesT)
library(dplyr)
library(tidyr)
library(knitr)
library(ggplot2)
ggplot2::theme_set(bayesplot::bayesplot_theme_get())
knitr::knit_hooks$set(pngquant = knitr::hook_pngquant)
knitr::opts_chunk$set(
  dev = "ragg_png",
  dpi = 72,
  fig.retina = 1.5,
  fig.width = 1.62*4,
  fig.height = 4,
  fig.align = "center",
  out.width = "100%",
  pngquant = "--speed=1 --quality=50"
  )

## ----SETTINGS-sampling, include=FALSE-----------------------------------------
## sampling settings used to build vignettes
## setup up fast sampling when run on CRAN
not_CRAN <- Sys.getenv("NOT_CRAN", "false") == "true" 
## NOTE: for running this vignette locally, please uncomment the
## following line:
## not_CRAN <- TRUE
.user_mc_options <- list()

if (!not_CRAN) {
    .user_mc_options <- options(OncoBayes2.MC.warmup=40, OncoBayes2.MC.iter=100, OncoBayes2.MC.chains=1, OncoBayes2.MC.save_warmup=FALSE, OncoBayes2.MC.control = list(adapt_delta=0.85), mc.cores=1)
} else {
    .user_mc_options <- options(OncoBayes2.MC.warmup=500, OncoBayes2.MC.iter=1000, OncoBayes2.MC.chains=4, OncoBayes2.MC.save_warmup=FALSE, OncoBayes2.MC.control = list(adapt_delta=0.99), mc.cores=1)
}
set.seed(6475863)

## -----------------------------------------------------------------------------
library(OncoBayes2)
library(RBesT)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)

## ----mtd_hist, echo = FALSE---------------------------------------------------
knitr::kable(
  data.frame(
    Drug = c("Drug A", "Drug B"),
    "Declared MTD" = c("100mg", "3mg"),
    check.names = FALSE
  )
)

## ----hist_A, echo = FALSE-----------------------------------------------------
knitr::kable(
  data.frame(
    "Drug A" = c(12.5, 25, 50, 80, 100, 150),
    "Evaluable patients" = c(1, 1, 3, 9, 23, 3),
    "DLTs" = c(0, 0, 0, 1, 4, 2),
    check.names = FALSE
  ),
  caption = "Historical data from single-agent study of Drug A"
)

## ----hist_B, echo = FALSE-----------------------------------------------------
knitr::kable(
  data.frame(
    "Drug B dose" = c(0.125, 0.25, 0.5, 1, 2, 2.5, 3, 4),
    "Evaluable patients" = c(2, 1, 2, 2, 3, 7, 12, 3),
    "DLTs" = c(0, 0, 0, 0, 1, 0, 0, 1),
    check.names = FALSE
  ),
  caption = "Historical data from single-agent study of Drug B"
)

## ----heterogeneity, echo = FALSE----------------------------------------------
knitr::kable(
  data.frame(
    "Heterogeneity degree" = c("small", "moderate", "substantial", "large", "very large"),
    "$\\tau_\\alpha$" = c(0.125, 0.25, 0.5, 1.0, 2.0),
    "$\\tau_\\beta$" = c(0.0625, 0.125, 0.25, 0.5, 1.0),
    check.names = FALSE
  ),
  caption = "Heterogeneity categorization for intercept $\\alpha$ and slope $\\beta$"
)

## ----map_model_A--------------------------------------------------------------
# historical data for drug A
hist_A <- data.frame(
  group_id = factor("trial_A"),
  drug_A = c(12.5, 25, 50, 80, 100, 150),
  num_patients = c(1, 1, 3, 9, 23, 3),
  num_toxicities = c(0, 0, 0, 1, 4, 2)
)

# reference dose for drug A
dref_A <- 80

# bivariate normal prior for (intercept, log-slope) 
prior_A <- mixmvnorm(c(1.0,
                       logit(0.2), 0,
                       diag(c(1, (log(4) / 1.96)^2 ))))

# between-trial heterogeneity is centered on "substantial"
tau_substantial_prior <- mixmvnorm(c(1.0,
                                     log(c(0.5, 0.25)),
                                     diag(rep(log(2) / 1.96, 2)^2)))

# fit the model to the historical data
map_model_A <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~ 
    1 + log(drug_A / dref_A) | 
    0 | 
    group_id,
  data = hist_A,
  prior_EX_mu_comp = list(prior_A),
  prior_EX_tau_comp = tau_substantial_prior,
  prior_EX_prob_comp = matrix(1, nrow = 1, ncol = 1),
  prior_tau_dist = 1,
  sample_map = TRUE
)

## ----map_model_B--------------------------------------------------------------
# historical data for drug B
hist_B <- data.frame(
  group_id = factor("trial_B"),
  drug_B = c(0.125, 0.25, 0.5, 1, 2, 2.5, 3, 4),
  num_patients = c(2, 1, 2, 2, 3, 7, 12, 3),
  num_toxicities = c(0, 0, 0, 0, 1, 0, 0, 1)
)

# reference dose
dref_B <- 1

# prior for (intercept, log-slope)
prior_B <- prior_A

# between-trial heterogeneity is centered on "small", but with a
# greater degree of uncertainty on this classification
tau_small_prior <- mixmvnorm(c(1.0,
                               log(c(0.125, 0.0625)),
                               diag(rep(log(4) / 1.96, 2)^2)))

# fit model to historical data
map_model_B <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~ 
    1 + log(drug_B / dref_B) | 
    0 | 
    group_id,
  data = hist_B,
  prior_EX_mu_comp = list(prior_B),
  prior_EX_tau_comp = tau_small_prior,
  prior_EX_prob_comp = matrix(1, nrow = 1, ncol = 1),
  prior_tau_dist = 1,
  sample_map = TRUE
)

## ----beta_group---------------------------------------------------------------
rvar_map_log_beta_A <- as_draws_rvars(map_model_A, variable="map_log_beta")$map_log_beta
rvar_map_log_beta_B <- as_draws_rvars(map_model_B, variable="map_log_beta")$map_log_beta

## ----beta_dimnames------------------------------------------------------------
rvar_map_log_beta_A[1,"log(drug_A/dref_A)",,drop=TRUE]
rvar_map_log_beta_B[1,"log(drug_B/dref_B)",,drop=TRUE]

## ----mixfit-------------------------------------------------------------------
map_mix_A <- automixfit(draws_of(rvar_map_log_beta_A[1,"log(drug_A/dref_A)",,drop=TRUE]),
                        type = "mvnorm")
map_mix_A

map_mix_B <- automixfit(draws_of(rvar_map_log_beta_B[1,"log(drug_B/dref_B)",,drop=TRUE]),
                        type = "mvnorm")
map_mix_B

## ----mixfit_plot_A------------------------------------------------------------
plot(map_mix_A)$mixpairs
# plot(map_mix_B)$mixpairs  # respective plot for drug B

## ----map_combo----------------------------------------------------------------
# one row of dummy data with number of patients set to zero in order
# to allow the model to be fit
dummy_data_AB <- data.frame(
  group_id = factor("trial_AB", levels = c("trial_AB")),
  drug_A = 1,
  drug_B = 1,
  num_patients = 0,
  num_toxicities = 0
)

interaction_prior <- mixmvnorm(c(1, 0, (log(9)/1.96)^2))

# fit a double-combination model using the MAP priors
map_combo <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
    1 + I(log(drug_A / dref_A)) |
    1 + I(log(drug_B / dref_B)) |
    0 + I(2 * (drug_A/dref_A * drug_B/dref_B) / (1 + drug_A/dref_A * drug_B/dref_B)) |
    group_id,
  data = dummy_data_AB,
  prior_EX_mu_comp = list(map_mix_A, map_mix_B), # BVN mixtures for (log-alpha, log-beta) for drugs A and B
  prior_EX_mu_inter = interaction_prior,
  
  # shut off hierarchical part of model --
  prior_tau_dist = NULL
)

## ----plot_map-----------------------------------------------------------------
doses <- expand_grid(
  group_id = "trial_AB",
  drug_A = c(25, 50, 80, 100),
  drug_B = c(0.5, 1, 3)
) %>%
  arrange(drug_B, drug_A)

plot_toxicity_curve(map_combo,
                    newdata = doses,
                    x = vars(drug_A),
                    group = ~ drug_B)  +
  theme(legend.position="bottom")

## ----summ_map-----------------------------------------------------------------
map_summ <- summary(map_combo, newdata = doses, interval_prob = c(0, 0.16, 0.33, 1))
kable(
  cbind(doses[c("drug_A", "drug_B")], map_summ[c("mean", "sd", "(0.33,1]")]),
  col.names = c("Drug A", "Drug B", "mean", "sd", "P(DLT rate > 0.33)"),
  digits = 2,
  caption = "Posterior summary statistics for P(DLT) by dose"
)

## ----pred_summ_map------------------------------------------------------------
map_pred_summ <- summary(map_combo,
                         newdata = mutate(doses, num_toxicities=0, num_patients=6),
                         predictive=TRUE, interval_prob = c(-1, 0, 1, 6))
kable(
  cbind(doses[c("drug_A", "drug_B")], map_pred_summ[c("(-1,0]", "(0,1]", "(1,6]")]),
  col.names = c("Drug A", "Drug B", "Pr(0 of 6 DLT)", "Pr(1 of 6 DLT)", "Pr(>=2 of 6 DLT)"),
  digits = 2,
  caption = "Posterior predictive summary for Pr(DLT) by dose"
)

## ----rmap---------------------------------------------------------------------
# weakly-informative mixture component is chosen here to be equal to
# the priors used in the first place
weak_A <- prior_A
weak_B <- prior_B

# adding this to the drug-A and drug-B MAP priors with weight 0.1
mix_A_robust <- mixcombine(map_mix_A, weak_A, weight = c(0.9, 0.1))
mix_B_robust <- mixcombine(map_mix_B, weak_B, weight = c(0.9, 0.1))

# robust MAP priors
mix_robust <- list(mix_A_robust, mix_B_robust)

# fit a double-combination model using the rMAP priors
robust_map_combo <- update(map_combo, prior_EX_mu_comp = mix_robust)


## ----rmap_plot, eval=FALSE----------------------------------------------------
#  # disabled plot here
#  plot_toxicity_curve(robust_map_combo,
#                      newdata = doses,
#                      x = vars(drug_A),
#                      group = ~ drug_B)  +
#    theme(legend.position="bottom")

## ----mac_groups---------------------------------------------------------------
groups <- c("trial_A", "trial_B", "trial_AB")
mac_data <- bind_rows_0(hist_A, hist_B)
mac_data$group_id <- factor(as.character(mac_data$group_id), levels = groups)

## ----mac_combo----------------------------------------------------------------
dummy_tau_inter <- mixmvnorm(c(1.0, log(1E-4), (log(2) / 1.96)^2))

num_comp   <- 2 # number of treatment components
num_inter  <- 1
num_groups <- length(groups)

mac_combo <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
    1 + I(log(drug_A / dref_A)) |
    1 + I(log(drug_B / dref_B)) |
    0 + I(2 * (drug_A/dref_A * drug_B/dref_B) / (1 + drug_A/dref_A * drug_B/dref_B)) |
    group_id,
  data = mac_data,
  prior_EX_mu_comp = list(prior_A, prior_B),
  prior_EX_tau_comp = list(tau_substantial_prior, tau_small_prior),
  prior_EX_mu_inter = interaction_prior,
  prior_EX_tau_inter = dummy_tau_inter,
  prior_is_EXNEX_comp = rep(FALSE, num_comp),
  prior_is_EXNEX_inter = rep(FALSE, num_inter),
  prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
  prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
  prior_tau_dist = 1
)

## ----map_mac_compare----------------------------------------------------------
mac_summ <- summary(mac_combo, newdata = doses, interval_prob = c(0, 0.16, 0.33, 1))
kable(
  cbind(doses[c("drug_A", "drug_B")],
        map_summ$mean, mac_summ$mean,
        abs(map_summ$mean - mac_summ$mean),
        map_summ[["(0.33,1]"]], mac_summ[["(0.33,1]"]]),
  col.names = c("Drug A", "Drug B",
                "MAP mean", "MAC mean",
                "|MAP mean - MAC mean|",
                "MAP P(DLT rate > 0.33)", "MAC P(DLT rate > 0.33)"),
  digits = 2,
  caption = "Posterior summary statistics for P(DLT) by dose"
)

## -----------------------------------------------------------------------------
map_mix_A

## -----------------------------------------------------------------------------
map_mix_B

## -----------------------------------------------------------------------------
sessionInfo()

## ----include=FALSE------------------------------------------------------------
## restore previous global user options
options(.user_mc_options)

