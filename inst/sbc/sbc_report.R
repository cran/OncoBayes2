#' ---
#' title: Simulation based calibration for OncoBayes2
#' author: ""
#' date: "`r date()`"
#' output: html_vignette
#' params:
#'   include_plots: FALSE
#' ---
#'
#+ include=FALSE
library(knitr)
library(tools)
library(assertthat)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
theme_set(theme_bw())
source("sbc_tools.R")

knitr::opts_chunk$set(
  fig.width = 1.62*4,
  fig.height = 4,
  cache=FALSE,
  echo=FALSE
)
#'
#' This report documents the results of a simulation based calibration
#' (SBC) run for `OncoBayes2`. TODO
#'
#' The calibration data presented here has been generated at and with
#' the `OncoBayes` git version as:
cat(readLines("calibration.md5"), sep="\n")
#'
#' The MD5 hash of the calibration data file presented here must match
#' the above listed MD5:
md5sum("calibration.rds")
#'
#' # Introduction
#'
#' Simulation based calibration (SBC) is a necessary condition which
#' must be met for any Bayesian analysis with proper priors. The
#' details are presented in Talts, et. al (see
#' https://arxiv.org/abs/1804.06788).
#'
#' Self-consistency of any Bayesian analysis with a proper prior:
#'
#' $$ p(\theta) = \iint \mbox{d}\tilde{y} \, \mbox{d}\tilde{\theta} \, p(\theta|\tilde{y}) \, p(\tilde{y}|\tilde{\theta}) \, p(\tilde{\theta}) $$
#' $$ \Leftrightarrow p(\theta) = \iint \mbox{d}\tilde{y} \, \mbox{d}\tilde{\theta} \, p(\theta,\tilde{y},\tilde{\theta}) $$
#'
#' SBC procedure:
#'
#' Repeat $s=1, ..., S$ times:
#'
#' 1. Sample from the prior $$\tilde{\theta} \sim p(\theta)$$
#'
#' 2. Sample fake data $$\tilde{y} \sim p(y|\tilde{\theta})$$
#'
#' 3. Obtain $L$ posterior samples $$\{\theta_1, ..., \theta_L\} \sim p(\tilde{\theta}|\tilde{y})$$
#'
#' 4. Calculate the *rank* $r_s$ of the prior draw $\tilde{\theta}$ wrt to
#' the posterior sample $\{\theta_1, ..., \theta_L\} \sim p(\tilde{\theta}|\tilde{y})$ which falls into the range $[0,L]$
#' out of the possible $L+1$ ranks. The rank is calculated as
#' $$r_s = \sum_{l=1}^L \mathbb{I}[ \theta_l < \tilde{\theta}]$$
#'
#' The $S$ ranks then form a uniform $0-1$ density and the count in
#' each bin has a binomial distribution with probability of
#' $$p(r \in \mbox{Any Bin}) =\frac{(L+1)}{S}.$$
#'
#' ## Model description TODO
#'
#' The fake data simulation function returns ... TODO. Please refer to
#' the `sbc_tools.R` and `make_reference_rankhist.R` R programs for the
#' implementation details.
#'
#' The reference runs are created with $L=1023$ posterior draws for
#' each replication and a total of $S=10^4$ replications are run per
#' case. For the evaluation here the results are reduced to
#' $B=L'+1=64$ bins to ensure a sufficiently large sample size per
#' bin.
#'

calibration <- readRDS("calibration.rds")


include_plots <- TRUE
if("params" %in% ls())
    include_plots <- params$include_plots

# The summary function we use here scales down the $L+1=1024$ bins to
# smaller number of rank bins. This improves the number of counts
# expected per rank bin ($S/(L+1)$) and is thus more robust in terms
# of large number laws. We choose $L=1023$ samples from the posterior
# such that we have $1024 = 2^10$ bins for the ranks. Thus any power
# of $2$ can be used to scale down the number of bins.

plot_binned <- function(cal_df) {

    pl <- NULL

    if(!include_plots)
        return(pl)
    
    if(!all(cal_df$count == 0)) {
        
        S <- calibration$S
        B <- calibration$B
        
        c95 <- qbinom(c(0.025, 0.5, 0.975), S, 1 / B)
        
        dd <- cal_df %>%
            arrange(param, bin) %>%
                group_by(param) %>%
                    mutate(ecdf = cumsum(count) / S, ecdf_ref = (bin + 1) / B) %>%
                        filter(!all(ecdf == 0))
        
        nparam <- length(unique(dd$param))
        if(unique(dd$partype) %in% c("mu_eta", "tau_eta")){
            nc <- nparam
        } else{
            nc <- 2
        }
        
        nr <- max(1, ceiling(nparam / nc))

        pl <- list()
        pl[["hist"]] <- ggplot(dd, aes(bin, count)) +
            facet_wrap(~ param, nrow = nr, ncol = nc) +
                geom_col() +
                    geom_hline(yintercept=c95[c(1,3)], linetype=I(2)) +
                        geom_hline(yintercept=c95[c(2)], linetype=I(3))
        pl[["ecdf_diff"]] <- ggplot(dd, aes(bin, ecdf-ecdf_ref)) +
            facet_wrap(~ param, nrow = nr, ncol = nc) +
                geom_step() +
                    geom_hline(yintercept=0, linetype=I(3))
        pl
    }

    return(pl)
}


B <- calibration$B
S <- calibration$S

bins_all <- calibration$data %>%
  tidyr::gather(key = "param", value = "count", - model, -bin) %>%
  mutate(partype = sapply(strsplit(param, "[[]"), '[', 1),
         group = interaction(model, partype))


cal_split <- split(bins_all, bins_all$group)

pl_split <- lapply(cal_split, function(cal_df) plot_binned(cal_df))


#' # SBC results
#'
#' ## Sampler Diagnostics Overview
#'
kable(calibration$sampler_diagnostics, digits=3)

chisq  <- bins_all %>%
  arrange(model, partype, param, bin) %>%
  group_by(model, partype, param) %>%
  mutate(allna = all(count == 0)) %>%
  filter(!allna) %>%
  do(tidy(chisq.test(.$count))[,c(1,3,2)] ) %>%
  rename(df = parameter) %>%
  ungroup()

#'
#' Large Rhat is defined as exceeding $1.2$.
#' 

#'
#' ## $\chi^2$ Statistic, Model 1: Single-agent logistic regression
#'

kable(chisq %>% filter(model == "log2bayes_EXNEX") %>% select(-model, -partype), digits=3)

#'
#' ## $\chi^2$ Statistic, Model 2: Double combination, fully exchangeable
#'


kable(chisq %>% filter(model == "combo2_EX") %>% select(-model, -partype), digits=3)

#'
#' ## $\chi^2$ Statistic, Model 3: Double combination, EXchangeable/NonEXchangeable model
#'

kable(chisq %>% filter(model == "combo2_EXNEX") %>% select(-model, -partype), digits=3)

#'
#' ## $\chi^2$ Statistic, Model 4: Triple combination, EX/NEX model
#'

kable(chisq %>% filter(model == "combo3_EXNEX") %>% select(-model, -partype), digits=3)

#+ results="asis", include=include_plots, eval=include_plots
spin_child("sbc_report_plots.R")

#'
#' ## Session Info
#'
sessionInfo()

