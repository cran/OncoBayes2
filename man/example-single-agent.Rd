% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_single_agent.R
\name{example-single-agent}
\alias{example-single-agent}
\title{Single Agent Example}
\description{
Example using a single experimental drug.
}
\details{
The single agent example is described in the reference
Neuenschwander, B. et al (2008). The data are described
in the help page for \code{hist_SA}. In this case, the data
come from only one study, with the treatment being only single
agent. Hence the model specified does not involve a hierarchical
prior for the intercept and log-slope parameters. The model
described in Neuenschwander, et al (2008) is adapted as follows:
\deqn{\text{logit}\, \pi(d) = \log\, \alpha + \beta \, \log\, \Bigl(\frac{d}{d^*}\Bigr),}
where \eqn{d^* = 250}, and the prior for
\eqn{\boldsymbol\theta = (\log\, \alpha, \log\, \beta)} is
\deqn{\boldsymbol\theta \sim \text{N}(\boldsymbol m, \boldsymbol S),}
and \eqn{\boldsymbol m = (\text{logit}\, 0.5, \log\, 1)} and
\eqn{\boldsymbol S = \text{diag}(2^2, 1^2)} are constants.

The above model is non-hierarchical. To disable the hierarchical
model structure of the \code{blrm_exnex} framework, the user can
specify the option \code{prior_tau_dist=NULL}. This will internally
set all the heterogeniety parameters (\eqn{\tau^2_\alpha} and
\eqn{\tau^2_\beta}) to zero.
}
\examples{
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

## Example from Neuenschwander, B., et al. (2009). Stats in Medicine

dref <- 50

## Since there is no prior information the hierarchical model
## is not used in this example by setting tau to (almost) 0.
blrmfit <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
      1 + log(drug_A / dref) |
      0 |
      group_id,
  data = hist_SA,
  prior_EX_mu_comp = mixmvnorm(c(1, logit(1 / 2), log(1), diag(c(2^2, 1)))),
  ## Setting prior_tau_dist=NULL disables the hierarchical prior which is
  ## not required in this example as we analyze a single trial.
  prior_tau_dist = NULL,
  prior_PD = FALSE
)
## Recover user set sampling defaults
options(.user_mc_options)

}
\references{
Neuenschwander, B., Branson, M., & Gsponer, T. (2008).
Critical aspects of the Bayesian approach to phase I cancer trials.
\emph{Statistics in medicine}, 27(13), 2420-2439.
}
