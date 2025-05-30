% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictive_interval.R
\name{predictive_interval.blrmfit}
\alias{predictive_interval.blrmfit}
\alias{predictive_interval}
\title{Posterior predictive intervals}
\usage{
\method{predictive_interval}{blrmfit}(object, prob = 0.95, newdata, ...)
}
\arguments{
\item{object}{fitted model object}

\item{prob}{central probability mass to report, i.e. the quantiles
0.5-prob/2 and 0.5+prob/2 are displayed. Multiple central
widths can be specified.}

\item{newdata}{optional data frame specifying for what to predict;
if missing, then the data of the input model \code{object} is
used}

\item{...}{not used in this function}
}
\value{
Matrix with as many rows as the input data set and two
columns which contain the lower and upper quantile
corresponding to the central probability mass \code{prob} for
the number of responses of the predictive distribution.
}
\description{
Posterior predictive intervals of the model.
}
\details{
Reports for each row of the input data set the predictive interval
according to the fitted model.
}
\examples{
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

example_model("single_agent", silent = TRUE)

predictive_interval(blrmfit)

## Recover user set sampling defaults
options(.user_mc_options)

}
