% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update.R
\name{update.blrm_trial}
\alias{update.blrm_trial}
\title{Update data and/or prior of a BLRM trial}
\usage{
\method{update}{blrm_trial}(object, ...)
}
\arguments{
\item{object}{blrm_trial object}

\item{...}{passed to default \code{update} command of \code{blrm_exnex}}
}
\description{
\itemize{
\item Adds data rows to a \code{\link{blrm_trial}} object (add_data argument)
\item Replaces data of a \code{\link{blrm_trial}} object (data argument)
\item Sets the prior of a \code{\link{blrm_trial}} object (... argument will be
passed to blrm_exnex)
}
}
\examples{
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)


# the combo2_trial example demonstrates the use of add_data of
# update.blrmfit
example_model("combo2_trial")

## Recover user set sampling defaults
options(.user_mc_options)

}
