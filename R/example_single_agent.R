#' @title Single Agent Example
#'
#' @name example-single-agent
#'
#' @description
#' Example using a single experimental drug.
#'
#' @details
#' The single agent example is described in the reference
#' Neuenschwander, B. et al (2008). The data are described
#' in the help page for \code{hist_SA}. In this case, the data
#' come from only one study, with the treatment being only single
#' agent. Hence the model specified does not involve a hierarchical
#' prior for the intercept and log-slope parameters. The model
#' described in Neuenschwander, et al (2008) is adapted as follows:
#' \deqn{\text{logit}\, \pi(d) = \log\, \alpha + \beta \, \log\, \Bigl(\frac{d}{d^*}\Bigr),}
#' where \eqn{d^* = 250}, and the prior for
#' \eqn{\boldsymbol\theta = (\log\, \alpha, \log\, \beta)} is
#' \deqn{\boldsymbol\theta \sim \text{N}(\boldsymbol m, \boldsymbol S),}
#' and \eqn{\boldsymbol m = (\text{logit}\, 0.5, \log\, 1)} and
#' \eqn{\boldsymbol S = \text{diag}(2^2, 1^2)} are constants.
#'
#' The above model is non-hierarchical. To disable the hierarchical
#' model structure of the \code{blrm_exnex} framework, the user can
#' specify the option \code{prior_tau_dist=NULL}. This will internally
#' set all the heterogeniety parameters (\eqn{\tau^2_\alpha} and
#' \eqn{\tau^2_\beta}) to zero.
#'
#' @template ref-critical_aspects
#'
#' @template start-example
#' @template example-single_agent
#' @template stop-example
#'
NULL
