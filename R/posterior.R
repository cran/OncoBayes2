# Code adapted from
# https://raw.githubusercontent.com/paul-buerkner/brms/refs/heads/master/R/posterior.R

#'
#' Transform \code{blrmfit} or \code{blrm_trial} to \code{draws} objects
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' 
#' Transform a \code{blrmfit} or \code{blrm_trial} object to a format supported by the
#' \pkg{posterior} package.
#'
#' @aliases as_draws as_draws_matrix as_draws_array as_draws_df
#' @aliases as_draws_rvars as_draws_list
#'
#' @param x A \code{blrmfit} or \code{blrm_trial} object.
#' @param variable A character vector providing the variables to
#'   extract.  By default, all variables are extracted.
#' @param regex Logical; Should variable be treated as a (vector of)
#'   regular expressions? Any variable in \code{x} matching at least
#'   one of the regular expressions will be selected. Defaults to
#'   \code{FALSE}.
#' @param inc_warmup Should warmup draws be included? Defaults to
#'   \code{FALSE}.
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @details To subset iterations, chains, or draws, use the
#'   \code{\link[posterior:subset_draws]{subset_draws}} method after
#'   transforming the input object to a \code{draws} object.
#'
#' The function is experimental as the set of exported posterior
#' variables are subject to updates.
#'
#' @seealso \code{\link[posterior:draws]{draws}}
#'   \code{\link[posterior:subset_draws]{subset_draws}}
#'
#' @template start-example
#' @examples
#' # fit an example model. See documentation for "combo2" example
#' example_model("combo2")
#'
#' post <- as_draws(blrmfit)
#'
#' @template stop-example
#'
#' @name draws-OncoBayes2
NULL

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws
#' @method as_draws blrmfit
#' @export
#' @export as_draws
as_draws.blrmfit <- function(x, variable = NULL, regex = FALSE,
                             inc_warmup = FALSE, ...) {
  # draws_list is the fastest format to convert to at the moment
  .as_draws_conversion(x$stanfit,
                       as_draws_list, variable = variable,
                       regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix blrmfit
#' @export
#' @export as_draws_matrix
as_draws_matrix.blrmfit <- function(x, variable = NULL, regex = FALSE,
                                    inc_warmup = FALSE, ...) {
  .as_draws_conversion(x$stanfit,
                       as_draws_matrix, variable = variable,
                       regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_array
#' @method as_draws_array blrmfit
#' @export
#' @export as_draws_array
as_draws_array.blrmfit <- function(x, variable = NULL, regex = FALSE,
                                   inc_warmup = FALSE, ...) {
  .as_draws_conversion(x$stanfit,
                       as_draws_array, variable = variable,
                       regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_df
#' @method as_draws_df blrmfit
#' @export
#' @export as_draws_df
as_draws_df.blrmfit <- function(x, variable = NULL, regex = FALSE,
                                inc_warmup = FALSE, ...) {
  .as_draws_conversion(x$stanfit,
                       as_draws_df, variable = variable,
                       regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_list
#' @method as_draws_list blrmfit
#' @export
#' @export as_draws_list
as_draws_list.blrmfit <- function(x, variable = NULL, regex = FALSE,
                                  inc_warmup = FALSE, ...) {
  .as_draws_conversion(x$stanfit,
                       as_draws_list, variable = variable,
                       regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars blrmfit
#' @export
#' @export as_draws_rvars
as_draws_rvars.blrmfit <- function(x, variable = NULL, regex = FALSE,
                                   inc_warmup = FALSE, ...) {
  .as_draws_conversion(x$stanfit,
                       as_draws_rvars, variable = variable,
                       regex = regex, inc_warmup = inc_warmup, ...)
}

# in stanfit objects draws are stored in a draws_list-like format
# so converting from there will be most efficient
# may be removed once rstan supports posterior natively
#' @keywords internal
.as_draws_conversion <- function(x, draws_converter, variable = NULL, regex = FALSE,
                                 inc_warmup = FALSE, ...) {
  stopifnot(.is.stanfit(x))
  inc_warmup <- .as_one_logical(inc_warmup)
  if (!length(x@sim$samples)) {
    .stop2("The model does not contain posterior draws.")
  }
  out <- draws_converter(x)
  # first subset variables then remove warmup as removing warmup
  # will take a lot of time when extracting many variables
  out <- subset_draws(out, variable = variable, regex = regex)
  if (nvariables(out) > 0 && !inc_warmup) {
    nwarmup <- x@sim$warmup2[1] %||% 0
    warmup_ids <- seq_len(nwarmup)
    iteration_ids <- posterior::iteration_ids(out)
    if (length(warmup_ids)) {
      iteration_ids <- iteration_ids[-warmup_ids]
    }
    out <- subset_draws(out, iteration = iteration_ids)
  }
  out
}

#' @keywords internal
.is.stanfit <- function(x) {
  inherits(x, "stanfit")
}

# coerce 'x' to a single logical value
#' @keywords internal
.as_one_logical <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.logical(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- .deparse0(s, max_char = 100L)
    .stop2("Cannot coerce '", s, "' to a single logical value.")
  }
  x
}

#' @keywords internal
.stop2 <- function(...) {
  stop(..., call. = FALSE)
}

#' @keywords internal
'%||%' <- function(x, y) {
  # NOTE: This exact same function is part of base R starting with
  # version R 4.5
  if (is.null(x)) x <- y
  x
}

# combine deparse lines into one string
# since R 4.0 we also have base::deparse1 for this purpose
#' @keywords internal
.deparse0 <- function(x, max_char = NULL, ...) {
  out <- collapse(deparse(x, ...))
  if (isTRUE(max_char > 0)) {
    out <- substr(out, 1L, max_char)
  }
  out
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws
#' @method as_draws blrm_trial
#' @export
#' @export as_draws
as_draws.blrm_trial <- function(x, variable = NULL, regex = FALSE,
                                inc_warmup = FALSE, ...) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws(x$blrmfit, variable = variable,
           regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix blrm_trial
#' @export
#' @export as_draws_matrix
as_draws_matrix.blrm_trial <- function(x, variable = NULL, regex = FALSE,
                                    inc_warmup = FALSE, ...) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_matrix(x$blrmfit, variable = variable,
                  regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_array
#' @method as_draws_array blrm_trial
#' @export
#' @export as_draws_array
as_draws_array.blrm_trial <- function(x, variable = NULL, regex = FALSE,
                                   inc_warmup = FALSE, ...) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_array(x$blrmfit, variable = variable,
                 regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_df
#' @method as_draws_df blrm_trial
#' @export
#' @export as_draws_df
as_draws_df.blrm_trial <- function(x, variable = NULL, regex = FALSE,
                                inc_warmup = FALSE, ...) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_df(x$blrmfit, variable = variable,
              regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_list
#' @method as_draws_list blrm_trial
#' @export
#' @export as_draws_list
as_draws_list.blrm_trial <- function(x, variable = NULL, regex = FALSE,
                                  inc_warmup = FALSE, ...) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_list(x$blrmfit, variable = variable,
           regex = regex, inc_warmup = inc_warmup, ...)
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars blrm_trial
#' @export
#' @export as_draws_rvars
as_draws_rvars.blrm_trial <- function(x, variable = NULL, regex = FALSE,
                                   inc_warmup = FALSE, ...) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_rvars(x$blrmfit, variable = variable,
           regex = regex, inc_warmup = inc_warmup, ...)
}

