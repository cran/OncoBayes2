if (!getOption("sbc.job.startup.complete", FALSE)) {
  here::here("inst/sbc/sbc_job_startup.R")

  pkg <- c("assertthat", "rstan", "mvtnorm", "checkmate", "Formula", "abind", "dplyr", "tidyr", "posterior", "here", "bayesplot")
  pkg_load <- sapply(pkg, library, character.only = TRUE)

  source(here("inst", "sbc", "sbc_tools.R"))
  source(here("inst", "sbc", "lkj.R"))

  assert_that(dir.exists(here("build", "installed", "OncoBayes2")))
  ob2_lib_dir <- here("build", "installed")

  load_OB2_dev(ob2_lib_dir)
  assert_that(.libPaths()[1] == ob2_lib_dir)

  source(here("inst", "sbc", "sbc_example_models.R"))
  options(sbc.job.startup.complete = TRUE)
  cat("SBC startup completed.\n")
} else {
  cat("SBC startup already completed, doing nothing.\n")
}
