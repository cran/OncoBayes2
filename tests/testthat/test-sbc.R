context("SBC validation")

test_that("blrm_exnex meets SBC requirements wrt to a Chi-Square statistic.", {
  require(dplyr)
  require(tidyr)
  sbc_chisq_test <- OncoBayes2:::calibration_data |>
    group_by(data_scenario, param) |>
    do(as.data.frame(chisq.test(.$count)[c("statistic", "p.value")])) |>
    ungroup()
  num_tests <- nrow(sbc_chisq_test)
  num_failed <- sum(sbc_chisq_test$p.value < 0.05)
  expect_integer(
    num_failed,
    lower = qbinom(0.025, num_tests, 0.05),
    upper = qbinom(0.975, num_tests, 0.05),
    any.missing = FALSE
  )
})

test_that("blrm_exnex meets SBC requirements per bin.", {
  require(dplyr)
  require(tidyr)
  B <- OncoBayes2:::calibration_meta$B
  S <- OncoBayes2:::calibration_meta$S
  alpha <- 0.2
  ptrue <- 1 / B
  crit_low <- qbinom(alpha / 2, S, ptrue)
  crit_high <- qbinom(1 - alpha / 2, S, ptrue)
  sbc_binom_test <- OncoBayes2:::calibration_data |>
    summarise(
      crit = sum(count < crit_low | count > crit_high),
      .by = c(data_scenario, param)
    ) |>
    mutate(
      pvalue = pbinom(crit, B, alpha),
      extreme = pvalue < 0.025 | pvalue > 0.975
    )
  num_tests <- nrow(sbc_binom_test)
  num_failed <- sum(sbc_binom_test$extreme)
  expect_integer(
    num_failed,
    lower = qbinom(0.025, num_tests, 0.05),
    upper = qbinom(0.975, num_tests, 0.05),
    any.missing = FALSE
  )
})

test_that("blrm_exnex data was up to date at package creation.", {
  calibration_datum <- OncoBayes2:::calibration_meta$created
  package_datum <- OncoBayes2:::pkg_create_date
  delta <- difftime(package_datum, calibration_datum, units = "weeks")
  expect_true(delta < 52. / 2.)
})
