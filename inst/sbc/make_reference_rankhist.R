#! /usr/bin/env Rscript

start_time <- Sys.time()

here::i_am("inst/sbc/make_reference_rankhist.R")
library(here)

## ensure that the current dev version of OncoBayes2 is build and
## loaded.  code to use if one wants the more conservative
## load_OB2_dev_install routine to be used
setwd(here())
if (system("make dev-install") != 0) {
  stop("make dev-install failed!")
}

setwd(here("inst", "sbc"))

# Always reload startup script when working interactively
options(sbc.job.startup.complete = FALSE)
source(here("inst", "sbc", "sbc_job_startup.R"))

library(clustermq)
library(knitr)
set.seed(453453)

## replications to use
S <- 1E4
##S <- 1E3
##S <- 2E3
##S <- 5E2
## S <- 10

examples_set <- names(example_models)
## temporarily reduced set of examples
##examples_set <- c("combo2_EX", "log2bayes_EXNEX")
##examples_set <- c("log2bayes_EXNEX")
##examples_set <- c("combo3_EXNEX", "log2bayes_EXNEX")

scenarios <- expand.grid(repl = 1:S, data_scenario = examples_set, stringsAsFactors = FALSE)
scenarios <- cbind(job.id = 1:nrow(scenarios), scenarios)

num_simulations <- nrow(scenarios)

cat("Total number of jobs to dispatch:", num_simulations, "\n")

scheduler <- getOption("clustermq.scheduler")

if (is.null(scheduler)) {
  ## in this case we enable the multiprocess option to leverage local CPUs
  options(clustermq.scheduler = "multiprocess")
}

## used for debugging
## options(clustermq.scheduler="LOCAL")

scheduler <- getOption("clustermq.scheduler")

n_jobs <- 1
if (scheduler == "multiprocess") {
  ## on a local machine we only use as many CPUs as available
  n_jobs <- parallelly::availableCores(methods = "system")
}
if (scheduler %in% c("LSF", "SGE", "SLURM", "PBS", "Torque")) {
  ## on a queuing enabled backend, we use a lot more parallelism
  n_jobs <- 200
  ## n_jobs  <- 10
}

# Only use as many jobs as there are simulations at most
n_jobs <- min(n_jobs, num_simulations)

cat("Using clustermq backend", scheduler, "with", n_jobs, "concurrent jobs.\n")


RNGkind("L'Ecuyer-CMRG")
set.seed(56969839)
rng_seeds <- setup_lecuyer_seeds(.Random.seed, num_simulations)

scenarios$rng_seed <- rng_seeds

S_warmup <- max(floor(0.1 * S), 1)
S_final <- S - S_warmup

idx_warmup <- scenarios$repl <= S_warmup

scenarios_warmup <- scenarios[idx_warmup, ]
scenarios_main <- scenarios[!idx_warmup, ]

worker <- clustermq::workers(n_jobs, reuse = TRUE, template = list(walltime = 300, job_name = "ob2_sbc", memory = 6 * 1024), log_worker = FALSE)

if (FALSE) {
  ## for debugging purposes
  res <- list()
  for (row in seq_len(nrow(scenarios_warmup))) {
    res[[row]] <- with(scenarios_warmup, run_sbc_case(job.id[row], repl[row], data_scenario[row], example_models, rng_seeds))
  }
  res[[2]]
}

sim_result_warmup <- Q_rows(
  scenarios_warmup,
  run_sbc_case,
  const = list(example_models_cmq = example_models),
  n_jobs = n_jobs,
  workers = worker
)


## TODO: crew alternative? Should be implemented as we then get much better
## and easier debugging!!! Easiest debugging is via targets...
# if(FALSE) {
#
#     library(crew)
#
#     dir.create(here("log"), FALSE)
#     controller <- crew_controller_local(
#         name = "sbc",
#         workers = 10,
#         seconds_idle = 10,
#         local_log_directory=here("log")
#     )
#
#     sim_result_warmup <- controller$map(sbc_tools$run_sbc_case(job.id, repl, data_scenario, base_scenarios, seeds),
#                                         iterate=scenarios_warmup,
#                                         data=list(base_scenarios=sbc_tools$example_models, seeds=rng_seeds),
#                                         globals=c(as.list(sbc_tools), list(sbc_tools=sbc_tools)),
#                                         library=.libPaths(),
#                                         packages=pkg)
#
#     controller$error
#     controller$log
# }

assert_that(sum(idx_warmup) == length(sim_result_warmup), msg = "Check if all warmup simulations were processed.")

#'
#' Collect warmup info
#'

warmup_info <- lapply(sim_result_warmup, function(run) {
  run[c("stepsize", "inv_metric", "draws")]
})

scenarios_warmup <- transform(scenarios_warmup, warmup_job.id = 1:nrow(scenarios_warmup))

warmup_info_by_model <- lapply(
  split(scenarios_warmup$warmup_job.id, scenarios_warmup$data_scenario),
  function(mj) list(warmup_info = warmup_info[mj])
)

#'
#' Schedule remaining runs with available warmup info
#'

example_models_with_warmup <- modifyList(example_models, warmup_info_by_model)

sim_result_main <- Q_rows(
  scenarios_main,
  run_sbc_case,
  const = list(example_models_cmq = example_models_with_warmup),
  n_jobs = n_jobs,
  workers = worker
)

assert_that(sum(!idx_warmup) == length(sim_result_main), msg = "Check if all main simulations were processed.")

## shut down cluster worker
## if proper cleanup is successful, cancel kill-on-exit
## see https://mschubert.github.io/clustermq/articles/technicaldocs.html

if (worker$cleanup()) {
  on.exit()
}

extract_results <- function(x) {
  c(
    rank = as.list(x$rank),
    list(
      job.id = x$job.id,
      time.running = x$time.running,
      n_divergent = x$n_divergent,
      min_Neff_bulk = x$min_Neff_bulk,
      min_Neff_tail = x$min_Neff_tail,
      max_Rhat = x$max_Rhat,
      lp_ess_bulk = x$lp_ess_bulk,
      lp_ess_tail = x$lp_ess_tail,
      stepsize = x$stepsize,
      accept_stat = x$accept_stat
    )
  )
}

calibration_data <- bind_rows(lapply(sim_result_warmup, extract_results), lapply(sim_result_main, extract_results)) %>%
  arrange(job.id) %>%
  merge(scenarios, by = "job.id")

## check that indeed all jobs have finished
assert_that(nrow(calibration_data) == num_simulations)

## collect sampler diagnostics
sampler_diagnostics_summary <- calibration_data %>%
  group_by(data_scenario) %>%
  summarize(
    N = n(),
    total_divergent = sum(n_divergent),
    min_ess_bulk = min(min_Neff_bulk),
    min_ess_tail = min(min_Neff_tail),
    max_Rhat = max(max_Rhat),
    total_large_Rhat = sum(max_Rhat > 1.1),
    min_lp_ess_bulk = min(lp_ess_bulk),
    min_lp_ess_tail = min(lp_ess_tail)
  )


cat("\nSampler diagnostics:\n\n")
kable(sampler_diagnostics_summary, digits = 3)
cat("\n")

if (sum(sampler_diagnostics_summary$total_divergent) != 0) {
  warning("There were some divergent transitions!")
}
if (any(sampler_diagnostics_summary$max_Rhat > 1.1)) {
  warning("There were some parameters with large Rhat!")
}

#' Bin raw data as used in the analysis.
B <- 1024L / 2^5

rank_params <- names(calibration_data)[grepl(names(calibration_data), pattern = "rank")]

calibration_data_binned <- calibration_data %>%
  mutate(across(starts_with("rank"), function(x) ceiling((x + 1) / (1024 / B) - 1)))

## head(calibration_data_binned)

names(calibration_data_binned) <- gsub(
  names(calibration_data_binned),
  pattern = "rank[.]",
  replacement = ""
)

params <- gsub(rank_params, pattern = "rank[.]", replacement = "")

calibration_binned <- calibration_data_binned %>%
  dplyr::select(-job.id, -n_divergent, -min_Neff_bulk, -min_Neff_tail, -rng_seed) %>%
  tidyr::gather(key = "param", value = "bin", -data_scenario) %>%
  group_by(data_scenario, param, bin) %>%
  tally() %>%
  right_join(
    expand.grid(
      data_scenario = unique(calibration_data_binned$data_scenario),
      param = params,
      bin = 0:(B - 1),
      stringsAsFactors = FALSE
    ),
    c("data_scenario", "param", "bin")
  ) %>%
  replace_na(list(n = 0)) %>%
  arrange(data_scenario, param, bin) %>%
  spread(key = param, value = n)



#' Save as data.frame to avoid data.table dependency.
calibration_data <- as.data.frame(calibration_data)
calibration_binned <- as.data.frame(calibration_binned)

#' Further identification and verification data of run
git_hash <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
created <- Sys.time()
created_str <- format(created, "%F %T %Z", tz = "UTC")

calibration <- list( ## raw = calibration_data, ## stop storing raw results, which are not needed for SBC reports
  data = calibration_binned,
  sampler_diagnostics = sampler_diagnostics_summary,
  S = S,
  B = B,
  git_hash = git_hash,
  created = created
)

saveRDS(calibration, file = here("inst", "sbc", "calibration.rds"))
saveRDS(calibration_data, file = here("inst", "sbc", "calibration_data.rds"))

library(tools)
md5 <- md5sum(here("inst", "sbc", "calibration.rds"))
cat(paste0("Created:  ", created_str, "\ngit hash: ", git_hash, "\nMD5:      ", md5, "\n"),
  file = here("inst", "sbc", "calibration.md5")
)


#'
#' Summarize execution time
#'

job_report <- calibration_data[c("job.id", "repl", "data_scenario", "time.running")]
job_report$time.running <- job_report$time.running / 60
job_report$phase <- factor(idx_warmup, c(TRUE, FALSE), c("warmup", "main"))

runtime_by_problem_phase <- job_report %>%
  group_by(data_scenario, phase) %>%
  summarize(total = sum(time.running), mean = mean(time.running), max = max(time.running)) %>%
  pivot_wider(id_cols = "data_scenario", names_from = "phase", values_from = c("total", "mean", "max"))

runtime_by_problem <- job_report %>%
  group_by(data_scenario) %>%
  summarize(total = sum(time.running), mean = mean(time.running), max = max(time.running))

cat("Summary on job runtime on cluster:\n\n")

cat("\nRuntime by data scenario and phase:\n")
print(kable(runtime_by_problem_phase, digits = 2))

cat("\nRuntime by data scenario:\n")
print(kable(runtime_by_problem, digits = 2))

end_time <- Sys.time()

total_runtime <- difftime(end_time, start_time)
units(total_runtime) <- "mins"

cat("\n\nTotal runtime (min):", as.numeric(total_runtime), "\n\n\n")


#' Session info
print(sessionInfo())
