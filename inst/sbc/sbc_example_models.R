#'
#' Utilities for SBC validation
#'

example_designs <- list(
  combo3_EXNEX = list(
    modeltype = "EXNEX",
    EXNEX_comp = c(TRUE, TRUE, FALSE),
    EX_prob_comp = matrix(c(
      0.5, 0.9, 1.0,
      1.0, 0.8, 1.0,
      0.5, 0.5, 1.0
    ), byrow = TRUE, nrow = 3, ncol = 3),
    design = expand.grid(
      stratum_id = "STRAT",
      group_id = LETTERS[1:3],
      drug_A = c(0.25, 0.5, 1),
      drug_B = c(0.25, 0.5, 1),
      drug_C = c(0.25, 0.5, 1),
      num_patients = 3
    ) |> arrange(stratum_id, group_id, drug_A, drug_B, drug_C),
    dref = c(1, 2, 4)
  ),
  combo2_EX = list(
    modeltype = "EX",
    EXNEX_comp = c(FALSE, FALSE),
    EX_prob_comp = matrix(c(1.0), nrow = 4, ncol = 2),
    design = bind_rows(expand.grid(
      stratum_id = "STRAT1",
      group_id = LETTERS[1:2],
      drug_A = c(0.25, 0.5, 1),
      drug_B = c(0.25, 0.5, 1),
      num_patients = 5
    ), expand.grid(
      stratum_id = "STRAT2",
      group_id = LETTERS[3:4],
      drug_A = c(0.25, 0.5, 1),
      drug_B = c(0.25, 0.5, 1),
      num_patients = 5
    )) |> arrange(stratum_id, group_id, drug_A, drug_B),
    dref = c(1, 2)
  ),
  combo2_EXNEX = list(
    modeltype = "EXNEX",
    EXNEX_comp = c(TRUE, TRUE),
    EX_prob_comp = matrix(c(0.8), nrow = 3, ncol = 2),
    design = expand.grid(
      stratum_id = "STRAT",
      group_id = LETTERS[1:3],
      drug_A = c(0.25, 0.5, 1),
      drug_B = c(0.25, 0.5, 1),
      num_patients = 5
    ) |> arrange(stratum_id, group_id, drug_A, drug_B),
    dref = c(1, 2)
  ),
  log2bayes_EXNEX = list(
    modeltype = "EXNEX",
    EXNEX_comp = c(TRUE),
    EX_prob_comp = matrix(c(0.1, 0.6, 1.0), nrow = 3, ncol = 1),
    design = expand.grid(
      stratum_id = c("STRAT1"),
      group_id = LETTERS[1:3],
      drug_A = c(0.0625, 0.125, 0.25, 0.5, 1),
      num_patients = 5
    ) |>
      arrange(stratum_id, group_id, drug_A),
    dref = c(1)
  ),
  log2bayes_EX = list(
    modeltype = "EX",
    EXNEX_comp = c(FALSE),
    EX_prob_comp = matrix(c(1.0, 1.0, 1.0), nrow = 3*2, ncol = 1),
    design = bind_rows(expand.grid(
      stratum_id = c("STRAT1"),
      group_id = LETTERS[1:3],
      drug_A = c(0.0625, 0.125, 0.25, 0.5, 1),
      num_patients = 5
    ), expand.grid(
      stratum_id = c("STRAT2"),
      group_id = LETTERS[4:6],
      drug_A = c(0.0625, 0.125, 0.25, 0.5, 1),
      num_patients = 5
    )) |>
      arrange(stratum_id, group_id, drug_A),
    dref = c(1)
  )
)

example_models <- lapply(
  ## example_designs[c("combo2_EX", "log2bayes_EXNEX")], ## faster set for testing only
  ## example_designs[c("log2bayes_EXNEX")], ## faster set for testing only
  ## example_designs[c("combo2_EXNEX", "log2bayes_EXNEX")], ## faster set for testing only
  example_designs,
  function(example) {
    design <- example$design
    dref <- example$dref
    modeltype <- example$modeltype
    ## is_exnex <- modeltype != "EX"
    ## p_exch <- switch(
    ##    modeltype,
    ##    "EX" = 1,
    ##    "NEX" = 1e-06,
    ##    "EXNEX" = 0.8
    ## )
    example_prior_EX_prob_comp <- example$EX_prob_comp
    example_prior_is_EXNEX_comp <- example$EXNEX_comp


    names(dref) <- grep(names(design), pattern = "drug_", value = TRUE)

    linear_formula <- blrm_formula_linear(dref, 3)

    num_comp <- linear_formula$num_components
    num_inter <- linear_formula$num_interaction_terms
    formula <- as.formula(linear_formula$blrm_formula)

    num_groups <- nlevels(design$group_id)
    num_strata <- nlevels(design$stratum)

    prior_EX_tau_comp_arg <- list()
    for(s in seq_len(num_strata)) {
      prior_EX_tau_comp_arg[[s]] <- replicate(num_comp, mixmvnorm(c(1,
                                                                    log(c(0.25, 0.125) / (2.0 * s)),
                                                                    diag(c(log(0.5 + s) / 1.96, log(0.5 + s) / 1.96)^2))), FALSE)  
    }

    blrm_args <- list(
      formula = formula,
      data = design,
      prior_EX_mu_comp = replicate(num_comp, mixmvnorm(c(1, logit(1 / 3), 0, diag(c(1, 0.5)^2))), FALSE),
      prior_EX_tau_comp = prior_EX_tau_comp_arg,
      prior_EX_corr_eta_comp = rep(2.0, num_comp),
      prior_EX_corr_eta_inter = 2.0,
      prior_EX_prob_comp = example_prior_EX_prob_comp,
      prior_EX_prob_inter = matrix(1.0, nrow = num_groups, ncol = num_inter),
      prior_is_EXNEX_comp = example_prior_is_EXNEX_comp,
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_tau_dist = 1,
      ## SW: to actually test NEX  make the distribution different!
      prior_NEX_mu_comp = replicate(num_comp, mixmvnorm(c(1, logit(2 / 3), 0, diag(c(1, 0.5)^2))), FALSE),
      iter = 1000 + 1000,
      warmup = 1000,
      ## iter = 200 + 600,
      ## warmup = 200,
      control = list(
        stepsize = 0.5,
        adapt_init_buffer = 75,
        adapt_window = 25,
        adapt_term_buffer = 4 * 50, ## run longer as normal terminal window
        adapt_delta = 0.8
      ),
      ## iter = 150,
      ## warmup = 50,
      thin = 1,
      init = 1.0,
      chains = 2,
      ## cores = 1, ## control via mc.cores option
      prior_PD = FALSE,
      save_warmup = FALSE
    )

    if(num_inter > 0) {
      blrm_args$prior_EX_mu_inter <- mixmvnorm(c(1, rep.int(0, num_inter), diag((log(2)/1.96)^2, num_inter, num_inter)))
      blrm_args$prior_EX_tau_inter <- replicate(num_strata, mixmvnorm(c(1, rep.int(log(0.25), num_inter), diag((log(2)/1.96)^2, num_inter, num_inter))), FALSE)
      blrm_args$prior_NEX_mu_inter <- mixmvnorm(c(1, rep.int(0, num_inter), diag((log(2)/1.96)^2, num_inter, num_inter)))
    }

    base_args <- blrm_args
    base_args$data <- base_args$data |> mutate(num_toxicities = 0)
    base_args$iter <- 2
    base_args$warmup <- 1
    base_fit <- do.call(blrm_exnex, base_args)
    ## blrm arguments used when warmup info of stepsize and the
    ## inverse metric is being provided
    blrm_args_with_warmup_info <- modifyList(
      blrm_args,
      list(
        warmup = 1000, iter = 1000 + 1000,
        init = NULL,
        control = list(
          adapt_init_buffer = 75, ## default
          adapt_window = 25, ## default
          adapt_term_buffer = 4 * 50, ## longer final terminal buffer
          adapt_delta = 0.8
        )
      )
    )

    return(
      list(
        dref = dref,
        num_strata = num_strata,
        num_groups = num_groups,
        num_comp = num_comp,
        num_inter = num_inter,
        blrm_args = blrm_args,
        blrm_args_with_warmup_info = blrm_args_with_warmup_info,
        base_fit = base_fit
      )
    )
  }
)
