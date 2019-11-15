## Examples used for blrm_trial tests

examples <- list(

  # Single-agent example -------------------------------------------------------
  single_agent = list(
    histdata = bind_rows(list(
      tibble(
        group_id = "hist_A",
        drug1 = 1:2,
        num_patients = 3,
        num_toxicities = 0
      ),
      tibble(
        group_id = "hist_B",
        drug1 = 1:2,
        num_patients = 3,
        num_toxicities = 0
      )
    )),
    dose_info = bind_rows(list(
      tibble(
        group_id = "cur_A",
        drug1 = 1:2
      ),
      tibble(
        group_id = "cur_B",
        drug1 = 1:2
      )
    )),
    drug_info = tibble(
      drug_name = "drug1",
      dose_ref  = 1,
      dose_unit = "ngogn",
      reference_p_dlt = 0.1
    )
  ),

  # Combo2 example -------------------------------------------------------
  combo2 = list(
    histdata = bind_rows(list(
      tibble(
        group_id = "hist_A",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        num_patients = 3,
        num_toxicities = 0
      ),
      tibble(
        group_id = "hist_B",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        num_patients = 3,
        num_toxicities = 0
      )
    )),
    dose_info = bind_rows(list(
      tibble(
        group_id = "cur_A",
        drug1 = 1:2,
        drug2 = 100 * (1:2)
      ),
      tibble(
        group_id = "cur_B",
        drug1 = 1:2,
        drug2 = 100 * (1:2)
      )
    )),
    drug_info = tibble(
      drug_name = c("drug1", "drug2"),
      dose_ref  = c(1, 100),
      dose_unit = c("ngogn", "potrzebie"),
      reference_p_dlt=0.3
    )
  ),

  # Combo3 example -------------------------------------------------------
  combo3 = list(
    histdata = bind_rows(list(
      tibble(
        group_id = "hist_A",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        drug3 = 1000 * (1:2),
        num_patients = 3,
        num_toxicities = 0
      ),
      tibble(
        group_id = "hist_B",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        drug3 = 1000 * (1:2),
        num_patients = 3,
        num_toxicities = 0
      )
    )),
    dose_info = bind_rows(list(
      tibble(
        group_id = "cur_A",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        drug3 = 1000 * (1:2)
      ),
      tibble(
        group_id = "cur_B",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        drug3 = 1000 * (1:2)
      )
    )),
    drug_info = tibble(
      drug_name = paste0("drug", 1:3),
      dose_ref  = 10 ^ c(0, 2, 3),
      dose_unit = c("ngogn", "potrzebie", "blintz")
    )
  ),

  single_drug_with_strata = list(
    histdata = bind_rows(list(
      tibble(
        group_id = "hist_A",
        stratum_id = "strat_A",
        drug1 = 1:2,
        num_patients = 3,
        num_toxicities = 0
      ),
      tibble(
        group_id = "hist_B",
        stratum_id = "strat_B",
        drug1 = 1:2,
        num_patients = 3,
        num_toxicities = 0
      )
    )),
    dose_info = bind_rows(list(
      tibble(
        group_id = "cur_A",
        stratum_id = "strat_A",
        drug1 = 1:2
      ),
      tibble(
        group_id = "cur_B",
        stratum_id = "strat_B",
        drug1 = 1:2
      )
    )),
    drug_info = tibble(
      drug_name = "drug1",
      dose_ref  = 1,
      dose_unit = "ngogn"
    )
  )
)