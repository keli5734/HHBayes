# Suppress R CMD check NOTEs for dplyr/data.table NSE variables
utils::globalVariables(c(
  ":=", "N_pop", "Primary_AR", "Reinf_Rate_Among_Infected", "Reinf_Rate_Pct",
  "alpha_val", "bin_start_idx", "cases", "count_in_role", "covariate", "ct_value",
  "date_bin", "date_formatted", "date_infection", "date_resolved_final", "day_idx",
  "day_index", "days_since_inf", "episode_id", "episode_start", "familyidstars",
  "first_inf", "group", "hh_id", "hh_id_int", "hh_sort_num", "imputed_val",
  "infection_time", "label", "lbl_x", "lbl_x_num", "lbl_y", "log_value", "med",
  "n_id", "n_infections", "n_primaries", "n_reinfections", "n_total_episodes",
  "n_unique_people", "ndist", "ndx", "ndy", "nfinal_x", "nfinal_y", "nmx", "nmy",
  "noffset", "nshift_x", "nshift_y", "nx1", "nx2", "ny1", "ny2", "p_id",
  "p_sort_num", "participantid", "pcr_sample", "person_id", "prob", "rank_in_role",
  "role", "role_group", "role_id", "role_name", "source_int", "source_role",
  "source_x", "source_y", "start_risk_date", "target", "target_x", "target_y",
  "test_result", "total", "type", "val", "valid", "value", "x1", "x2",
  "y1", "y2", "y_plot"
))
