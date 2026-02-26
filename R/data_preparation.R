#' Prepare Data for Stan Model
#'
#' @param df_clean Dataframe with observation data (must contain 'episode_id' from simulation or data).
#' @param surveillance_df Dataframe with columns 'date' and 'cases'.
#' @param role_levels Character vector. Role categories to use. Defaults to \code{c("adult", "infant", "toddler", "elderly")}.
#' @param study_start_date Date. Start date of the study period. Defaults to 2024-07-01.
#' @param study_end_date Date. End date of the study period. Defaults to 2025-07-01.
#' @param seasonal_forcing_list List of numeric vectors (one per role) for seasonal forcing. NULL for no forcing.
#' @param use_vl_data Logical. Whether to include viral load data. Defaults to TRUE.
#' @param covariates_susceptibility Vector of column names to use as covariates for susceptibility.
#' @param covariates_infectivity Vector of column names to use as covariates for infectivity.
#' @param priors List of flexible priors (dist, params).
#' @param recovery_params List of Gamma parameters (shape, scale) for the immunity tail by role.
#' @param imputation_params List of parameters for mechanistic viral curves (Cpeak, r, d, t_peak) by role.
#' @param model_type String, either "empirical" or "ODE".
#' @param ODE_params_list List of ODE parameters (beta, delta, etc.) by role.
#' @param delta Numeric. Household size scaling parameter.
#'   When delta > 0, transmission rates are scaled by (1/max(household_size, 1))^delta.
#'   Should match the value used in simulation. Defaults to 0 (no scaling).
#' @param role_mixing_matrix 4x4 Matrix defining contact weights between roles.
#' @param seed Integer. Random seed for reproducibility. Defaults to 123.
#'
#' @return A named list formatted for input to the Stan model.
#' @export
prepare_stan_data <- function(df_clean,
                              surveillance_df = NULL,
                              role_levels = c("adult", "infant", "toddler", "elderly"),
                              study_start_date = as.Date("2024-07-01"),
                              study_end_date = as.Date("2025-07-01"),
                              seasonal_forcing_list = NULL,
                              use_vl_data = TRUE,
                              use_curve_logic = FALSE,
                              # --- COVARIATE ARGUMENTS ---
                              covariates_susceptibility = NULL,
                              covariates_infectivity = NULL,

                              model_type = "empirical", # "empirical" or "ODE"
                              ODE_params_list = NULL,   # Required if model_type="ODE"

                              delta = 0,

                              # --- CONTACT ARGUMENT ---
                              role_mixing_matrix = NULL,

                              recovery_params = NULL,
                              imputation_params = NULL,
                              priors = list(),
                              seed = 123) {

  if(!is.null(seed)) set.seed(seed)
  T_max <- as.integer(study_end_date - study_start_date) + 1
  if (T_max <= 0) stop("study_end_date must be after study_start_date")

  max_obs_val <- max(df_clean$pcr_sample, na.rm=TRUE)
  is_ct_data <- (is.finite(max_obs_val) && max_obs_val > 15)
  detected_vl_type <- if(is_ct_data) 0 else 1
  default_val <- if(is_ct_data) 45.0 else 0.0

  # =========================================================
  # 1. PARSE FLEXIBLE PRIORS
  # =========================================================

  parse_prior <- function(p_list, def_type, def_params) {
    if(is.null(p_list)) return(list(type=def_type, params=def_params))
    type_int <- 1
    if(!is.null(p_list$dist)) {
      dist <- tolower(p_list$dist)
      if(dist == "normal") type_int <- 1
      else if(dist == "uniform") type_int <- 2
      else if(dist == "lognormal") type_int <- 3
    }
    params <- if(!is.null(p_list$params)) p_list$params else def_params
    list(type=type_int, params=params)
  }

  p_beta1 <- parse_prior(priors$beta1, 1, c(-5, 1))
  p_beta2 <- parse_prior(priors$beta2, 1, c(-5, 1))
  p_alpha <- parse_prior(priors$alpha, 1, c(-6, 2))
  p_cov   <- parse_prior(priors$covariates, 1, c(0, 1))
  p_shape <- parse_prior(priors$gen_shape, 3, c(log(3.0), 0.2))
  p_rate  <- parse_prior(priors$gen_rate,  3, c(log(0.5), 0.2))

  default_vl_midpoint_prior <- if (detected_vl_type == 0) {
    list(dist = "normal", params = c(33.0, 2.0))   # Ct scale
  } else {
    list(dist = "normal", params = c(4.0,  1.0))   # Log10 VL scale
  }

  default_vl_slope_prior <- if (detected_vl_type == 0) {
    list(dist = "normal",    params = c(4.0,  2.0))       # Logistic steepness
  } else {
    list(dist = "lognormal", params = c(log(1.0), 0.5))   # Hill exponent
  }

  p_vl_midpoint <- parse_prior(priors$vl_midpoint,
                               default_vl_midpoint_prior$type,
                               default_vl_midpoint_prior$params)
  p_vl_slope    <- parse_prior(priors$vl_slope,
                               default_vl_slope_prior$type,
                               default_vl_slope_prior$params)


  # =========================================================
  # 2. STANDARDIZE COLUMN NAMES
  # =========================================================
  if("hh_id" %in% names(df_clean) && !"familyidstars" %in% names(df_clean)) df_clean <- df_clean %>% dplyr::rename(familyidstars = hh_id)
  if("person_id" %in% names(df_clean)) df_clean <- df_clean %>% dplyr::mutate(participantid = paste(familyidstars, person_id, sep = "_"))
  if("role" %in% names(df_clean) && !"role_name" %in% names(df_clean)) df_clean <- df_clean %>% dplyr::rename(role_name = role)
  if("pcr_sample" %in% names(df_clean) && !"ct_value" %in% names(df_clean)) df_clean <- df_clean %>% dplyr::rename(ct_value = pcr_sample)
  if("test_result" %in% names(df_clean) && !"is_in_episode" %in% names(df_clean)) df_clean <- df_clean %>% dplyr::rename(is_in_episode = test_result)

  if("day_index" %in% names(df_clean) && !"date" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(date = study_start_date + (day_index - 1))
  }

  # =========================================================
  # 3. IMPUTATION STEP (Interval Sampling + Recovery Tail)
  # =========================================================
  if (!is.null(recovery_params)) {
    params_resolve <- recovery_params
  } else {
    params_resolve <- list(
      adult   = list(shape=2, scale=3),
      infant  = list(shape=2, scale=3),
      toddler = list(shape=2, scale=3),
      elderly = list(shape=2, scale=3)
    )
  }

  participants <- unique(df_clean$participantid)
  imputed_episodes <- list()

  for (pid in participants) {
    p_data <- df_clean %>% dplyr::filter(participantid == pid) %>% dplyr::arrange(date)
    p_role <- p_data$role_name[1]
    if(is.na(p_role) || !p_role %in% role_levels) p_role <- "adult"
    cur_resolve <- params_resolve[[p_role]]
    if(is.null(cur_resolve)) cur_resolve <- list(shape=2, scale=3)

    episode_ids <- if ("episode_id" %in% names(p_data)) unique(p_data$episode_id[p_data$episode_id > 0]) else integer(0)

    for (eid in episode_ids) {
      ep_rows <- p_data %>% dplyr::filter(episode_id == eid)
      if (nrow(ep_rows) == 0) next

      first_pos_date <- min(ep_rows$date)
      last_pos_date  <- max(ep_rows$date)

      prev_neg_dates <- p_data %>% dplyr::filter(date < first_pos_date, episode_id == 0) %>% dplyr::pull(date)
      lower_bound <- if(length(prev_neg_dates) > 0) max(prev_neg_dates) else (first_pos_date - 7)

      next_neg_dates <- p_data %>% dplyr::filter(date > last_pos_date, episode_id == 0) %>% dplyr::pull(date)
      upper_bound <- if(length(next_neg_dates) > 0) min(next_neg_dates) else (last_pos_date + 7)

      possible_starts <- seq(from = lower_bound + 1, to = first_pos_date, by = "day")
      T_inf_new <- if(length(possible_starts) > 1) sample(possible_starts, 1) else first_pos_date

      possible_ends <- seq(from = last_pos_date, to = upper_bound - 1, by = "day")
      T_resolved_new <- if(length(possible_ends) > 1) sample(possible_ends, 1) else last_pos_date

      if(T_resolved_new < T_inf_new) T_resolved_new <- T_inf_new

      resolve_draw <- max(1, round(stats::rgamma(1, shape = cur_resolve$shape, scale = cur_resolve$scale)))
      T_resolved_final <- T_resolved_new + resolve_draw

      if(eid < max(episode_ids)) {
        next_ep_rows <- p_data %>% dplyr::filter(episode_id == (eid+1))
        if(nrow(next_ep_rows) > 0) {
          next_ep_raw_start <- min(next_ep_rows$date)
          if (T_resolved_final >= next_ep_raw_start) {
            T_resolved_final <- next_ep_raw_start - 1
          }
        }
      }

      imputed_episodes[[length(imputed_episodes) + 1]] <- data.frame(
        participantid = pid, episode_id = eid,
        date_infection = T_inf_new,
        date_infectious_start = T_inf_new,
        date_infectious_end = T_resolved_new,
        date_resolved = T_resolved_new,
        date_resolved_final = T_resolved_final,
        stringsAsFactors = FALSE
      )
    }
  }

  if(length(imputed_episodes) > 0) {
    df_imputed <- do.call(rbind, imputed_episodes)
  } else {
    df_imputed <- data.frame(
      participantid = character(), episode_id = integer(),
      date_infection = as.Date(character()), date_infectious_start = as.Date(character()),
      date_infectious_end = as.Date(character()), date_resolved = as.Date(character()),
      date_resolved_final = as.Date(character()), stringsAsFactors = FALSE
    )
  }

  # =========================================================
  # 4. METADATA & MATRIX BUILDING
  # =========================================================

  cols_to_keep <- c("participantid", "familyidstars", "role_name", "person_id")
  all_covs <- unique(c(covariates_susceptibility, covariates_infectivity))
  if (!is.null(all_covs)) {
    missing <- setdiff(all_covs, names(df_clean))
    if(length(missing) > 0) stop(paste("Covariates missing in df_clean:", paste(missing, collapse=", ")))
    cols_to_keep <- c(cols_to_keep, all_covs)
  }


  # --- STEP A: SAFE SORTING (NUMERIC AWARE) ---
  df_meta_all <- df_clean %>%
    dplyr::select(dplyr::all_of(cols_to_keep)) %>%
    dplyr::distinct() %>%
    dplyr::rename(hh_id = familyidstars, role = role_name) %>%
    dplyr::filter(!is.na(role)) %>%
    # Extract the hidden numbers for sorting
    dplyr::mutate(
      hh_sort_num = suppressWarnings(as.numeric(gsub("\\D", "", hh_id))),
      p_sort_num  = suppressWarnings(as.numeric(person_id))
    ) %>%
    # Sort by the extracted numbers first.
    # (If no number exists, it falls back to alphabetical sorting)
    dplyr::arrange(hh_sort_num, hh_id, p_sort_num, person_id) %>%
    # Clean up the temporary sorting columns
    dplyr::select(-hh_sort_num, -p_sort_num)


  df_model_full <- df_meta_all %>%
    dplyr::left_join(df_imputed, by = "participantid")

  if(!"date_resolved_final" %in% names(df_model_full)) {
    df_model_full$date_resolved_final <- as.Date(NA)
    df_model_full$episode_id <- NA
    df_model_full$date_infection <- as.Date(NA)
    df_model_full$date_infectious_end <- as.Date(NA)
  }

  # df_model_full <- df_model_full %>%
  #   dplyr::group_by(participantid) %>%
  #   dplyr::mutate(
  #     prev_resolved_final = dplyr::lag(date_resolved_final, default = study_start_date - 1),
  #     start_risk_date = dplyr::case_when(
  #       is.na(episode_id) | episode_id == 1 ~ study_start_date,
  #       TRUE ~ prev_resolved_final + 1
  #     )
  #   ) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::mutate(
  #     start_risk = as.integer(start_risk_date - study_start_date) + 1,
  #     i_idx = dplyr::row_number()
  #   )
  # =========================================================
  # STRICT SORTING & RISK WINDOW CALCULATION
  # =========================================================

  df_model_full <- df_model_full %>%
    # 1. Extract pure numbers for strict sorting (HH10 comes after HH2)
    dplyr::mutate(
      hh_sort_num = suppressWarnings(as.numeric(gsub("\\D", "", hh_id))),
      p_sort_num  = suppressWarnings(as.numeric(gsub("\\D", "", person_id)))
    ) %>%
    # 2. Sort explicitly: Household -> Person -> Infection Date
    dplyr::arrange(hh_sort_num, hh_id, p_sort_num, person_id, date_infection) %>%

    # 3. Calculate Risk Windows
    dplyr::group_by(participantid) %>%
    dplyr::mutate(
      prev_resolved_final = dplyr::lag(date_resolved_final, default = study_start_date - 1),
      start_risk_date = dplyr::case_when(
        is.na(episode_id) | episode_id == 1 ~ study_start_date,
        TRUE ~ prev_resolved_final + 1
      )
    ) %>%
    dplyr::ungroup() %>%

    # 4. Re-affirm the sort (just in case group_by changed it)
    dplyr::arrange(hh_sort_num, hh_id, p_sort_num, person_id, start_risk_date) %>%

    # 5. GENERATE THE MASTER INDEX (i_idx)
    dplyr::mutate(
      start_risk = as.integer(start_risk_date - study_start_date) + 1,
      i_idx = dplyr::row_number() # <-- THIS IS THE ANCHOR FOR EVERY MATRIX
    ) %>%
    # Clean up
    dplyr::select(-hh_sort_num, -p_sort_num)



  N <- nrow(df_model_full)
  unique_hh <- unique(df_model_full$hh_id)
  H <- length(unique_hh)

  df_model_full$hh_id_int <- as.integer(factor(df_model_full$hh_id, levels = unique_hh))
  df_model_full$p_id_int  <- as.integer(as.factor(df_model_full$participantid))

  # --- COVARIATE MATRICES ---
  if (!is.null(covariates_susceptibility) && length(covariates_susceptibility) > 0) {
    X_susc <- as.matrix(df_model_full[, covariates_susceptibility, drop = FALSE])
    K_susc <- ncol(X_susc)
  } else {
    X_susc <- matrix(0, nrow = N, ncol = 0)
    K_susc <- 0
  }

  if (!is.null(covariates_infectivity) && length(covariates_infectivity) > 0) {
    X_inf <- as.matrix(df_model_full[, covariates_infectivity, drop = FALSE])
    K_inf <- ncol(X_inf)
  } else {
    X_inf <- matrix(0, nrow = N, ncol = 0)
    K_inf <- 0
  }


  I <- matrix(0L, N, T_max)
  Y <- matrix(0L, N, T_max)
  V <- matrix(default_val, N, T_max)

  # Empirical Helper
  calc_curve_val <- function(t, p, is_ct) {
    if(is_ct) {
      return(ifelse(t <= p$t_peak, p$Cpeak + p$r * (p$t_peak - t), p$Cpeak + p$d * (t - p$t_peak)))
    } else {
      val = 2 * 10^p$v_p / (exp(-p$lambda_g * (t - p$t_p)) + exp(p$lambda_d * (t - p$t_p)))
      return(log10(val))
    }
  }

  for (i in 1:N) {
    row <- df_model_full[i, ]

    if (!is.na(row$episode_id)) {
      idx_inf <- as.integer(row$date_infection - study_start_date) + 1
      idx_end <- as.integer(row$date_infectious_end - study_start_date) + 1

      if (idx_inf >= 1 && idx_inf <= T_max) I[i, idx_inf] <- 1L
      y_start <- max(1, idx_inf); y_end <- min(T_max, idx_end)
      if (y_end >= y_start) Y[i, y_start:y_end] <- 1L

      # MECHANISTIC IMPUTATION LOGIC
      p_role <- row$role
      days_seq <- seq(y_start, y_end)

      if(length(days_seq) > 0) {
        t_vals <- days_seq - idx_inf
        imputed_vals <- NULL

        # A. Ct Data
        if (is_ct_data) {
          p_params <- imputation_params[[p_role]]
          if(is.null(p_params)) p_params <- list(Cpeak=33, r=1.5, d=1.2, t_peak=5)
          imputed_vals <- sapply(t_vals, function(t) calc_curve_val(t, p_params, TRUE))

          # B. Viral Load Data
        } else {

          if (model_type == "ODE") {
            # 1. Prepare Parameters (Fix logic here: Check passed params first)
            current_ode_params <- ODE_params_list
            if(is.null(current_ode_params)) {
              # --- CORRECTED DEFAULTS TO MATCH YOUR SIMULATION ---
              base_pars <- list(
                beta = 5e-7, phi_R = 0.2, rho_R = 0.34, k1 = 4.0, delta0 = 0.5, delta1 = 3.04,
                p_I = 200, theta_F = 1e5, delta_V = 10, delta_F = 0.4,
                a = 1.0, b0 = 1.0, sigma_b = 0, tau_adaptive = 7.5
              )
              current_ode_params <- list(adult=base_pars, infant=base_pars, toddler=base_pars, elderly=base_pars)
            }

            p_vec <- unlist(current_ode_params[[p_role]])
            max_t <- max(t_vals)
            traj <- solve_ode_trajectory(duration_days = max_t + 2, params = p_vec)

            traj_indices <- t_vals + 1
            traj_indices <- traj_indices[traj_indices <= length(traj)]

            if(length(traj_indices) > 0) {
              imputed_vals <- traj[traj_indices]
            }

          } else {
            # --- EMPIRICAL PATH ---
            p_params <- imputation_params[[p_role]]
            if(is.null(p_params)) p_params <- list(v_p=5, t_p=4, lambda_g=2.8, lambda_d=1.0)
            imputed_vals <- sapply(t_vals, function(t) calc_curve_val(t, p_params, FALSE))
          }
        }

        # Fill Matrix
        if(!is.null(imputed_vals) && length(imputed_vals) == length(days_seq)) {
          V[i, days_seq] <- imputed_vals
        }

        # 3. OVERWRITE WITH OBSERVED DATA
        obs <- df_clean %>% dplyr::filter(participantid == row$participantid, !is.na(ct_value))
        for(k in seq_len(nrow(obs))) {
          d_idx <- as.integer(obs$date[k] - study_start_date) + 1
          if(d_idx >= 1 && d_idx <= T_max) V[i, d_idx] <- obs$ct_value[k]
        }
      }
    }
  }

  # --- GENERATE CONTACT MATRICES FOR STAN ---
  # Stan needs a flattened/sparse representation

  # Helper to generate matrix for a set of roles
  generate_contact_matrix_from_roles <- function(current_roles, role_mixing_matrix) {
    n <- length(current_roles)
    mat <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          weight <- tryCatch(role_mixing_matrix[current_roles[i], current_roles[j]], error=function(e) 1.0)
          mat[i, j] <- weight
        }
      }
    }
    return(mat)
  }

  hh_ids <- unique(df_model_full$hh_id_int)
  contact_src <- integer()
  contact_tgt <- integer()
  contact_w   <- numeric()


  for(h in hh_ids) {
    # Get members of this household
    members <- df_model_full[df_model_full$hh_id_int == h, ]
    member_indices <- members$i_idx # Global index (1 to N)
    member_roles   <- members$role

    # NEW: Get actual physical person IDs to detect clones
    member_pids    <- members$participantid

    # Generate Matrix based on user input or default
    if(!is.null(role_mixing_matrix)) {
      mat <- generate_contact_matrix_from_roles(member_roles, role_mixing_matrix)
    } else {
      # Default: Homogeneous Mixing
      n_h <- length(member_indices)
      mat <- matrix(1, n_h, n_h); diag(mat) <- 0
    }

    # --- NEW: ZERO OUT SELF-CONTACT ACROSS EPISODES ---
    for(r in 1:nrow(mat)) {
      for(c in 1:ncol(mat)) {
        # If the target (row) and source (col) are the EXACT SAME physical person
        if (member_pids[r] == member_pids[c]) {
          mat[r, c] <- 0.0
        }
      }
    }
    # --------------------------------------------------

    # Flatten to Sparse Format
    for(r in 1:nrow(mat)) {
      for(c in 1:ncol(mat)) {
        if(mat[r,c] > 0) {
          contact_tgt <- c(contact_tgt, member_indices[r]) # Row = Target
          contact_src <- c(contact_src, member_indices[c]) # Col = Source
          contact_w   <- c(contact_w,   mat[r,c])
        }
      }
    }
  }

  # --- HOUSEHOLD STRUCTURE ---
  hh_members_list <- split(1:N, df_model_full$hh_id_int)
  hh_size_eps <- sapply(hh_members_list, length)
  hh_max_size <- max(hh_size_eps)
  hh_members <- matrix(0L, nrow = H, ncol = hh_max_size)
  for (h in 1:H) hh_members[h, 1:hh_size_eps[h]] <- hh_members_list[[h]]

  # hh_size_df <- df_clean %>%
  #   dplyr::group_by(familyidstars) %>%
  #   dplyr::summarise(n = dplyr::n_distinct(participantid))
  # hh_size_people <- hh_size_df$n[match(unique_hh, hh_size_df$familyidstars)]

  hh_size_df <- df_model_full %>%
    dplyr::group_by(hh_id_int) %>%
    dplyr::summarise(n = dplyr::n_distinct(participantid)) %>% # <-- FIXED HERE
    dplyr::arrange(hh_id_int)

  hh_size_people <- hh_size_df$n


  # --- SEASONALITY ---
  seasonal_forcing_mat <- matrix(1.0, nrow=T_max, ncol=4)
  if (!is.null(surveillance_df)) {
    if (!all(c("date", "cases") %in% names(surveillance_df))) stop("surveillance_df must have 'date' and 'cases'")
    surveillance_df$date <- as.Date(surveillance_df$date)
    target_dates <- seq(from = study_start_date, to = study_end_date, by = "day")
    interp_res <- stats::approx(x = surveillance_df$date, y = surveillance_df$cases, xout = target_dates, method = "linear", rule = 2)
    daily_vec <- interp_res$y
    if(max(daily_vec, na.rm=TRUE) > 0) daily_vec <- daily_vec / max(daily_vec, na.rm=TRUE)
    daily_vec[is.na(daily_vec)] <- 0
    for(k in 1:4) seasonal_forcing_mat[,k] <- daily_vec
  } else if(!is.null(seasonal_forcing_list)) {
    # Assumes valid list
  }

  # =========================================================
  # 5. RETURN LIST
  # =========================================================

  list(
    N = N, T = T_max, H = H, R = length(role_levels), delta = delta,
    hh_id = df_model_full$hh_id_int,
    role_id = match(df_model_full$role, role_levels),
    I = I, Y = Y, V = V,
    start_risk = df_model_full$start_risk,

    # --- CONTACT DATA ---
    n_contacts  = length(contact_w),
    contact_src = contact_src,
    contact_tgt = contact_tgt,
    contact_w   = contact_w,
    # --------------------

    p_id = df_model_full$p_id_int,
    hh_size_people = as.integer(hh_size_people),
    hh_max_size = as.integer(hh_max_size),
    hh_members = hh_members,
    seasonal_forcing_mat = seasonal_forcing_mat,
    reference_phi = 1.0, reference_kappa = 1.0,
    use_vl_data = as.integer(use_vl_data),

    vl_type = as.integer(detected_vl_type),
    use_curve_logic = as.integer(use_curve_logic),

    K_susc = K_susc, X_susc = X_susc,
    K_inf  = K_inf,  X_inf  = X_inf,

    prior_beta1_type = p_beta1$type, prior_beta1_params = p_beta1$params,
    prior_beta2_type = p_beta2$type, prior_beta2_params = p_beta2$params,
    prior_alpha_type = p_alpha$type, prior_alpha_params = p_alpha$params,
    prior_cov_type   = p_cov$type,   prior_cov_params   = p_cov$params,

    prior_shape_type = p_shape$type, prior_shape_params = p_shape$params,
    prior_rate_type  = p_rate$type,  prior_rate_params  = p_rate$params,
    prior_vl_midpoint_type   = p_vl_midpoint$type,   prior_vl_midpoint_params = p_vl_midpoint$params,
    prior_vl_slope_type      = p_vl_slope$type,       prior_vl_slope_params    = p_vl_slope$params

  )
}
