
#' @import deSolve
NULL

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

simulate_viral_load_trajectory <- function(t, v_p, t_p, lambda_g, lambda_d) {
  vt = 2 * 10^v_p / (exp(-lambda_g * (t - t_p)) + exp(lambda_d * (t - t_p)))
  return(log10(vt))
}

draw_random_VL_params <- function(role, params) {
  p <- params[[role]]
  if (is.null(p)) stop("Unknown role in VL_params: ", role)
  list(v_p = p$v_p, t_p = p$t_p, lambda_g = p$lambda_g, lambda_d = p$lambda_d)
}

simulate_Ct_trajectory <- function(t, Cpeak, r, d, t_peak) {
  ct = ifelse(t <= t_peak, Cpeak + r * (t_peak - t), Cpeak + d * (t - t_peak))
  return(ct)
}

draw_random_Ct_params <- function(role, params) {
  p <- params[[role]]
  if (is.null(p)) stop("Unknown role in Ct_params: ", role)
  list(Cpeak = p$Cpeak, r = p$r, d = p$d, t_peak = p$t_peak)
}

#' ODE System for Within-Host Dynamics
#'
#' @param t Numeric. Current time point.
#' @param state Numeric vector. Current state variables.
#' @param parms List. Model parameters for the ODE system.
#'
#' @return A list containing the derivatives of the state variables.
ode_system_func <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # 1. Time-Dependent Parameters
    current_b <- if (t < tau_adaptive) b0 else b0 * exp(-sigma_b * (t - tau_adaptive))
    current_delta_I <- if (t < tau_adaptive) delta0 else delta1

    # 2. Enforce Non-Negativity
    T_cell <- max(0, T); R_cell <- max(0, R); I1 <- max(0, I1); I2 <- max(0, I2)
    V <- max(0, V); F_cyt <- max(0, F)

    # 3. Algebraic Quantities
    VI <- a * (V ^ current_b)
    inhibition <- 1 - (F_cyt / (F_cyt + theta_F))

    # 4. Differential Equations
    dT  <- -beta * VI * T_cell - phi_R * F_cyt * T_cell + rho_R * R_cell
    dR  <- phi_R * F_cyt * T_cell - rho_R * R_cell
    dI1 <- beta * VI * T_cell - k1 * I1
    dI2 <- k1 * I1 - current_delta_I * I2
    dV  <- p_I * I2 * inhibition - delta_V * V
    dF  <- I2 - delta_F * F_cyt

    list(c(dT, dR, dI1, dI2, dV, dF))
  })
}

solve_ode_trajectory <- function(duration_days, params) {
  # Initial Conditions (Assumed T=1e7, V=0 seed, I2=1/30 to start)
  y_init <- c(T = 1e7, R = 0.0, I1 = 0.0, I2 = 1/30, V = 0, F = 0.0)

  # Solve with high resolution
  times <- seq(0, duration_days, by = 0.1)

  out <- deSolve::ode(y = y_init, times = times, func = ode_system_func, parms = params)
  out_df <- as.data.frame(out)

  # Downsample to Daily Resolution (t = 0, 1, 2...) for simulator
  daily_df <- out_df[out_df$time %% 1 == 0, ]

  # Convert to Log10 for compatibility with main loop
  log_V <- log10(pmax(daily_df$V, 1))

  return(log_V)
}

draw_random_ODE_params <- function(role, params_list) {
  p <- params_list[[role]]
  if (is.null(p)) stop("Unknown role in ODE_params: ", role)
  unlist(p)
}

generate_household_roles <- function(country_profile) {
  base_profile <- list(
    prob_adults    = c(0.0, 0, 1),
    prob_infant    = 1.0,
    prob_siblings  = c(0.1, 0.7, 0.2),
    prob_elderly   = c(0.7, 0.15, 0.15)
  )

  if(is.null(country_profile)) country_profile <- list()
  profile <- utils::modifyList(base_profile, country_profile)

  roles <- c()
  n_adults <- sample(0:2, 1, prob = profile$prob_adults)
  if(n_adults > 0) roles <- c(roles, rep("adult", n_adults))

  has_infant <- rbinom(1, 1, profile$prob_infant)
  if(has_infant == 1) roles <- c(roles, "infant")

  n_toddlers <- sample(0:2, 1, prob = profile$prob_siblings)
  if(n_toddlers > 0) roles <- c(roles, rep("toddler", n_toddlers))

  n_elderly <- sample(0:2, 1, prob = profile$prob_elderly)
  if(n_elderly > 0) roles <- c(roles, rep("elderly", n_elderly))

  # ADD this at the end, before return:
  if (length(roles) == 0) {
    warning("Generated empty household - forcing at least one adult")
    roles <- c("adult")
  }

  return(roles)
}

#' Generate Household Contact Matrix from Role Rules
#' @param current_roles Vector of strings (e.g. c("adult", "adult", "infant"))
#' @param role_mixing_matrix 4x4 Matrix with row/col names: infant, toddler, adult, elderly
#' @return N x N matrix where N is length(current_roles)
generate_contact_matrix_from_roles <- function(current_roles, role_mixing_matrix) {
  n <- length(current_roles)
  mat <- matrix(0, nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) { # Skip diagonal
        role_i <- current_roles[i]
        role_j <- current_roles[j]

        weight <- tryCatch(
          role_mixing_matrix[role_i, role_j],
          error = function(e) 1.0 # Default to 1.0 if role not found
        )
        mat[i, j] <- weight
      }
    }
  }
  return(mat)
}

# ==============================================================================
# 2. CORE SIMULATION ENGINE
# ==============================================================================

simulate_one_household_comm <- function(hh_id,
                                        roles,
                                        alpha_comm_by_role,
                                        beta1, beta2, delta,
                                        phi_by_role, kappa_by_role,
                                        infectious_shape, infectious_scale,
                                        waning_shape, waning_scale,
                                        peak_day, width,
                                        max_days,
                                        test_weekly_before_detection,
                                        perfect_detection,
                                        contact_mat = NULL,
                                        verbose,
                                        seasonal_forcing_list,
                                        detect_threshold_log10,
                                        detect_threshold_Ct,
                                        surveillance_interval,
                                        test_daily,
                                        viral_testing,
                                        V_ref, V_rho,
                                        Ct_50, Ct_delta,
                                        VL_params_input,
                                        Ct_params_input,
                                        model_type = "empirical",
                                        ODE_params_input = NULL,
                                        susc_modifiers_vec = NULL,
                                        inf_modifiers_vec = NULL,
                                        covariate_data = NULL,
                                        max_infections = Inf) {

  n <- length(roles)

  # Guard against empty households
  if (n == 0) {
    return(list(
      hh_df = data.frame(
        hh_id = character(0), person_id = integer(0), role = character(0),
        infection_time = as.Date(character(0)), infectious_end = as.Date(character(0)),
        resolved_time = as.Date(character(0)), stringsAsFactors = FALSE
      ),
      diagnostic_df = data.frame(
        hh_id = character(0), person_id = integer(0), role = character(0),
        day_index = integer(0), pcr_sample = numeric(0), test_result = integer(0),
        episode_id = integer(0), stringsAsFactors = FALSE
      )
    ))
  }


  infection_history       <- vector("list", n)
  infectious_end_history <- vector("list", n)
  immunity_end_history   <- vector("list", n)

  infection_counts  <- integer(n)
  current_status    <- integer(n) # 0=S, 1=I, 2=R
  time_next_state   <- rep(NA_integer_, n)

  current_vl_traj   <- vector("list", n)
  current_inf_start <- rep(NA_integer_, n)
  detection_time    <- rep(NA_integer_, n)

  # Defaults
  if(is.null(susc_modifiers_vec)) susc_modifiers_vec <- rep(1.0, n)
  if(is.null(inf_modifiers_vec))  inf_modifiers_vec  <- rep(1.0, n)

  scaling_n <- (1.0 / max(n, 1))^delta
  phi_vec   <- phi_by_role[roles]
  kappa_vec <- kappa_by_role[roles]

  if (is.null(contact_mat)) {
    contact_mat <- matrix(1, n, n); diag(contact_mat) <- 0
  }

  household_detected <- FALSE

  # ==========================
  # MAIN TIME LOOP
  # ==========================
  for (t in 1:max_days) {

    # --- A. Update States ---
    for(i in seq_len(n)) {
      if(current_status[i] != 0 && !is.na(time_next_state[i]) && t >= time_next_state[i]) {
        st <- current_status[i]
        if(st == 1) { # I -> R
          current_status[i] <- 2
          infectious_end_history[[i]] <- c(infectious_end_history[[i]], t)
          dur <- pmax(1, ceiling(stats::rgamma(1, shape = waning_shape, scale = waning_scale)))
          time_next_state[i] <- t + dur
          immunity_end_history[[i]] <- c(immunity_end_history[[i]], t + dur)
          current_vl_traj[i] <- list(NULL)
          current_inf_start[i] <- NA
        } else if(st == 2) { # R -> S
          current_status[i] <- 0
          time_next_state[i] <- NA
        }
      }
    }

    is_scheduled_day <- ((t - 1) %% surveillance_interval == 0)
    test_today <- if(household_detected && test_daily) TRUE else is_scheduled_day

    # --- B. Calculate Force ---
    total_hh_force <- rep(0.0, n)
    infectors <- which(current_status == 1)

    if(length(infectors) > 0) {
      infectivity_values <- numeric(length(infectors))
      for(k in seq_along(infectors)) {
        idx <- infectors[k]
        rel_day <- t - current_inf_start[idx] + 1
        val <- NA_real_
        traj <- current_vl_traj[[idx]]
        if(!is.null(traj) && rel_day >= 1 && rel_day <= length(traj)) val <- traj[rel_day]

        term1 <- beta1 * 1.0
        term2 <- 0.0
        if(!is.na(val)) {
          if(viral_testing == "viral load") {
            val_clean <- max(0, val); term2 <- beta2 * (val_clean / V_ref)^V_rho
          } else {
            val_clean <- ifelse(val > 45, 45, val); term2 <- beta2 * 1.0 / (1.0 + exp((val_clean - Ct_50) / Ct_delta))
          }
        }
        infectivity_values[k] <- scaling_n * kappa_vec[idx] * inf_modifiers_vec[idx] * (term1 + term2)
      }
      for(k in seq_along(infectors)) {
        src <- infectors[k]; force <- infectivity_values[k]
        contacts <- which(contact_mat[, src] > 0)
        total_hh_force[contacts] <- total_hh_force[contacts] + (force * contact_mat[contacts, src])
      }
    }

    # --- C. New Infections ---
    targets <- which(current_status == 0 & infection_counts < max_infections)
    if(length(targets) > 0) {
      alpha_comm_val <- numeric(length(targets))
      for(k in seq_along(targets)) {
        tgt <- targets[k]; role_name <- roles[tgt]
        season_val <- seasonal_forcing_list[[role_name]][t]
        alpha_comm_val[k] <- alpha_comm_by_role * season_val
      }

      lambda_vec <- (phi_vec[targets] * susc_modifiers_vec[targets]) * (alpha_comm_val + total_hh_force[targets])
      lambda_vec <- pmin(lambda_vec, 1e6)
      prob_inf <- 1.0 - exp(-lambda_vec)

      is_infected <- stats::runif(length(targets)) < prob_inf
      infected_indices <- targets[is_infected]

      for(j in infected_indices) {
        infection_history[[j]] <- c(infection_history[[j]], t)
        infection_counts[j] <- infection_counts[j] + 1
        current_status[j] <- 1

        dur_I <- pmax(1, ceiling(stats::rgamma(1, shape = infectious_shape, scale = infectious_scale)))
        time_next_state[j] <- t + dur_I
        current_inf_start[j] <- t

        t_seq <- 0:dur_I

        if (viral_testing == "viral load") {
          if (model_type == "ODE") {
            p_vec <- draw_random_ODE_params(roles[j], ODE_params_input)
            traj <- solve_ode_trajectory(duration_days = dur_I + 2, params = p_vec)
            if(length(traj) > dur_I) traj <- traj[1:dur_I]
          } else {
            p <- draw_random_VL_params(roles[j], VL_params_input)
            traj <- simulate_viral_load_trajectory(t_seq, p$v_p, p$t_p, p$lambda_g, p$lambda_d)
          }
        } else {
          p <- draw_random_Ct_params(roles[j], Ct_params_input)
          traj <- simulate_Ct_trajectory(t_seq, p$Cpeak, p$r, p$d, p$t_peak)
        }
        current_vl_traj[[j]] <- traj
      }
    }

    # --- D. Testing ---
    if (test_today) {
      shedders <- which(current_status == 1)
      for(i in shedders) {
        rel_day <- t - current_inf_start[i] + 1
        traj <- current_vl_traj[[i]]
        if(!is.null(traj) && rel_day >= 1 && rel_day <= length(traj)) {
          val <- traj[rel_day]
          is_pos <- FALSE
          if(viral_testing == "viral load") { if(val >= detect_threshold_log10) is_pos <- TRUE }
          else { if(val <= detect_threshold_Ct) is_pos <- TRUE }
          if(is_pos && is.na(detection_time[i])) {
            detection_time[i] <- t; household_detected <- TRUE
          }
        }
      }
    }
  }

  # ==============================================================================
  # 3. EXPORT RESULTS
  # ==============================================================================

  hh_rows <- list()
  for(i in seq_len(n)) {
    hist <- infection_history[[i]]
    base_row <- data.frame(hh_id=hh_id, person_id=i, role=roles[i], stringsAsFactors=F)

    if(!is.null(covariate_data)) {
      covs_clean <- covariate_data[i, !names(covariate_data) %in% "person_id", drop=FALSE]
      base_row <- cbind(base_row, covs_clean)
    }

    if(is.null(hist)) {
      base_row$infection_time <- NA; base_row$infectious_end <- NA; base_row$resolved_time <- NA
      hh_rows[[length(hh_rows)+1]] <- base_row
    } else {
      for(k in seq_along(hist)) {
        row <- base_row
        row$infection_time <- hist[k]
        row$infectious_end <- if(k <= length(infectious_end_history[[i]])) infectious_end_history[[i]][k] else NA
        row$resolved_time  <- if(k <= length(immunity_end_history[[i]])) immunity_end_history[[i]][k] else NA
        hh_rows[[length(hh_rows)+1]] <- row
      }
    }
  }
  hh_df <- do.call(rbind, hh_rows)

  # --- Diagnostic DF ---
  test_days <- seq(1, max_days, by=surveillance_interval)
  if(any(!is.na(detection_time)) && test_daily) {
    first_det <- min(detection_time, na.rm=T)
    if(!is.infinite(first_det)) test_days <- unique(sort(c(test_days, seq(first_det, max_days))))
  }
  n_tests <- length(test_days)
  res_matrix <- matrix(if(viral_testing=="viral load") 0 else 45, nrow=n, ncol=n_tests)
  ep_matrix <- matrix(0, nrow=n, ncol=n_tests)

  for(i in seq_len(n)) {
    hist <- infection_history[[i]]
    if(!is.null(hist)) {
      for(k in seq_along(hist)) {
        inf_t   <- hist[k]
        inf_end <- if(k <= length(infectious_end_history[[i]])) infectious_end_history[[i]][k] else (inf_t + 10)
        relevant_indices <- which(test_days >= inf_t & test_days <= inf_end)

        if(length(relevant_indices) > 0) {
          d_vals <- test_days[relevant_indices]; rel_vals <- d_vals - inf_t + 1

          if(viral_testing == "viral load") {
            if(model_type == "ODE") {
              p_vec <- draw_random_ODE_params(roles[i], ODE_params_input)
              max_rel <- max(rel_vals)
              traj <- solve_ode_trajectory(max_rel + 2, p_vec)
              valid_rel <- rel_vals[rel_vals <= length(traj)]
              valid_idx <- relevant_indices[rel_vals <= length(traj)]
              v_vals <- traj[valid_rel]
              res_matrix[i, valid_idx] <- pmax(res_matrix[i, valid_idx], v_vals)
            } else {
              p <- draw_random_VL_params(roles[i], VL_params_input)
              v_vals <- simulate_viral_load_trajectory(rel_vals, p$v_p, p$t_p, p$lambda_g, p$lambda_d)
              res_matrix[i, relevant_indices] <- pmax(res_matrix[i, relevant_indices], v_vals)
            }
          } else {
            p <- draw_random_Ct_params(roles[i], Ct_params_input)
            ct_vals <- simulate_Ct_trajectory(rel_vals, p$Cpeak, p$r, p$d, p$t_peak)
            res_matrix[i, relevant_indices] <- pmin(res_matrix[i, relevant_indices], ct_vals)
          }
        }
        if(length(relevant_indices) > 0) ep_matrix[i, relevant_indices] <- k
      }
    }
  }

  diag_df_list <- vector("list", n)
  for(i in seq_len(n)) {
    vals <- res_matrix[i, ]
    eps  <- ep_matrix[i, ]
    results <- if(viral_testing=="viral load") as.integer(vals >= detect_threshold_log10) else as.integer(vals <= detect_threshold_Ct)
    diag_df_list[[i]] <- data.frame(
      hh_id=hh_id, person_id=i, role=roles[i],
      day_index=test_days, pcr_sample=vals, test_result=results,
      episode_id=eps,
      stringsAsFactors=FALSE
    )
  }
  diagnostic_df <- do.call(rbind, diag_df_list)

  list(hh_df = hh_df, diagnostic_df = diagnostic_df)
}

# ==============================================================================
# 3. EXPORTED WRAPPER
# ==============================================================================

#' Simulate Household Transmission
#'
#' Simulates infection dynamics across multiple households with community forcing.
#'
#' @param n_households Integer. Number of households to simulate. Defaults to 50.
#' @param surveillance_df Dataframe with columns 'date' and 'cases' for community forcing. NULL for none.
#' @param start_date Character. Simulation start date ("YYYY-MM-DD"). Defaults to "2024-07-01".
#' @param end_date Character. Simulation end date ("YYYY-MM-DD"). Defaults to "2025-06-30".
#' @param alpha_comm_by_role Numeric or named vector. Community acquisition rate by role. Defaults to 5e-4.
#' @param beta1 Numeric. Within-household transmission rate (pathway 1). Defaults to 8e-3.
#' @param beta2 Numeric. Within-household transmission rate (pathway 2). Defaults to 8e-3.
#' @param delta Numeric. Co-infection parameter. Defaults to 0.
#' @param phi_by_role Named numeric vector. Susceptibility multipliers by role.
#' @param kappa_by_role Named numeric vector. Infectivity multipliers by role.
#' @param infectious_shape Numeric. Shape parameter for the Gamma infectious period. Defaults to 3.
#' @param infectious_scale Numeric. Scale parameter for the Gamma infectious period. Defaults to 1.
#' @param waning_shape Numeric. Shape parameter for the Gamma immunity waning period. Defaults to 16.
#' @param waning_scale Numeric. Scale parameter for the Gamma immunity waning period. Defaults to 10.
#' @param peak_day Numeric. Day of peak infectiousness. Defaults to 1.
#' @param width Numeric. Width of the infectiousness peak. Defaults to 4.
#' @param verbose Logical. Print progress messages. Defaults to FALSE.
#' @param seasonal_forcing_list List of numeric vectors for seasonal forcing by role. NULL for none.
#' @param detect_threshold_log10 Numeric. Detection threshold on log10 viral load scale. Defaults to 1e-6.
#' @param detect_threshold_Ct Numeric. Detection threshold on Ct scale. Defaults to 99.
#' @param surveillance_interval Integer. Days between surveillance tests. Defaults to 1.
#' @param test_daily Logical. Whether to test daily. Defaults to FALSE.
#' @param viral_testing Character. Type of viral testing: "viral load" or "Ct". Defaults to "viral load".
#' @param V_ref Numeric. Reference viral load for transmission scaling. Defaults to 3.0.
#' @param V_rho Numeric. Exponent for viral load scaling. Defaults to 2.5.
#' @param Ct_50 Numeric. Ct value at 50 percent detection probability. Defaults to 40.
#' @param Ct_delta Numeric. Steepness of Ct detection curve. Defaults to 2.
#' @param VL_params_list List of viral load trajectory parameters by role. NULL for defaults.
#' @param Ct_params_list List of Ct trajectory parameters by role. NULL for defaults.
#' @param household_profile_list List defining household composition probabilities. NULL for defaults.
#' @param perfect_detection Logical. Whether detection is perfect. Defaults to TRUE.
#' @param contact_mat Matrix. Custom contact matrix between individuals. NULL for default.
#' @param role_mixing_matrix Matrix. Contact weights between roles. 4x4 Matrix where element (i,j) represents
#'   contact weight FROM role j TO role i. For asymmetric patterns,
#'   role_mixing_matrix(adult, infant) != role_mixing_matrix(infant, adult).
#'   Use dimnames(role_mixing_matrix) <- list(role_levels, role_levels). NULL for default.
#' @param model_type Character. Either "empirical" or "ODE". Defaults to "empirical".
#' @param ODE_params_list List of ODE parameters by role. NULL for defaults.
#' @param covariates_config List defining covariate configurations. NULL for none.
#' @param seed Integer. Random seed. NULL for no seed.
#' @param max_infections Numeric. Maximum infections per person. Defaults to Inf.
#'
#' @return A list with two elements: hh_df (household-level results) and diagnostic_df (testing results).
#' @export
simulate_multiple_households_comm <- function(n_households = 50,
                                              surveillance_df = NULL,
                                              start_date = "2024-07-01",
                                              end_date = "2025-06-30",
                                              alpha_comm_by_role = 5e-4,
                                              beta1 = 8e-3, beta2 = 8e-3, delta = 0,
                                              phi_by_role = c(adult = 1, infant = 4, toddler = 5, elderly = 1),
                                              kappa_by_role = c(adult = 1, infant = 1, toddler = 1.2, elderly = 1),
                                              infectious_shape = 3, infectious_scale = 1,
                                              waning_shape = 16, waning_scale = 10,
                                              peak_day = 1, width = 4,
                                              verbose = FALSE,
                                              seasonal_forcing_list = NULL,
                                              detect_threshold_log10 = 1e-6,
                                              detect_threshold_Ct = 99,
                                              surveillance_interval = 1,
                                              test_daily = FALSE,
                                              viral_testing = "viral load",
                                              V_ref = 3.0, V_rho = 2.5,
                                              Ct_50 = 40, Ct_delta = 2,
                                              VL_params_list = NULL,
                                              Ct_params_list = NULL,
                                              household_profile_list = NULL,
                                              perfect_detection = TRUE,

                                              # --- CONTACT ARGS ---
                                              contact_mat = NULL,
                                              role_mixing_matrix = NULL,
                                              # --------------------

                                              model_type = "empirical",
                                              ODE_params_list = NULL,
                                              covariates_config = NULL,
                                              seed = NULL,
                                              max_infections = Inf) {

  if (!is.null(seed)) {
    set.seed(seed)
  }
  d_start <- as.Date(start_date); d_end <- as.Date(end_date)
  max_days <- as.integer(d_end - d_start) + 1
  if (max_days <= 0) stop("end_date must be after start_date")

  final_forcing_list <- NULL
  if (!is.null(surveillance_df)) {
    target_dates <- seq(from = d_start, to = d_end, by = "day")
    interp_res <- stats::approx(x=as.Date(surveillance_df$date), y=surveillance_df$cases, xout=target_dates, rule=2)
    daily_vec <- interp_res$y; if(max(daily_vec, na.rm=T)>0) daily_vec <- daily_vec/max(daily_vec, na.rm=T); daily_vec[is.na(daily_vec)]<-0
    final_forcing_list <- list(adult=daily_vec, infant=daily_vec, toddler=daily_vec, elderly=daily_vec)
  } else if (!is.null(seasonal_forcing_list)) {
    final_forcing_list <- seasonal_forcing_list
  }
  if(is.null(final_forcing_list)) final_forcing_list <- list(adult=rep(1,max_days), infant=rep(1,max_days), toddler=rep(1,max_days), elderly=rep(1,max_days))

  if(is.null(household_profile_list)) household_profile_list <- list(prob_single_parent=0, prob_siblings=c(0.10, 0.50, .40), prob_elderly=c(0.9, 0.08, 0.02))
  if(is.null(VL_params_list)) VL_params_list <- list(adult=list(v_p=4.14, t_p=5.09, lambda_g=2.31, lambda_d=2.71), infant=list(v_p=5.84, t_p=4.09, lambda_g=2.82, lambda_d=1.01), toddler=list(v_p=5.84, t_p=4.09, lambda_g=2.82, lambda_d=1.01), elderly=list(v_p=2.95, t_p=5.1, lambda_g=3.15, lambda_d=0.87))
  if(is.null(Ct_params_list)) Ct_params_list <- list(infant=list(Cpeak=33.3, r=2.11, d=1.38, t_peak=5.06), toddler=list(Cpeak=34, r=1.26, d=1.27, t_peak=4.75), adult=list(Cpeak=33, r=1.49, d=1.22, t_peak=5.14), elderly=list(Cpeak=33, r=1.49, d=1.22, t_peak=5.14))

  if(model_type == "ODE" && is.null(ODE_params_list)) {
    base_pars <- list(
      beta = 5e-7, phi_R = 0.2, rho_R = 0.34, k1 = 4.0, delta0 = 0.5, delta1 = 3.04,
      p_I = 200, theta_F = 1e5, delta_V = 10, delta_F = 0.4,
      a = 1.0, b0 = 1.0, sigma_b = 0, tau_adaptive = 7.5
    )
    ODE_params_list <- list(adult=base_pars, infant=base_pars, toddler=base_pars, elderly=base_pars)
  }

  households <- vector("list", n_households)

  for (h in seq_len(n_households)) {
    roles <- generate_household_roles(household_profile_list)
    n_hh <- length(roles)

    # --- CHANGED: Logic to build contact matrix from roles ---
    hh_contact_mat <- NULL

    if (!is.null(role_mixing_matrix)) {
      hh_contact_mat <- generate_contact_matrix_from_roles(roles, role_mixing_matrix)
    } else if (!is.null(contact_mat)) {
      if(nrow(contact_mat) == n_hh) hh_contact_mat <- contact_mat
    }
    # -------------------------------------------------------

    hh_covariates <- data.frame(person_id = 1:n_hh)
    hh_susc_modifiers <- rep(1.0, n_hh)
    hh_inf_modifiers  <- rep(1.0, n_hh)

    if (!is.null(covariates_config)) {
      for (cov in covariates_config) {
        col_vec <- numeric(n_hh)
        for (i in seq_len(n_hh)) {
          prob <- cov$coverage[[roles[i]]]
          if (is.null(prob)) prob <- 0
          col_vec[i] <- rbinom(1, 1, prob)
        }
        hh_covariates[[cov$name]] <- col_vec
        eff_type <- tolower(cov$effect_on)
        multiplier <- (1.0 - (col_vec * cov$efficacy))
        if (eff_type %in% c("susceptibility", "both")) hh_susc_modifiers <- hh_susc_modifiers * multiplier
        if (eff_type %in% c("infectivity", "both")) hh_inf_modifiers <- hh_inf_modifiers * multiplier
      }
    }

    hh <- simulate_one_household_comm(hh_id = paste0("HH", h),
                                      roles = roles,
                                      alpha_comm_by_role = alpha_comm_by_role,
                                      beta1 = beta1, beta2 = beta2, delta = delta,
                                      phi_by_role = phi_by_role, kappa_by_role = kappa_by_role,
                                      infectious_shape = infectious_shape, infectious_scale = infectious_scale,
                                      waning_shape = waning_shape, waning_scale = waning_scale,
                                      peak_day = peak_day, width = width,
                                      max_days = max_days,
                                      verbose = verbose,
                                      seasonal_forcing_list = final_forcing_list,
                                      detect_threshold_log10 = detect_threshold_log10,
                                      detect_threshold_Ct = detect_threshold_Ct,
                                      surveillance_interval = surveillance_interval,
                                      test_daily = test_daily,
                                      viral_testing = viral_testing,
                                      V_ref = V_ref, V_rho = V_rho,
                                      Ct_50 = Ct_50, Ct_delta = Ct_delta,
                                      VL_params_input = VL_params_list,
                                      Ct_params_input = Ct_params_list,
                                      perfect_detection = perfect_detection,

                                      # --- CHANGED: Passing the generated matrix ---
                                      contact_mat = hh_contact_mat,
                                      # ---------------------------------------------

                                      model_type = model_type,
                                      ODE_params_input = ODE_params_list,
                                      test_weekly_before_detection = TRUE,
                                      susc_modifiers_vec = hh_susc_modifiers,
                                      inf_modifiers_vec  = hh_inf_modifiers,
                                      covariate_data = hh_covariates,
                                      max_infections = max_infections)
    households[[h]] <- hh
  }

  hh_list <- lapply(households, function(x) x$hh_df)
  hh_df <- dplyr::bind_rows(hh_list)

  diag_list <- lapply(households, function(x) x$diagnostic_df)
  diag_list <- diag_list[!sapply(diag_list, is.null)]
  diagnostic_df <- if(length(diag_list) > 0) do.call(rbind, diag_list) else data.frame()

  list(hh_df = hh_df, diagnostic_df = diagnostic_df)
}
