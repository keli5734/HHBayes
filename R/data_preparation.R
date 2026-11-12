# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Resolve the viral scale (log10 viral load vs Ct)
#'
#' @param df Data frame with a \code{pcr_sample} or \code{ct_value} column.
#' @param viral_scale One of "auto", "log10", "Ct". "auto" keeps the historical
#'   heuristic (values above 15 are taken to be Ct), which is fragile; set the
#'   scale explicitly whenever you know it.
#' @return list(is_ct, vl_type, default_val, scale)
#' @keywords internal
.resolve_viral_scale <- function(df, viral_scale = c("auto", "log10", "Ct"), verbose = TRUE) {
  viral_scale <- match.arg(viral_scale)
  col <- intersect(c("pcr_sample", "ct_value"), names(df))[1]
  max_obs <- if (is.na(col)) NA_real_ else suppressWarnings(max(df[[col]], na.rm = TRUE))
  looks_ct <- is.finite(max_obs) && max_obs > 15

  if (viral_scale == "auto") {
    is_ct <- looks_ct
    if (verbose) message(sprintf(
      "viral_scale = 'auto': max observed value %.2f => treated as %s. Set viral_scale explicitly to avoid this guess.",
      ifelse(is.finite(max_obs), max_obs, NA_real_), if (is_ct) "Ct" else "log10 viral load"))
  } else {
    is_ct <- (viral_scale == "Ct")
    if (is.finite(max_obs) && is_ct != looks_ct) {
      warning(sprintf("viral_scale = '%s' but the largest observed value is %.2f, which looks like %s data.",
                      viral_scale, max_obs, if (looks_ct) "Ct" else "log10 viral load"), call. = FALSE)
    }
  }
  list(is_ct = is_ct, vl_type = if (is_ct) 0L else 1L,
       default_val = if (is_ct) 45.0 else 0.0, scale = if (is_ct) "Ct" else "log10")
}

#' Standardize column names of the observation data frame
#'
#' Shared by \code{prepare_stan_data()} and \code{prepare_kinetics_data()}.
#' Keeps the original positivity flag (\code{test_result}) in \code{.test_pos}
#' before it is renamed, so downstream code can tell positives from negatives.
#' @keywords internal
.standardize_clean <- function(df_clean, study_start_date) {
  df <- as.data.frame(df_clean, stringsAsFactors = FALSE)
  if ("test_result" %in% names(df) && !".test_pos" %in% names(df)) {
    df$.test_pos <- as.integer(df$test_result) == 1L
  }
  if ("hh_id" %in% names(df) && !"familyidstars" %in% names(df)) names(df)[names(df) == "hh_id"] <- "familyidstars"
  if ("person_id" %in% names(df)) df$participantid <- paste(df$familyidstars, df$person_id, sep = "_")
  if ("role" %in% names(df) && !"role_name" %in% names(df)) names(df)[names(df) == "role"] <- "role_name"
  if ("pcr_sample" %in% names(df) && !"ct_value" %in% names(df)) names(df)[names(df) == "pcr_sample"] <- "ct_value"
  if ("test_result" %in% names(df) && !"is_in_episode" %in% names(df)) names(df)[names(df) == "test_result"] <- "is_in_episode"
  if ("day_index" %in% names(df) && !"date" %in% names(df)) {
    df$date <- as.Date(study_start_date) + (df$day_index - 1)
  }
  if ("date" %in% names(df)) df$date <- as.Date(df$date)
  df
}

# =============================================================================
# prepare_stan_data()
# =============================================================================

#' Prepare Data for Stan Model
#'
#' Builds the data list for the household transmission model. Infection and
#' infectious-period timing can come from two sources:
#' \itemize{
#'   \item \strong{Random imputation} (default, \code{episode_timelines = NULL}): the infection
#'     day is drawn uniformly between the last negative and the first positive sample and the
#'     end of infectiousness between the last positive and the next negative.
#'   \item \strong{Stage 1 draw} (\code{episode_timelines} supplied): timing comes from a posterior
#'     draw of the viral kinetics model (see \code{prepare_kinetics_data()},
#'     \code{fit_viral_kinetics()} and \code{kinetics_draw_to_timelines()}), so the negatives
#'     before the first and after the last positive inform onset directly. Fit the transmission
#'     model once per draw and pool the posterior draws to propagate Stage 1 uncertainty.
#' }
#' In both cases \code{recovery_params} still adds the post-infectious immunity tail that
#' determines when a person is at risk again.
#'
#' @section At-risk windows:
#' With \code{risk_window = "observed"} (default) a person is at risk from their first sample
#' (or the day after the previous episode's resolution) until the infection day, or until their
#' \emph{last} sample if never infected. Days outside observation are not counted as
#' susceptible person-time. Infections dated before observation began (prevalent infections)
#' still transmit but contribute no at-risk time. Reinfections cannot precede the resolution of
#' the previous episode. \code{risk_window = "study"} reproduces the old behaviour (everyone at
#' risk over the whole study window).
#'
#' @param df_clean Dataframe with observation data (must contain 'episode_id' from simulation or data).
#' @param surveillance_df Dataframe with columns 'date' and 'cases'.
#' @param role_levels Character vector. Role categories to use. Defaults to \code{c("adult", "infant", "toddler", "elderly")}.
#' @param study_start_date Date. Start date of the study period. Defaults to 2024-07-01.
#' @param study_end_date Date. End date of the study period. Defaults to 2025-07-01.
#' @param seasonal_forcing_list List of numeric vectors (one per role, length = study days) for seasonal forcing. NULL for no forcing.
#' @param use_vl_data Logical. Whether to include viral load data. Defaults to TRUE.
#' @param use_curve_logic Logical. Whether to use the Gamma generation interval curve
#'   as a fallback for infectiousness when \code{use_vl_data = FALSE}.
#'   When \code{use_vl_data = TRUE} this argument is ignored. Defaults to FALSE.
#' @param covariates_susceptibility Vector of column names to use as covariates for susceptibility.
#' @param covariates_infectivity Vector of column names to use as covariates for infectivity.
#' @param priors List of flexible priors (dist, params).
#' @param recovery_params List of Gamma parameters (shape, scale) for the immunity tail by role.
#' @param imputation_params List of parameters for mechanistic viral curves by role (used when no
#'   Stage 1 curve is available): Ct needs \code{Cpeak, r, d, t_peak}; log10 viral load needs
#'   \code{v_p, t_p, lambda_g, lambda_d}.
#' @param model_type String, either "empirical" or "ODE".
#' @param ODE_params_list List of ODE parameters (beta, delta, etc.) by role.
#' @param delta Numeric. Household size scaling parameter.
#'   When delta > 0, transmission rates are scaled by (1/max(household_size, 1))^delta.
#'   Should match the value used in simulation. Defaults to 0 (no scaling).
#' @param role_mixing_matrix 4x4 Matrix defining contact weights between roles.
#' @param viral_scale "auto", "log10" or "Ct". The scale of the viral measurements. "auto"
#'   guesses from the data (largest value above 15 means Ct); set it explicitly when known.
#' @param episode_timelines Optional data frame with one row per episode:
#'   \code{participantid}, \code{episode_id}, \code{date_infection}, and optionally
#'   \code{date_infectious_start}, \code{date_infectious_end} (see
#'   \code{kinetics_draw_to_timelines()}). Episodes missing from the table fall back to the
#'   random imputation (with a message).
#' @param episode_curves Optional long data frame (\code{participantid}, \code{episode_id},
#'   \code{date}, \code{value}) with the expected viral value on each infectious day, from the
#'   same Stage 1 draw. Observed values still overwrite the curve on tested days.
#' @param risk_window "observed" (default) or "study"; see the section on at-risk windows.
#' @param seed Integer. Random seed for reproducibility. Defaults to 123.
#'
#' @return A named list formatted for input to the Stan model. Besides the Stan data it carries
#'   \code{inert_pars} and a \code{units} data frame (one row per unit with its dates and
#'   at-risk window) for diagnostics; both are dropped by \code{fit_household_model()}.
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

                              # --- NEW: viral scale, Stage 1 timing, at-risk windows ---
                              viral_scale = c("auto", "log10", "Ct"),
                              episode_timelines = NULL,
                              episode_curves = NULL,
                              risk_window = c("observed", "study"),

                              seed = 123) {

  if(!is.null(seed)) set.seed(seed)
  risk_window <- match.arg(risk_window)
  study_start_date <- as.Date(study_start_date)
  study_end_date   <- as.Date(study_end_date)

  T_max <- as.integer(study_end_date - study_start_date) + 1
  if (T_max <= 0) stop("study_end_date must be after study_start_date")

  sc <- .resolve_viral_scale(df_clean, viral_scale)
  is_ct_data       <- sc$is_ct
  detected_vl_type <- sc$vl_type
  default_val      <- sc$default_val

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
  p_phi_role   <- parse_prior(priors$phi_role,   1, c(0, 1))
  p_kappa_role <- parse_prior(priors$kappa_role, 1, c(0, 1))

  # =========================================================
  # 1b. VALIDATE PRIOR SCALE (see manuscript prior table)
  # =========================================================
  # phi_role, kappa_role, beta1, beta2, alpha are sampled on the LOG scale.
  #   Normal(mu, sigma) here  ==  LogNormal(mu, sigma) on the biological quantity.
  #   Uniform(a, b) here      ==  log-uniform on the biological quantity.
  # A 'lognormal' prior on a log-scale parameter is one-sided and NOT a lognormal
  # on the biological scale, so it is disallowed.
  log_scale_priors <- c("phi_role", "kappa_role", "beta1", "beta2", "alpha")
  for (nm in log_scale_priors) {
    d <- priors[[nm]]$dist
    if (!is.null(d) && tolower(d) == "lognormal") {
      stop(sprintf(paste0(
        "Prior for '%s' is applied on the LOG scale, so 'lognormal' is one-sided and\n",
        "  does not give a lognormal on the biological scale.\n",
        "  * For LogNormal(mu, sigma) on the biological %s, use dist='normal', params=c(mu, sigma).\n",
        "  * Otherwise use dist='normal' or dist='uniform' (interpreted on the log scale)."),
        nm, nm))
    }
  }

  # Covariate-effect priors: the Stan model implements only normal/uniform;
  # anything else is silently treated as uniform, so reject it explicitly.
  d_cov <- priors$covariates$dist
  if (!is.null(d_cov) && !tolower(d_cov) %in% c("normal", "uniform")) {
    stop("Covariate-effect priors support only 'normal' or 'uniform' (log-linear coefficients).")
  }

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

  # Defaults are lists with a `dist` string, so convert them with parse_prior() too
  # (the integer codes are what the Stan model expects). A prior supplied by the user
  # takes precedence.
  p_vl_midpoint <- parse_prior(if (is.null(priors$vl_midpoint)) default_vl_midpoint_prior else priors$vl_midpoint,
                               1, default_vl_midpoint_prior$params)
  p_vl_slope    <- parse_prior(if (is.null(priors$vl_slope)) default_vl_slope_prior else priors$vl_slope,
                               1, default_vl_slope_prior$params)

  # =========================================================
  # 2. STANDARDIZE COLUMN NAMES
  # =========================================================
  df_clean <- .standardize_clean(df_clean, study_start_date)

  # Per-person observation span (first / last sample) -> at-risk windows
  person_span <- df_clean %>%
    dplyr::filter(!is.na(date)) %>%
    dplyr::group_by(participantid) %>%
    dplyr::summarise(first_obs_date = min(date), last_obs_date = max(date), .groups = "drop") %>%
    as.data.frame(stringsAsFactors = FALSE)
  if (risk_window == "study") {
    person_span$first_obs_date <- study_start_date
    person_span$last_obs_date  <- study_end_date
  }

  # =========================================================
  # 3. TIMING STEP: Stage 1 draw OR random interval imputation, + recovery tail
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

  # --- optional Stage 1 timelines / curves ---
  tl_lookup <- NULL
  if (!is.null(episode_timelines)) {
    tl <- as.data.frame(episode_timelines, stringsAsFactors = FALSE)
    miss <- setdiff(c("participantid", "episode_id", "date_infection"), names(tl))
    if (length(miss) > 0) stop("episode_timelines is missing column(s): ", paste(miss, collapse = ", "))
    tl$date_infection <- as.Date(tl$date_infection)
    if (!"date_infectious_start" %in% names(tl)) tl$date_infectious_start <- tl$date_infection
    if (!"date_infectious_end"   %in% names(tl)) tl$date_infectious_end   <- as.Date(NA)
    tl$date_infectious_start <- as.Date(tl$date_infectious_start)
    tl$date_infectious_end   <- as.Date(tl$date_infectious_end)
    tl_lookup <- split(tl, paste(tl$participantid, tl$episode_id, sep = "||"))
  }
  curve_lookup <- NULL
  if (!is.null(episode_curves)) {
    cv <- as.data.frame(episode_curves, stringsAsFactors = FALSE)
    miss <- setdiff(c("participantid", "episode_id", "date", "value"), names(cv))
    if (length(miss) > 0) stop("episode_curves is missing column(s): ", paste(miss, collapse = ", "))
    cv$date <- as.Date(cv$date)
    curve_lookup <- split(cv, paste(cv$participantid, cv$episode_id, sep = "||"))
  }

  df_by_pid    <- split(df_clean, df_clean$participantid)
  participants <- unique(df_clean$participantid)
  imputed_episodes <- list()
  n_fallback <- 0L

  for (pid in participants) {
    p_data <- df_by_pid[[pid]]
    p_data <- p_data[order(p_data$date), , drop = FALSE]
    p_role <- p_data$role_name[1]
    if(is.na(p_role) || !p_role %in% role_levels) p_role <- "adult"
    cur_resolve <- params_resolve[[p_role]]
    if(is.null(cur_resolve)) cur_resolve <- list(shape=2, scale=3)

    episode_ids <- if ("episode_id" %in% names(p_data)) unique(p_data$episode_id[which(p_data$episode_id > 0)]) else integer(0)
    prev_resolved_final <- as.Date(NA)

    for (k_ep in seq_along(episode_ids)) {
      eid <- episode_ids[k_ep]
      ep_rows <- p_data[which(p_data$episode_id == eid), , drop = FALSE]
      if (nrow(ep_rows) == 0) next

      first_pos_date <- min(ep_rows$date)
      last_pos_date  <- max(ep_rows$date)

      tl_row <- if (!is.null(tl_lookup)) tl_lookup[[paste(pid, eid, sep = "||")]] else NULL

      if (!is.null(tl_row)) {
        # ---- Stage 1 timing ----
        T_inf_new      <- tl_row$date_infection[1]
        T_inf_start    <- tl_row$date_infectious_start[1]
        T_resolved_new <- tl_row$date_infectious_end[1]
        if (is.na(T_inf_start))    T_inf_start    <- T_inf_new
        if (is.na(T_resolved_new)) T_resolved_new <- last_pos_date
      } else {
        # ---- random interval imputation (original method) ----
        if (!is.null(tl_lookup)) n_fallback <- n_fallback + 1L
        prev_neg_dates <- p_data$date[which(p_data$date < first_pos_date & p_data$episode_id == 0)]
        lower_bound <- if(length(prev_neg_dates) > 0) max(prev_neg_dates) else (first_pos_date - 7)

        next_neg_dates <- p_data$date[which(p_data$date > last_pos_date & p_data$episode_id == 0)]
        upper_bound <- if(length(next_neg_dates) > 0) min(next_neg_dates) else (last_pos_date + 7)

        possible_starts <- seq(from = lower_bound + 1, to = first_pos_date, by = "day")
        T_inf_new <- if(length(possible_starts) > 1) sample(possible_starts, 1) else first_pos_date

        possible_ends <- seq(from = last_pos_date, to = upper_bound - 1, by = "day")
        T_resolved_new <- if(length(possible_ends) > 1) sample(possible_ends, 1) else last_pos_date

        T_inf_start <- T_inf_new
      }

      # A reinfection cannot precede recovery from the previous episode. (An infection
      # dated before the unit's at-risk window opens would otherwise be silently dropped
      # by the Stan model.)
      if (!is.na(prev_resolved_final) && T_inf_new <= prev_resolved_final) {
        T_inf_new <- prev_resolved_final + 1
      }
      if (T_inf_new > first_pos_date) T_inf_new <- first_pos_date
      if (T_inf_start < T_inf_new) T_inf_start <- T_inf_new
      if (T_resolved_new < T_inf_start) T_resolved_new <- T_inf_start

      # Post-infectious immunity tail (recovery_params): the person is not at risk again
      # until T_resolved_final.
      resolve_draw <- max(1, round(stats::rgamma(1, shape = cur_resolve$shape, scale = cur_resolve$scale)))
      T_resolved_final <- T_resolved_new + resolve_draw

      next_eid <- if (k_ep < length(episode_ids)) episode_ids[k_ep + 1] else NA
      if (!is.na(next_eid)) {
        next_ep_raw_start <- min(p_data$date[which(p_data$episode_id == next_eid)])
        if (T_resolved_final >= next_ep_raw_start) T_resolved_final <- next_ep_raw_start - 1
      }
      prev_resolved_final <- T_resolved_final

      imputed_episodes[[length(imputed_episodes) + 1]] <- data.frame(
        participantid = pid, episode_id = eid,
        date_infection = T_inf_new,
        date_infectious_start = T_inf_start,
        date_infectious_end = T_resolved_new,
        date_resolved = T_resolved_new,
        date_resolved_final = T_resolved_final,
        stringsAsFactors = FALSE
      )
    }
  }
  if (n_fallback > 0) {
    message(sprintf("%d episode(s) were not in episode_timelines; used the random interval imputation for them.", n_fallback))
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
    dplyr::mutate(
      hh_sort_num = suppressWarnings(as.numeric(gsub("\\D", "", hh_id))),
      p_sort_num  = suppressWarnings(as.numeric(person_id))
    ) %>%
    dplyr::arrange(hh_sort_num, hh_id, p_sort_num, person_id) %>%
    dplyr::select(-hh_sort_num, -p_sort_num) %>%
    dplyr::left_join(person_span, by = "participantid")

  df_model_full <- df_meta_all %>%
    dplyr::left_join(df_imputed, by = "participantid")

  # =========================================================
  # STRICT SORTING & RISK WINDOW CALCULATION
  # =========================================================
  df_model_full <- df_model_full %>%
    dplyr::mutate(
      hh_sort_num = suppressWarnings(as.numeric(gsub("\\D", "", hh_id))),
      p_sort_num  = suppressWarnings(as.numeric(gsub("\\D", "", person_id)))
    ) %>%
    dplyr::arrange(hh_sort_num, hh_id, p_sort_num, person_id, date_infection) %>%
    dplyr::group_by(participantid) %>%
    dplyr::mutate(
      # the first unit of a person opens at the first sample; later ones one day after
      # the previous episode's resolution (immunity tail included)
      prev_resolved_final = dplyr::lag(date_resolved_final,
                                       default = dplyr::first(first_obs_date) - 1),
      start_risk_date = dplyr::case_when(
        is.na(episode_id) ~ first_obs_date,
        TRUE ~ prev_resolved_final + 1
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(hh_sort_num, hh_id, p_sort_num, person_id, start_risk_date) %>%
    dplyr::mutate(
      start_idx_raw = as.integer(start_risk_date - study_start_date) + 1L,
      end_idx_raw   = as.integer(last_obs_date  - study_start_date) + 1L,
      inf_idx_raw   = as.integer(date_infection - study_start_date) + 1L,
      i_idx = dplyr::row_number() # <-- THIS IS THE ANCHOR FOR EVERY MATRIX
    ) %>%
    dplyr::select(-hh_sort_num, -p_sort_num)

  N <- nrow(df_model_full)
  unique_hh <- unique(df_model_full$hh_id)
  H <- length(unique_hh)

  df_model_full$hh_id_int <- as.integer(factor(df_model_full$hh_id, levels = unique_hh))
  df_model_full$p_id_int  <- as.integer(as.factor(df_model_full$participantid))

  # --- AT-RISK WINDOW (start_risk / end_risk) ---
  start_raw <- df_model_full$start_idx_raw
  end_raw   <- df_model_full$end_idx_raw
  start_risk <- pmin(pmax(start_raw, 1L), T_max)
  end_risk   <- pmin(end_raw, T_max)
  end_risk[end_raw < 1L] <- 0L
  # window opens after the study ends -> not at risk
  end_risk[start_raw > T_max] <- pmax(start_risk[start_raw > T_max] - 1L, 0L)
  # prevalent infection: infected before the unit's window opened -> EMPTY window
  prevalent <- !is.na(df_model_full$inf_idx_raw) & df_model_full$inf_idx_raw < start_raw
  end_risk[prevalent] <- pmax(start_risk[prevalent] - 1L, 0L)
  df_model_full$start_risk <- as.integer(start_risk)
  df_model_full$end_risk   <- as.integer(end_risk)
  df_model_full$prevalent  <- prevalent

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

  # observed viral values, by person (written over the imputed curve on tested days)
  obs_df <- df_clean[!is.na(df_clean$ct_value) & !is.na(df_clean$date),
                     c("participantid", "date", "ct_value")]
  obs_by_pid <- split(obs_df, obs_df$participantid)

  # Role-level curve helper (single source of truth: vl_comput_fn.R)
  calc_curve_val <- function(t, p, is_ct) {
    if (is_ct) get_theoretical_ct(t, p) else get_theoretical_vl(t, p)
  }

  for (i in seq_len(N)) {
    row <- df_model_full[i, ]

    if (!is.na(row$episode_id)) {
      idx_inf <- as.integer(row$date_infection - study_start_date) + 1
      idx_ys  <- as.integer(row$date_infectious_start - study_start_date) + 1
      idx_end <- as.integer(row$date_infectious_end - study_start_date) + 1

      if (idx_inf >= 1 && idx_inf <= T_max) I[i, idx_inf] <- 1L
      y_start <- max(1, idx_ys); y_end <- min(T_max, idx_end)
      if (y_end < y_start) next
      Y[i, y_start:y_end] <- 1L

      # ---- viral values on the infectious days ----
      p_role <- row$role
      days_seq <- y_start:y_end
      t_vals <- days_seq - idx_inf          # days since infection
      imputed_vals <- NULL

      # A. Ct Data
      if (is_ct_data) {
        p_params <- imputation_params[[p_role]]
        if(is.null(p_params)) p_params <- list(Cpeak=33, r=1.5, d=1.2, t_peak=5)
        imputed_vals <- calc_curve_val(t_vals, p_params, TRUE)

        # B. Viral Load Data
      } else {
        if (model_type == "ODE") {
          current_ode_params <- ODE_params_list
          if(is.null(current_ode_params)) {
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
          imputed_vals <- rep(NA_real_, length(t_vals))
          ok <- traj_indices >= 1 & traj_indices <= length(traj)
          imputed_vals[ok] <- traj[traj_indices[ok]]
        } else {
          p_params <- imputation_params[[p_role]]
          if(is.null(p_params)) p_params <- list(v_p=5, t_p=4, lambda_g=2.8, lambda_d=1.0)
          imputed_vals <- calc_curve_val(t_vals, p_params, FALSE)
        }
      }

      # ---- Stage 1 episode-specific curve overrides the role-level curve ----
      if (!is.null(curve_lookup)) {
        cv <- curve_lookup[[paste(row$participantid, row$episode_id, sep = "||")]]
        if (!is.null(cv)) {
          cv_idx <- as.integer(cv$date - study_start_date) + 1L
          m <- match(days_seq, cv_idx)
          has <- !is.na(m)
          if (any(has)) {
            if (is.null(imputed_vals)) imputed_vals <- rep(NA_real_, length(days_seq))
            imputed_vals[has] <- cv$value[m[has]]
          }
        }
      }

      # ---- Fill matrix ----
      if (!is.null(imputed_vals) && length(imputed_vals) == length(days_seq)) {
        fill <- !is.na(imputed_vals)
        V[i, days_seq[fill]] <- imputed_vals[fill]
      }

      # ---- OVERWRITE WITH OBSERVED DATA ----
      obs <- obs_by_pid[[row$participantid]]
      if (!is.null(obs) && nrow(obs) > 0) {
        d_idx <- as.integer(obs$date - study_start_date) + 1L
        ok <- d_idx >= 1 & d_idx <= T_max
        V[i, d_idx[ok]] <- obs$ct_value[ok]
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

    # Actual physical person IDs, to detect several episodes of the same person
    member_pids    <- members$participantid

    # Generate Matrix based on user input or default
    if(!is.null(role_mixing_matrix)) {
      mat <- generate_contact_matrix_from_roles(member_roles, role_mixing_matrix)
    } else {
      # Default: Homogeneous Mixing
      n_h <- length(member_indices)
      mat <- matrix(1, n_h, n_h); diag(mat) <- 0
    }

    # --- ZERO OUT SELF-CONTACT ACROSS EPISODES ---
    for(r in 1:nrow(mat)) {
      for(c in 1:ncol(mat)) {
        # If the target (row) and source (col) are the EXACT SAME physical person
        if (member_pids[r] == member_pids[c]) {
          mat[r, c] <- 0.0
        }
      }
    }

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

  hh_size_df <- df_model_full %>%
    dplyr::group_by(hh_id_int) %>%
    dplyr::summarise(n = dplyr::n_distinct(participantid)) %>%
    dplyr::arrange(hh_id_int)

  hh_size_people <- hh_size_df$n

  # --- SEASONALITY ---
  R_roles <- length(role_levels)
  seasonal_forcing_mat <- matrix(1.0, nrow = T_max, ncol = R_roles)
  if (!is.null(surveillance_df)) {
    if (!all(c("date", "cases") %in% names(surveillance_df))) stop("surveillance_df must have 'date' and 'cases'")
    surveillance_df$date <- as.Date(surveillance_df$date)
    target_dates <- seq(from = study_start_date, to = study_end_date, by = "day")
    interp_res <- stats::approx(x = surveillance_df$date, y = surveillance_df$cases, xout = target_dates, method = "linear", rule = 2)
    daily_vec <- interp_res$y
    if(max(daily_vec, na.rm=TRUE) > 0) daily_vec <- daily_vec / max(daily_vec, na.rm=TRUE)
    daily_vec[is.na(daily_vec)] <- 0
    for(k in seq_len(R_roles)) seasonal_forcing_mat[,k] <- daily_vec
  } else if(!is.null(seasonal_forcing_list)) {
    # one numeric vector (length = number of study days) per role; matched by name when
    # the list is named, otherwise by position (order of role_levels)
    for (k in seq_len(R_roles)) {
      v <- if (!is.null(names(seasonal_forcing_list)) && role_levels[k] %in% names(seasonal_forcing_list)) {
        seasonal_forcing_list[[role_levels[k]]]
      } else if (k <= length(seasonal_forcing_list)) seasonal_forcing_list[[k]] else NULL
      if (!is.null(v)) {
        if (length(v) != T_max) stop(sprintf("seasonal_forcing_list entry for '%s' has length %d, expected %d (study days).",
                                             role_levels[k], length(v), T_max))
        seasonal_forcing_mat[, k] <- v
      }
    }
  }

  # =========================================================
  # 4b. FLAG INERT (PRIOR-ONLY) PARAMETERS
  # =========================================================
  # Stan declares every parameter unconditionally, so parameters not used under
  # the current configuration are still sampled from their prior and would show
  # up in print(fit) as if estimated. Record which are inert so that
  # fit_household_model() can exclude them from the returned stanfit.
  if (as.integer(use_vl_data) == 1) {
    inert_pars <- c("gen_shape", "gen_rate")
  } else {
    inert_pars <- c("vl_midpoint", "vl_slope")
    if (as.integer(use_curve_logic) == 0) {
      inert_pars <- c(inert_pars, "gen_shape", "gen_rate")
    }
  }

  # Per-unit table for diagnostics (not passed to Stan)
  units <- data.frame(
    i_idx = df_model_full$i_idx, participantid = df_model_full$participantid,
    hh_id = df_model_full$hh_id, role = df_model_full$role,
    episode_id = df_model_full$episode_id,
    date_infection = df_model_full$date_infection,
    date_infectious_start = df_model_full$date_infectious_start,
    date_infectious_end = df_model_full$date_infectious_end,
    date_resolved_final = df_model_full$date_resolved_final,
    start_risk = df_model_full$start_risk, end_risk = df_model_full$end_risk,
    prevalent = df_model_full$prevalent,
    stringsAsFactors = FALSE
  )

  # =========================================================
  # 5. RETURN LIST
  # =========================================================
  list(
    N = N, T = T_max, H = H, R = R_roles, delta = delta,
    hh_id = df_model_full$hh_id_int,
    role_id = match(df_model_full$role, role_levels),
    I = I, Y = Y, V = V,
    start_risk = df_model_full$start_risk,
    end_risk   = df_model_full$end_risk,

    # --- CONTACT DATA ---
    n_contacts  = length(contact_w),
    contact_src = contact_src,
    contact_tgt = contact_tgt,
    contact_w   = contact_w,

    p_id = df_model_full$p_id_int,
    hh_size_people = as.integer(hh_size_people),
    hh_max_size = as.integer(hh_max_size),
    hh_members = hh_members,
    seasonal_forcing_mat = seasonal_forcing_mat,
    reference_phi = 1.0, reference_kappa = 1.0,
    use_vl_data = as.integer(use_vl_data),

    vl_type = as.integer(detected_vl_type),
    use_curve_logic = as.integer(use_curve_logic),

    inert_pars = inert_pars,
    units = units,

    K_susc = K_susc, X_susc = X_susc,
    K_inf  = K_inf,  X_inf  = X_inf,

    prior_beta1_type = p_beta1$type, prior_beta1_params = p_beta1$params,
    prior_beta2_type = p_beta2$type, prior_beta2_params = p_beta2$params,
    prior_alpha_type = p_alpha$type, prior_alpha_params = p_alpha$params,
    prior_cov_type   = p_cov$type,   prior_cov_params   = p_cov$params,

    prior_shape_type = p_shape$type, prior_shape_params = p_shape$params,
    prior_rate_type  = p_rate$type,  prior_rate_params  = p_rate$params,
    prior_vl_midpoint_type   = p_vl_midpoint$type,   prior_vl_midpoint_params = p_vl_midpoint$params,
    prior_vl_slope_type      = p_vl_slope$type,       prior_vl_slope_params    = p_vl_slope$params,
    prior_phi_role_type     = p_phi_role$type,   prior_phi_role_params   = p_phi_role$params,
    prior_kappa_role_type   = p_kappa_role$type, prior_kappa_role_params = p_kappa_role$params
  )
}
