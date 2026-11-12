# =============================================================================
# STAGE 1: viral kinetics + infection timing
#
#   prepare_kinetics_data()        build the Stage 1 Stan data from the observation table
#   fit_viral_kinetics()           fit inst/stan/viral_kinetics.stan
#   summarize_kinetics()           group-level kinetic parameters (median, 95% CrI)
#   kinetics_draw_to_timelines()   one posterior draw -> episode_timelines + episode_curves
#                                  for prepare_stan_data() (Stage 2)
#
# Workflow (sequential multiple imputation):
#   kin  <- prepare_kinetics_data(df, viral_scale = "Ct", lod = 45, ...)
#   kfit <- fit_viral_kinetics(kin)
#   ex   <- kinetics_extract(kfit)
#   for (m in draw_ids) {
#     tl <- kinetics_draw_to_timelines(kfit, m, ex)
#     si <- prepare_stan_data(df, ..., episode_timelines = tl$episode_timelines,
#                                      episode_curves    = tl$episode_curves)
#     fit_m <- fit_household_model(si)          # then pool the posterior draws over m
#   }
# =============================================================================

#' Expected viral value at s days since infection
#'
#' Mirrors \code{mean_curve()} in \code{viral_kinetics.stan}.
#' @param s Days since infection (numeric vector).
#' @param cp Curve parameters. Ct (\code{vl_type = 0}): \code{(Cpeak, r, d, t_peak)}.
#'   log10 viral load (\code{vl_type = 1}): \code{(v_p, t_p, lambda_g, lambda_d)}.
#' @param vl_type 0 = Ct, 1 = log10 viral load.
#' @return Numeric vector of expected values on the native scale.
#' @export
kinetics_curve <- function(s, cp, vl_type) {
  cp <- as.numeric(cp)
  if (vl_type == 0) {
    ifelse(s <= cp[4], cp[1] + cp[2] * (cp[4] - s), cp[1] + cp[3] * (s - cp[4]))
  } else {
    a <- -cp[3] * (s - cp[2]); b <- cp[4] * (s - cp[2])
    m <- pmax(a, b)
    log10(2) + cp[1] - (m + log(exp(a - m) + exp(b - m))) / log(10)
  }
}

#' Infectious window from an expected viral curve
#'
#' The rule used to turn a Stage 1 trajectory into the days a person is infectious
#' (all days relative to the day of the first positive sample, which is day 0):
#' \itemize{
#'   \item the infectious period \strong{starts} the first day the expected curve is detectable
#'     (Ct <= lod, or log10 VL >= lod), but not before the day after the last negative sample
#'     before the first positive, and never after the first positive;
#'   \item it \strong{ends} the last day the expected curve is detectable, but not after the day
#'     before the first negative sample following the last positive, and never before the last
#'     positive (observed positives are always included);
#'   \item infection is \code{max(min_incub_days, round(tinc))} days before the first positive.
#'     Infectiousness starts at least one day after infection, as in the simulator with
#'     \code{infectious_period = "viral"}; the days in between are infected but not infectious.
#' }
#' @param cp,tinc,vl_type,lod Curve parameters, days from infection to first positive, scale, threshold.
#' @param last_pos Day of the last positive sample (relative to the first positive).
#' @param prev_neg,next_neg Days of the last negative before the first positive and the first negative
#'   after the last positive (relative to the first positive); \code{NA} if none.
#' @param min_incub_days Minimum days between infection and first positive (default 1, which matches
#'   the simulator with \code{infectious_period = "viral"}).
#' @param horizon_after Days after the last positive to search for the end of detectability.
#' @return list(infection, start, end), days relative to the first positive.
#' @export
infectious_window_from_curve <- function(cp, tinc, vl_type, lod, last_pos = 0,
                                         prev_neg = NA, next_neg = NA,
                                         min_incub_days = 1, horizon_after = 60) {
  grid <- seq(floor(-tinc) - 2, ceiling(last_pos) + horizon_after)
  mu   <- kinetics_curve(grid + tinc, cp, vl_type)
  det  <- if (vl_type == 0) mu <= lod else mu >= lod

  if (any(det)) {
    start_c <- min(grid[det]); end_c <- max(grid[det])
  } else {                      # curve never detectable: fall back to the observed positives
    start_c <- 0; end_c <- last_pos
  }
  lower_ok <- if (is.na(prev_neg)) -Inf else prev_neg + 1
  upper_ok <- if (is.na(next_neg))  Inf else next_neg - 1

  inf_rel <- -max(min_incub_days, round(tinc))
  start <- min(max(start_c, lower_ok), 0)     # include the first positive
  start <- min(max(start, inf_rel + 1), 0)    # infectious from at least the day after infection
  end   <- max(min(end_c, upper_ok), last_pos)  # include the last positive
  list(infection = inf_rel, start = start, end = end)
}

#' Prepare the Stage 1 (viral kinetics) data
#'
#' Builds one record per episode from the same observation table used for
#' \code{prepare_stan_data()}. For each episode the samples used are its own positives plus the
#' negatives before the first positive (up to \code{pre_window} days, after the previous episode)
#' and after the last positive (up to \code{post_window} days, before the next episode). Negatives
#' enter as \emph{censored} observations, which is what lets the data inform infection timing.
#'
#' @param df_clean Observation data (same format as for \code{prepare_stan_data()}); must have
#'   \code{episode_id}, a date or \code{day_index}, the viral value, and ideally \code{test_result}.
#' @param study_start_date Date of \code{day_index == 1} (only used when \code{day_index} is given).
#' @param viral_scale "auto", "log10" or "Ct".
#' @param lod Detection threshold on the native scale. A sample is negative when the Ct is above
#'   \code{lod} (Ct) or the log10 viral load is below it. Defaults to the simulator's thresholds
#'   (40 for Ct, 1e-6 for log10); use e.g. 45 for Ct data with a 45-cycle cutoff. If the data have a
#'   \code{test_result} column it defines positivity and \code{lod} only enters the censoring.
#' @param role_levels Roles in the data.
#' @param groups Optional named character vector mapping roles to groups, e.g.
#'   \code{c(elderly = "adult")} to pool elderly with adults. Unmapped roles form their own group.
#' @param pre_window,post_window Days of negatives kept before the first / after the last positive.
#' @param min_samples,min_positive Episodes with fewer samples / positives are not fitted in Stage 1;
#'   they use the group-level trajectory in \code{kinetics_draw_to_timelines()}.
#' @param symptom_onset Optional data frame with \code{participantid} and \code{symptom_date}.
#' @param sx_window Use a symptom date only if it is within this many days of the first positive.
#' @param sx_incubation Mean and SD (days) of the incubation period to symptom onset.
#' @param nu Student-t degrees of freedom of the measurement error (<= 0 for a normal likelihood).
#' @param sigma_floor,prior_sigma_y Lower bound and half-normal scale of the measurement SD.
#' @param priors Optional overrides, each a numeric vector of length 5, ordered
#'   (level = Cpeak or v_p, log rise rate, log decay rate, c = peak day relative to the first
#'   positive sample, log tinc): \code{mu} and \code{mu_sd} (population means), \code{tau_grp}
#'   and \code{tau_ep} (half-normal scales of the group- and episode-level SDs).
#' @param tinc_bounds,peak_bounds Soft bounds (days) on infection-to-first-positive and
#'   infection-to-peak. The default lower bound of 1 day matches the simulator with \code{infectious_period = "viral"}.
#' @param tend_max,tend_soft Soft upper bound (days from infection) on the end of the detectable
#'   period, and its softness.
#' @param bound_soft Softness (days) of the tinc / peak bounds.
#' @param no_prior_neg_sd_mult Multiplier on the tinc spread for episodes with no earlier negative.
#' @return An object of class \code{hh_kinetics_data}.
#' @export
prepare_kinetics_data <- function(df_clean,
                                  study_start_date = as.Date("2024-07-01"),
                                  viral_scale = c("auto", "log10", "Ct"),
                                  lod = NULL,
                                  role_levels = c("adult", "infant", "toddler", "elderly"),
                                  groups = NULL,
                                  pre_window = 21, post_window = 28,
                                  min_samples = 3, min_positive = 2,
                                  symptom_onset = NULL, sx_window = 14, sx_incubation = c(4, 1),
                                  nu = 4, sigma_floor = 0.1, prior_sigma_y = 2,
                                  priors = list(),
                                  tinc_bounds = c(1, 21), peak_bounds = c(1, 14),
                                  tend_max = 28, tend_soft = 7, bound_soft = 0.25,
                                  no_prior_neg_sd_mult = 1.5) {

  sc <- .resolve_viral_scale(df_clean, viral_scale)
  is_ct <- sc$is_ct; vl_type <- sc$vl_type
  if (is.null(lod)) {
    lod <- if (is_ct) 40 else 1e-6
    message(sprintf("lod not supplied: using %s (the simulator's default threshold). Set lod explicitly for real data.",
                    format(lod)))
  }

  df <- .standardize_clean(df_clean, study_start_date)
  miss <- setdiff(c("participantid", "role_name", "date", "ct_value", "episode_id"), names(df))
  if (length(miss) > 0) stop("df_clean is missing column(s): ", paste(miss, collapse = ", "))

  vals <- df$ct_value
  pos <- if (".test_pos" %in% names(df)) df$.test_pos else if (is_ct) vals <= lod else vals >= lod
  df$.pos <- !is.na(pos) & pos & !is.na(vals)

  map_group <- function(role) {
    if (is.na(role) || !role %in% role_levels) role <- "adult"
    if (!is.null(groups) && role %in% names(groups)) groups[[role]] else role
  }

  ep_rows <- list(); obs_rows <- list()
  df_by_pid <- split(df, df$participantid)

  for (pid in names(df_by_pid)) {
    p <- df_by_pid[[pid]]
    p <- p[order(p$date), , drop = FALSE]
    eids <- unique(p$episode_id[which(p$episode_id > 0)])
    if (length(eids) == 0) next
    role <- p$role_name[1]; grp <- map_group(role)
    hh <- if ("familyidstars" %in% names(p)) p$familyidstars[1] else NA_character_

    info <- lapply(eids, function(e) {
      r  <- p[which(p$episode_id == e), , drop = FALSE]
      pr <- r[r$.pos, , drop = FALSE]
      if (nrow(pr) == 0) pr <- r          # no quantified positive: use the episode span
      list(eid = e, first = min(pr$date), last = max(pr$date))
    })
    info <- info[order(vapply(info, function(z) as.numeric(z$first), numeric(1)))]
    negs <- p[!p$.pos & !is.na(p$date), , drop = FALSE]

    for (k in seq_along(info)) {
      z <- info[[k]]
      prev_last  <- if (k > 1) info[[k - 1]]$last else as.Date(NA)
      next_first <- if (k < length(info)) info[[k + 1]]$first else as.Date(NA)

      pn <- negs$date[negs$date < z$first]
      nn <- negs$date[negs$date > z$last]
      prev_neg <- if (length(pn)) max(pn) else as.Date(NA)
      next_neg <- if (length(nn)) min(nn) else as.Date(NA)

      lo <- z$first - pre_window
      if (!is.na(prev_last)) lo <- max(lo, prev_last + 1)
      hi <- z$last + post_window
      if (!is.na(next_first)) hi <- min(hi, next_first - 1)
      w <- p[which(p$date >= lo & p$date <= hi), , drop = FALSE]

      n_obs <- nrow(w); n_pos <- sum(w$.pos)
      has_pre_neg <- any(!w$.pos & w$date < z$first)
      included <- (n_obs >= min_samples) && (n_pos >= min_positive)

      key <- paste(pid, z$eid, sep = "||")
      ep_rows[[length(ep_rows) + 1]] <- data.frame(
        key = key, participantid = pid, hh_id = hh, episode_id = z$eid, role = role, group = grp,
        first_pos_date = z$first, last_pos_date = z$last,
        prev_neg_date = prev_neg, next_neg_date = next_neg,
        n_obs = n_obs, n_pos = n_pos, has_pre_neg = has_pre_neg, included = included,
        stringsAsFactors = FALSE)

      if (included) {
        obs_rows[[length(obs_rows) + 1]] <- data.frame(
          key = key,
          t_rel = as.numeric(w$date - z$first),
          y = ifelse(w$.pos, w$ct_value, lod),
          cens = as.integer(!w$.pos),
          stringsAsFactors = FALSE)
      }
    }
  }

  if (length(ep_rows) == 0) stop("No episodes found: check that df_clean has episode_id > 0 rows.")
  episodes <- do.call(rbind, ep_rows)
  rownames(episodes) <- NULL
  if (!any(episodes$included)) {
    stop(sprintf("No episode has >= %d samples and >= %d positives; relax min_samples / min_positive.",
                 min_samples, min_positive))
  }
  episodes$ep_idx <- NA_integer_
  episodes$ep_idx[episodes$included] <- seq_len(sum(episodes$included))

  obs <- do.call(rbind, obs_rows)
  obs$ep_idx <- episodes$ep_idx[match(obs$key, episodes$key)]
  obs <- obs[order(obs$ep_idx, obs$t_rel), ]

  # groups (in role order), for ALL episodes so excluded ones can use group trajectories
  group_levels <- unique(vapply(role_levels, map_group, character(1)))
  group_levels <- unique(c(group_levels, setdiff(unique(episodes$group), group_levels)))
  inc <- episodes[episodes$included, ]
  grp_idx <- match(inc$group, group_levels)

  # optional symptom-onset information
  has_sx <- integer(nrow(inc)); sx_lag <- numeric(nrow(inc))
  if (!is.null(symptom_onset)) {
    so <- as.data.frame(symptom_onset, stringsAsFactors = FALSE)
    if (!all(c("participantid", "symptom_date") %in% names(so))) {
      stop("symptom_onset needs columns 'participantid' and 'symptom_date'.")
    }
    so$symptom_date <- as.Date(so$symptom_date)
    for (j in seq_len(nrow(inc))) {
      sd_j <- so$symptom_date[so$participantid == inc$participantid[j] & !is.na(so$symptom_date)]
      if (length(sd_j) == 0) next
      lag_j <- as.numeric(inc$first_pos_date[j] - sd_j)
      best <- which.min(abs(lag_j))
      if (abs(lag_j[best]) <= sx_window) { has_sx[j] <- 1L; sx_lag[j] <- lag_j[best] }
    }
  }

  # default priors by scale. Order of the five entries:
  #   (level [Cpeak | v_p], log rise rate [r | lambda_g], log decay rate [d | lambda_d],
  #    c = peak day relative to the first positive, log tinc)
  def <- if (is_ct) {
    list(mu = c(33, log(1.5), log(1.2), 2, log(4)),
         mu_sd = c(3, 0.5, 0.5, 3, 0.6),
         tau_grp = c(2, 0.3, 0.3, 1, 0.3), tau_ep = c(2, 0.3, 0.3, 1, 0.3))
  } else {
    list(mu = c(5, log(2.5), log(1.5), 2, log(3)),
         mu_sd = c(1.5, 0.5, 0.5, 3, 0.6),
         tau_grp = c(1, 0.3, 0.3, 1, 0.3), tau_ep = c(1, 0.3, 0.3, 1, 0.3))
  }
  pr <- def
  for (nm in names(def)) if (!is.null(priors[[nm]])) {
    if (length(priors[[nm]]) != 5) stop("priors$", nm, " must have length 5.")
    pr[[nm]] <- as.numeric(priors[[nm]])
  }

  stan_data <- list(
    N_ep = nrow(inc), N_obs = nrow(obs), G = length(group_levels),
    vl_type = vl_type, lod = lod,
    ep = as.array(as.integer(obs$ep_idx)), t_rel = as.array(obs$t_rel),
    y = as.array(obs$y), cens = as.array(as.integer(obs$cens)),
    grp = as.array(as.integer(grp_idx)),
    tinc_sd_mult = as.array(ifelse(inc$has_pre_neg, 1, no_prior_neg_sd_mult)),
    has_sx = as.array(as.integer(has_sx)), sx_lag = as.array(sx_lag),
    sx_mean = sx_incubation[1], sx_sd = sx_incubation[2],
    nu = nu, sigma_floor = sigma_floor, prior_sigma_y = prior_sigma_y,
    prior_mu = pr$mu, prior_mu_sd = pr$mu_sd,
    prior_tau_grp = pr$tau_grp, prior_tau_ep = pr$tau_ep,
    tinc_bounds = tinc_bounds, peak_bounds = peak_bounds,
    tend_max = tend_max, tend_soft = tend_soft, bound_soft = bound_soft
  )

  message(sprintf("Stage 1 data: %d episodes fitted (%d excluded -> group trajectory), %d samples (%d negative/censored), scale = %s, lod = %s.",
                  nrow(inc), sum(!episodes$included), nrow(obs), sum(obs$cens), sc$scale, format(lod)))

  structure(list(stan_data = stan_data, episodes = episodes, scale = sc$scale,
                 vl_type = vl_type, lod = lod, study_start_date = as.Date(study_start_date),
                 group_levels = group_levels, role_levels = role_levels),
            class = "hh_kinetics_data")
}

#' Fit the Stage 1 viral kinetics model
#'
#' @param kin_data Output of \code{prepare_kinetics_data()}.
#' @param iter,chains,warmup Sampler settings.
#' @param init_fun Optional initial values; defaults to a small jitter around the priors.
#' @param ... Passed to \code{rstan::sampling()} (\code{cores}, \code{seed}, \code{control}, ...).
#' @return An object of class \code{hh_kinetics_fit} (a list with the \code{stanfit} and the data).
#' @export
fit_viral_kinetics <- function(kin_data, iter = 2000, chains = 4, warmup = 1000,
                               init_fun = NULL, ...) {
  sd <- kin_data$stan_data
  if (is.null(init_fun)) {
    init_fun <- function() {
      list(mu_pop  = sd$prior_mu + stats::rnorm(5, 0, 0.02),
           tau_grp = rep(0.1, 5), tau_ep = rep(0.1, 5),
           z_grp   = matrix(stats::rnorm(sd$G * 5, 0, 0.05), sd$G, 5),
           z_ep    = matrix(stats::rnorm(sd$N_ep * 5, 0, 0.05), sd$N_ep, 5),
           sigma_y = max(sd$sigma_floor * 2, 1))
    }
  }
  call_args <- list(
    object = stanmodels$viral_kinetics, data = sd,
    iter = iter, chains = chains, warmup = warmup, init = init_fun,
    pars = c("z_grp", "z_ep", "theta"), include = FALSE,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  call_args <- utils::modifyList(call_args, list(...))
  fit <- do.call(rstan::sampling, call_args)
  structure(list(fit = fit, data = kin_data), class = "hh_kinetics_fit")
}

#' Extract the posterior arrays needed to build timelines
#' @param kin_fit Output of \code{fit_viral_kinetics()}.
#' @return list(cp, tinc, grp_curve, grp_tinc, sigma_y) of posterior draws.
#' @export
kinetics_extract <- function(kin_fit) {
  rstan::extract(kin_fit$fit, pars = c("cp", "tinc", "grp_curve", "grp_tinc", "sigma_y"))
}

#' Group-level kinetic parameters from the Stage 1 fit
#'
#' @param kin_fit Output of \code{fit_viral_kinetics()}.
#' @param probs Quantiles (lower, median, upper).
#' @return Data frame with group, parameter, median, lwr, upr.
#' @export
summarize_kinetics <- function(kin_fit, probs = c(0.025, 0.5, 0.975)) {
  ex <- rstan::extract(kin_fit$fit, pars = c("grp_curve", "grp_tinc"))
  kd <- kin_fit$data
  pn <- if (kd$vl_type == 0) c("Cpeak", "r", "d", "t_peak") else c("v_p", "t_p", "lambda_g", "lambda_d")
  rows <- list()
  for (g in seq_along(kd$group_levels)) {
    for (k in 1:4) {
      q <- stats::quantile(ex$grp_curve[, g, k], probs, names = FALSE)
      rows[[length(rows) + 1]] <- data.frame(group = kd$group_levels[g], parameter = pn[k],
                                             median = q[2], lwr = q[1], upr = q[3])
    }
    q <- stats::quantile(ex$grp_tinc[, g], probs, names = FALSE)
    rows[[length(rows) + 1]] <- data.frame(group = kd$group_levels[g], parameter = "t_incub",
                                           median = q[2], lwr = q[1], upr = q[3])
  }
  do.call(rbind, rows)
}

#' One Stage 1 posterior draw -> timelines for Stage 2
#'
#' For draw \code{draw}, every episode gets its infection day, infectious start and end, and the
#' expected viral value on each infectious day. Episodes not fitted in Stage 1 (too few samples)
#' use the group-level trajectory and incubation from the same draw. Pass the result to
#' \code{prepare_stan_data(episode_timelines = , episode_curves = )}.
#'
#' @param kin_fit Output of \code{fit_viral_kinetics()}.
#' @param draw Index of the posterior draw (1..number of post-warmup draws).
#' @param ex Optional output of \code{kinetics_extract()} (pass it when looping over many draws).
#' @param min_incub_days Minimum days between infection and the first positive (default 1, which
#'   matches the simulator with \code{infectious_period = "viral"}).
#' @return list(episode_timelines, episode_curves).
#' @export
kinetics_draw_to_timelines <- function(kin_fit, draw, ex = NULL, min_incub_days = 1) {
  kd <- kin_fit$data
  if (is.null(ex)) ex <- kinetics_extract(kin_fit)
  if (draw < 1 || draw > dim(ex$cp)[1]) stop("draw must be between 1 and ", dim(ex$cp)[1])
  eps <- kd$episodes

  tl <- vector("list", nrow(eps)); cv <- vector("list", nrow(eps))
  for (r in seq_len(nrow(eps))) {
    e <- eps[r, ]
    if (isTRUE(e$included)) {
      cp <- ex$cp[draw, e$ep_idx, ]; tinc <- ex$tinc[draw, e$ep_idx]
    } else {
      g <- match(e$group, kd$group_levels)
      cp <- ex$grp_curve[draw, g, ]; tinc <- ex$grp_tinc[draw, g]
    }
    last_rel <- as.numeric(e$last_pos_date - e$first_pos_date)
    prev_rel <- if (is.na(e$prev_neg_date)) NA else as.numeric(e$prev_neg_date - e$first_pos_date)
    next_rel <- if (is.na(e$next_neg_date)) NA else as.numeric(e$next_neg_date - e$first_pos_date)
    w <- infectious_window_from_curve(cp, tinc, kd$vl_type, kd$lod, last_rel, prev_rel, next_rel,
                                      min_incub_days)
    tl[[r]] <- data.frame(
      participantid = e$participantid, episode_id = e$episode_id,
      date_infection        = e$first_pos_date + w$infection,
      date_infectious_start = e$first_pos_date + w$start,
      date_infectious_end   = e$first_pos_date + w$end,
      stringsAsFactors = FALSE)
    d <- seq(w$start, w$end)
    cv[[r]] <- data.frame(
      participantid = e$participantid, episode_id = e$episode_id,
      date = e$first_pos_date + d,
      value = kinetics_curve(d + tinc, cp, kd$vl_type),
      stringsAsFactors = FALSE)
  }
  list(episode_timelines = do.call(rbind, tl), episode_curves = do.call(rbind, cv))
}
