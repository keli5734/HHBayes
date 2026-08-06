#' Fit Household Transmission Model
#'
#' Fits the compiled Stan model to the prepared household data.
#'
#' @param stan_data A list of data formatted by \code{prepare_stan_data}.
#' @param iter Integer. Number of iterations per chain (including warmup). Defaults to 2000.
#' @param chains Integer. Number of Markov chains. Defaults to 4.
#' @param warmup Integer. Number of warmup iterations. Defaults to 1000.
#' @param init_fun Function or List. Initial values for the sampler. If NULL, uses robust defaults tailored for this model.
#' @param drop_inert Logical. If TRUE (default), viral-component parameters that
#'   are inert under the current data configuration -- i.e. sampled only from
#'   their prior and never entering the likelihood -- are excluded from the
#'   returned \code{stanfit}, together with internal bookkeeping quantities.
#'   They are still sampled internally; this only removes them from the stored
#'   and printed object so they cannot be mistaken for estimates. Which set is
#'   inert is taken from \code{stan_data$inert_pars} (attached by
#'   \code{prepare_stan_data}); if that field is absent it is derived from the
#'   \code{use_vl_data} / \code{use_curve_logic} flags. Set FALSE to keep all
#'   parameters in the output.
#' @param ... Additional arguments passed to \code{rstan::sampling} (e.g.,
#'   \code{cores}, \code{seed}, \code{control}, or an explicit \code{pars} /
#'   \code{include}, which override the automatic dropping).
#'
#' @return A \code{stanfit} object containing the posterior samples.
#' @import rstan
#' @import methods
#' @export
fit_household_model <- function(stan_data,
                                iter = 2000,
                                chains = 4,
                                warmup = 1000,
                                init_fun = NULL,
                                drop_inert = TRUE,
                                ...) {

  # ---------------------------------------------------------------------------
  # 0. Separate non-Stan metadata from the data list.
  #    `inert_pars` is bookkeeping used here; it is NOT declared in the Stan
  #    `data` block, so it must be stripped before the data list is passed on.
  # ---------------------------------------------------------------------------
  inert_pars <- stan_data$inert_pars
  stan_data$inert_pars <- NULL

  # Backward compatibility: if an older prepare_stan_data() did not attach
  # $inert_pars, derive the inert set from the configuration flags.
  #   use_vl_data == 1                        -> gen_shape / gen_rate inert
  #   use_vl_data == 0 & use_curve_logic == 1 -> vl_midpoint / vl_slope inert
  #   use_vl_data == 0 & use_curve_logic == 0 -> all four inert (constant infectiousness)
  if (is.null(inert_pars)) {
    uvd <- if (!is.null(stan_data$use_vl_data)) stan_data$use_vl_data else 1L
    ucl <- if (!is.null(stan_data$use_curve_logic)) stan_data$use_curve_logic else 0L
    if (uvd == 1) {
      inert_pars <- c("gen_shape", "gen_rate")
    } else {
      inert_pars <- c("vl_midpoint", "vl_slope")
      if (ucl == 0) inert_pars <- c(inert_pars, "gen_shape", "gen_rate")
    }
  }

  # ---------------------------------------------------------------------------
  # 1. Default initial values.
  #    Crucial for convergence in transmission models: prevents the sampler
  #    from starting in a region with zero likelihood. (Providing inits for
  #    inert parameters is harmless.)
  # ---------------------------------------------------------------------------
  if (is.null(init_fun)) {
    init_fun <- function() {
      R      <- stan_data$R
      K_susc <- stan_data$K_susc
      K_inf  <- stan_data$K_inf

      vl_midpoint_init <- if (stan_data$vl_type == 0) 33.0 else 6.0
      vl_slope_init    <- if (stan_data$vl_type == 0) 4.0  else 1.0

      init_list <- list(
        # Core transmission parameters (log scale for rates)
        log_beta1      = log(0.008),   # ~8e-3 baseline transmission
        log_beta2      = log(0.008),   # ~8e-3 VL-dependent transmission
        log_alpha_comm = log(5e-4),    # ~5e-4 community rate (CRITICAL!)

        # Role-specific effects (R-1 parameters, relative to reference).
        # Start just off the reference level (0 on log scale = 1.0 natural).
        log_phi_by_role_raw   = rep(0.1, max(1, R - 1)),
        log_kappa_by_role_raw = rep(0.1, max(1, R - 1)),

        # Viral dynamics parameters (must respect Stan bounds)
        gen_shape = 3.0,               # within [1.0, 20.0]
        gen_rate  = 1.0,               # within [0.1, 5.0]

        vl_midpoint = vl_midpoint_init,  # data-aware
        vl_slope    = vl_slope_init      # data-aware
      )

      # Conditional covariate effects (only add if covariates are present)
      if (K_susc > 0) init_list$beta_susc <- array(0.0, K_susc)  # start at no effect
      if (K_inf  > 0) init_list$beta_inf  <- array(0.0, K_inf)   # start at no effect

      return(init_list)
    }
  }

  # ---------------------------------------------------------------------------
  # 2. Assemble the arguments for rstan::sampling().
  #    Using do.call() lets user-supplied ... args cleanly override the
  #    defaults (e.g. a custom control list, seed, cores) without triggering
  #    "duplicated argument" errors.
  # ---------------------------------------------------------------------------
  dots <- list(...)

  call_args <- list(
    object  = stanmodels$household_transmission,
    data    = stan_data,
    iter    = iter,
    chains  = chains,
    warmup  = warmup,
    init    = init_fun,
    control = list(adapt_delta = 0.95, max_treedepth = 15)
  )

  # Exclude inert (prior-only) parameters and internal bookkeeping quantities
  # from the OUTPUT, unless the user has supplied their own pars/include.
  if (isTRUE(drop_inert) && !any(c("pars", "include") %in% names(dots))) {
    call_args$pars <- unique(c(
      inert_pars,
      # internal bookkeeping the README already strips manually:
      "log_phi_by_role_raw", "log_kappa_by_role_raw",
      "log_beta1", "log_beta2", "log_alpha_comm",
      "g_curve_est", "V_term_calc"
    ))
    call_args$include <- FALSE
  }

  # User-supplied ... arguments take precedence over the defaults above.
  call_args <- utils::modifyList(call_args, dots)

  # ---------------------------------------------------------------------------
  # 3. Run the sampler.
  #    We use the pre-compiled 'stanmodels$household_transmission' object built
  #    by the package; we do NOT refer to the .stan file path directly.
  # ---------------------------------------------------------------------------
  out <- do.call(rstan::sampling, call_args)

  return(out)
}
