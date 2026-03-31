#' Fit Household Transmission Model
#'
#' Fits the compiled Stan model to the prepared household data.
#'
#' @param stan_data A list of data formatted by \code{prepare_stan_data}.
#' @param iter Integer. Number of iterations per chain (including warmup). Defaults to 2000.
#' @param chains Integer. Number of Markov chains. Defaults to 4.
#' @param warmup Integer. Number of warmup iterations. Defaults to 1000.
#' @param init_fun Function or List. Initial values for the sampler. If NULL, uses robust defaults tailored for this model.
#' @param ... Additional arguments passed to \code{rstan::sampling} (e.g., \code{cores}, \code{seed}).
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
                                ...) {

  # 1. Define Default Initial Values
  # These defaults are crucial for convergence in transmission models to prevent
  # the sampler from starting in a region with 0 likelihood.
  if (is.null(init_fun)) {
    init_fun <- function() {
      # Get dimensions from stan_data
      R <- stan_data$R
      K_susc <- stan_data$K_susc
      K_inf <- stan_data$K_inf

      vl_midpoint_init <- if (stan_data$vl_type == 0) 33.0 else 6.0
      vl_slope_init    <- if (stan_data$vl_type == 0) 4.0  else 1.0


      init_list <- list(
        # Core transmission parameters (log scale for rates)
        log_beta1 = log(0.008),           # ~8e-3 baseline transmission
        log_beta2 = log(0.008),           # ~8e-3 VL-dependent transmission
        log_alpha_comm = log(5e-4),       # ~5e-4 community rate (CRITICAL!)

        # Role-specific effects (R-1 parameters, relative to reference)
        # Start at reference level (0 on log scale = 1.0 on natural scale)
        log_phi_by_role_raw = rep(0.1, max(1, R - 1)),     # Small positive to avoid exact 0
        log_kappa_by_role_raw = rep(0.1, max(1, R - 1)),   # Small positive to avoid exact 0

        # Viral dynamics parameters (must respect Stan bounds)
        gen_shape = 3.0,                  # Within bounds [1.0, 20.0]
        gen_rate = 1.0,                   # Within bounds [0.1, 5.0]

        vl_midpoint = vl_midpoint_init,   # RENAMED + data-aware
        vl_slope    = vl_slope_init       # RENAMED + data-aware
      )

      # Conditional covariate effects (only add if covariates are present)
      if (K_susc > 0) {
        init_list$beta_susc <- array(0.0, K_susc)  # Start at no effect
      }
      if (K_inf > 0) {
        init_list$beta_inf <- array(0.0, K_inf)    # Start at no effect
      }

      return(init_list)
    }
  }

  # 2. Run Sampler
  # NOTE: We use 'stanmodels$household_transmission' which is the pre-compiled object
  # created by the package. We do NOT refer to the .stan file path directly.
  out <- rstan::sampling(
    object = stanmodels$household_transmission,
    data = stan_data,
    iter = iter,
    chains = chains,
    warmup = warmup,
    init = init_fun,
    control = list(adapt_delta = 0.95, max_treedepth = 15),
    # We default to including all parameters.
    # Users can filter the output object later if needed.
    ...
  )

  return(out)
}
