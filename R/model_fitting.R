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
  if (is.null(init_fun)) {
    init_fun <- function() {
      R <- stan_data$R
      K_susc <- stan_data$K_susc
      K_inf <- stan_data$K_inf

      init_list <- list(
        # Core transmission parameters
        log_beta1 = log(0.008),
        log_beta2 = log(0.008),
        log_alpha_comm = log(5e-4),

        # Role-specific effects
        log_phi_by_role_raw = rep(0.1, max(1, R - 1)),
        log_kappa_by_role_raw = rep(0.1, max(1, R - 1)),

        # ALL viral dynamics parameters (always declared in Stan)
        gen_shape = 3.0,
        gen_rate = 1.0,
        Ct50 = 35.0,
        slope_ct = 2.0,
        V_ref = 3.0,                      # ← ADD THIS
        rho = 2.5                         # ← ADD THIS
      )

      # Conditional covariate effects
      if (K_susc > 0) init_list$beta_susc <- rep(0.0, K_susc)
      if (K_inf > 0) init_list$beta_inf <- rep(0.0, K_inf)

      return(init_list)
    }
  }

  # Run sampler (unchanged)
  out <- rstan::sampling(
    object = stanmodels$household_transmission,
    data = stan_data,
    iter = iter,
    chains = chains,
    warmup = warmup,
    init = init_fun,
    control = list(adapt_delta = 0.95, max_treedepth = 15),
    ...
  )

  return(out)
}
