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
        # Always needed
        log_beta1 = log(0.008),
        log_beta2 = log(0.008),
        log_alpha_comm = log(5e-4),
        log_phi_by_role_raw = rep(0.1, max(1, R - 1)),
        log_kappa_by_role_raw = rep(0.1, max(1, R - 1))
      )

      # Add covariate effects if present
      if (K_susc > 0) init_list$beta_susc <- rep(0.0, K_susc)
      if (K_inf > 0) init_list$beta_inf <- rep(0.0, K_inf)

      # CONDITIONAL INITIALIZATION
      if (stan_data$use_curve_logic == 1) {
        init_list$gen_shape_raw <- array(3.0, dim = 1)
        init_list$gen_rate_raw <- array(1.0, dim = 1)
      }

      if (stan_data$use_vl_data == 1 && stan_data$vl_type == 0) {
        init_list$Ct50_raw <- array(35.0, dim = 1)
        init_list$slope_ct_raw <- array(2.0, dim = 1)
      }

      if (stan_data$use_vl_data == 1 && stan_data$vl_type == 1) {
        init_list$V_ref_raw <- array(3.0, dim = 1)
        init_list$rho_raw <- array(2.5, dim = 1)
      }

      return(init_list)
    }
  }

  # Run sampler
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
