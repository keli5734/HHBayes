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
      list(
        log_beta1 = -5.3,  # Approx log(0.005)
        log_beta2 = -5.3,

        # Initialize multipliers near 1.0 (log scale 0) to start neutral
        # We check the size of R to ensure the vector length is correct
        log_phi_by_role_raw = rep(0.1, max(1, stan_data$R - 1)),
        log_kappa_by_role_raw = rep(0.1, max(1, stan_data$R - 1)),

        # Viral Load parameters (if used)
        V_ref = 3.0,
        V_rho = 2.5,

        # Community infection rate (if estimated)
        log_beta3 = 0, # Placeholder for antibody effects if active
        log_beta4 = 0
      )
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
