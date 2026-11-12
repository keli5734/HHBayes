#' @details
#' ## Funding
#' This work was supported by a grant from the National Institutes of Health
#' (R01AI137093). The content is solely the responsibility of the authors and
#' does not necessarily represent the official views of the National Institutes
#' of Health.
#' @keywords internal
"_PACKAGE"

# Package-level imports. Kept HERE (not in R/stanmodels.R) because
# rstantools::rstan_config() regenerates stanmodels.R; without useDynLib the
# compiled Stan models are never loaded ("no such symbol _rcpp_module_boot_...").
#' @useDynLib HHBayes, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
NULL
