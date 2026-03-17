
# Suppress R CMD CHECK "no visible binding" NOTEs for dplyr NSE column names.
utils::globalVariables(c(
  "susceptible_at_start", "first_inf_time", "infector_id",
  "has_index", "n_index", "is_index", "is_secondary",
  "n_total", "n_index_role", "n_suscept", "n_secondary",
  "sar_role", "sar", "sar_sd", "sar_se"
))


#' Calculate Secondary Attack Rate (Robust)
#'
#' Computes household secondary attack rates (SAR) from simulation output,
#' with fixes for generation conflation, co-index ambiguity, inconsistent
#' pooling, community importation, infinite first-infection time, and
#' susceptibility at study entry.
#'
#' @param sim_result A named list containing at minimum a \code{hh_df} data
#'   frame with household simulation output. Required columns: \code{hh_id},
#'   \code{person_id}, \code{role}, \code{infection_time}, \code{infector_id},
#'   \code{susceptible_at_start}.
#' @param generation Character string controlling which infections count toward
#'   the SAR numerator. One of:
#'   \describe{
#'     \item{\code{"strict_secondary"}}{(default) Only generation-2 infections
#'       (i.e., directly infected by the index case).}
#'     \item{\code{"all_subsequent"}}{All non-index infections, regardless of
#'       generation (original behavior).}
#'   }
#' @param pooling Character string controlling how household-level SARs are
#'   aggregated. One of:
#'   \describe{
#'     \item{\code{"mean_of_hh"}}{(default) Each household receives equal
#'       weight (recommended).}
#'     \item{\code{"pooled"}}{Individual-level pooling; larger households
#'       receive more weight.}
#'   }
#' @param alpha_flag Logical. If \code{TRUE} (default), emits a warning when
#'   post-index infections with no recorded within-household infector are
#'   detected, suggesting possible community importation.
#' @param verbose Logical. If \code{TRUE} (default), prints diagnostic
#'   messages via \code{rlang::inform()}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{overall}}{Numeric scalar: the overall household SAR.}
#'     \item{\code{by_age}}{Data frame of role-specific SARs, including
#'       \code{n_households}, \code{sar}, \code{sar_sd}, \code{sar_se},
#'       \code{sar_lo}, \code{sar_hi}, \code{n_secondary}, \code{n_suscept},
#'       and \code{overall_sar}.}
#'   }
#'
#' @examples
#' \dontrun{
#'   res <- calculate_sar_robust(
#'     sim_result = my_sim,
#'     generation = "strict_secondary",
#'     pooling    = "mean_of_hh",
#'     alpha_flag = TRUE,
#'     verbose    = TRUE
#'   )
#'   res$overall
#'   res$by_age
#' }
#'
#' @importFrom dplyr filter group_by summarize mutate left_join n_distinct
#'   select case_when n
#' @importFrom rlang abort warn inform
#' @importFrom stats sd
#' @importFrom utils globalVariables
#' @export
calculate_sar_robust <- function(sim_result,
                                 generation = c("strict_secondary", "all_subsequent"),
                                 pooling    = c("mean_of_hh", "pooled"),
                                 alpha_flag = TRUE,
                                 verbose    = TRUE) {

  # ------------------------------------------------------------------
  # INPUT VALIDATION
  # ------------------------------------------------------------------
  if (!is.list(sim_result) || !"hh_df" %in% names(sim_result))
    rlang::abort("`sim_result` must be a named list containing a `hh_df` element.")
  if (!is.data.frame(sim_result$hh_df))
    rlang::abort("`sim_result$hh_df` must be a data frame.")
  if (!is.logical(alpha_flag) || length(alpha_flag) != 1)
    rlang::abort("`alpha_flag` must be a single logical value (TRUE or FALSE).")
  if (!is.logical(verbose) || length(verbose) != 1)
    rlang::abort("`verbose` must be a single logical value (TRUE or FALSE).")

  generation <- match.arg(generation)
  pooling    <- match.arg(pooling)

  hh_df <- sim_result$hh_df

  # ------------------------------------------------------------------
  # GUARD: required columns
  # ------------------------------------------------------------------
  required_cols <- c("hh_id", "person_id", "role", "infection_time",
                     "infector_id", "susceptible_at_start")
  missing_cols  <- setdiff(required_cols, names(hh_df))

  if (length(missing_cols) > 0) {
    rlang::warn(sprintf(
      "Missing columns: [%s]. Falling back to defaults (may affect accuracy).",
      paste(missing_cols, collapse = ", ")
    ))
    if (!"infector_id"          %in% names(hh_df)) hh_df$infector_id          <- NA
    if (!"susceptible_at_start" %in% names(hh_df)) hh_df$susceptible_at_start <- TRUE
  }

  # ------------------------------------------------------------------
  # FIX 6: Susceptibility at enrollment
  # ------------------------------------------------------------------
  n_immune_entry <- sum(!hh_df$susceptible_at_start, na.rm = TRUE)
  if (verbose && n_immune_entry > 0)
    rlang::inform(sprintf(
      "[SAR] Excluding %d individual(s) immune at study entry.", n_immune_entry
    ))

  hh_df <- dplyr::filter(hh_df, susceptible_at_start == TRUE | is.na(susceptible_at_start))

  # ------------------------------------------------------------------
  # STEP 1: Per-household index identification
  # FIX 5: Zero-infection households -> first_inf_time = Inf
  # FIX 2: Co-index ties flagged
  # ------------------------------------------------------------------
  hh_summary <- hh_df %>%
    dplyr::group_by(hh_id) %>%
    dplyr::summarize(
      household_size = dplyr::n_distinct(person_id),
      n_infected     = sum(!is.na(infection_time)),
      first_inf_time = {
        inf_times <- infection_time[!is.na(infection_time)]
        if (length(inf_times) == 0) Inf else min(inf_times)
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      has_index = !is.infinite(first_inf_time),
      n_index   = mapply(function(hid, ft) {
        if (is.infinite(ft)) return(0L)
        sum(
          hh_df$hh_id == hid &
            !is.na(hh_df$infection_time) &
            hh_df$infection_time == ft
        )
      }, hh_id, first_inf_time)
    )

  n_ties <- sum(hh_summary$n_index > 1, na.rm = TRUE)
  if (verbose && n_ties > 0)
    rlang::inform(sprintf(
      "[SAR] %d household(s) have co-index ties (simultaneous first infection).", n_ties
    ))

  # ------------------------------------------------------------------
  # FIX 4: Community importation flag
  # ------------------------------------------------------------------
  if (alpha_flag && "infector_id" %in% names(hh_df)) {
    post_index_community <- hh_df %>%
      dplyr::left_join(
        dplyr::select(hh_summary, hh_id, first_inf_time), by = "hh_id"
      ) %>%
      dplyr::filter(
        !is.na(infection_time),
        infection_time > first_inf_time,
        is.na(infector_id)
      ) %>%
      nrow()

    if (verbose && post_index_community > 0)
      rlang::inform(sprintf(
        "[SAR] %d post-index infection(s) with no recorded HH infector (possible community import).",
        post_index_community
      ))
  }

  # ------------------------------------------------------------------
  # FIX 1: Generation tracking
  # ------------------------------------------------------------------
  hh_df <- hh_df %>%
    dplyr::left_join(
      dplyr::select(hh_summary, hh_id, first_inf_time, n_index, has_index),
      by = "hh_id"
    ) %>%
    dplyr::mutate(
      is_index = !is.na(infection_time) &
        has_index &
        infection_time == first_inf_time,

      is_secondary = dplyr::case_when(
        generation == "strict_secondary" ~
          !is.na(infection_time) &
          !is_index &
          (infector_id %in% person_id[is_index]),
        generation == "all_subsequent" ~
          !is.na(infection_time) & !is_index,
        TRUE ~
          !is.na(infection_time) & !is_index
      )
    )

  if (generation == "strict_secondary" && all(is.na(hh_df$infector_id))) {
    rlang::warn(
      "[SAR] `infector_id` is all NA -- `strict_secondary` falls back to `all_subsequent`."
    )
    hh_df <- dplyr::mutate(hh_df, is_secondary = !is.na(infection_time) & !is_index)
  }

  # ------------------------------------------------------------------
  # FIX 3: Consistent pooling for overall and age-specific SAR
  # ------------------------------------------------------------------

  # Household-level SAR per role
  hh_role_sar <- hh_df %>%
    dplyr::filter(has_index) %>%
    dplyr::group_by(hh_id, role) %>%
    dplyr::summarize(
      n_total      = dplyr::n(),
      n_index_role = sum(is_index,     na.rm = TRUE),
      n_secondary  = sum(is_secondary, na.rm = TRUE),
      n_suscept    = n_total - n_index_role,
      sar_role     = ifelse(n_suscept > 0, n_secondary / n_suscept, NA_real_),
      .groups = "drop"
    )

  # Household-level overall SAR
  hh_overall_sar <- hh_df %>%
    dplyr::filter(has_index) %>%
    dplyr::group_by(hh_id) %>%
    dplyr::summarize(
      n_index     = sum(is_index,     na.rm = TRUE),
      n_secondary = sum(is_secondary, na.rm = TRUE),
      n_suscept   = dplyr::n() - sum(is_index, na.rm = TRUE),
      sar_hh      = ifelse(n_suscept > 0, n_secondary / n_suscept, NA_real_),
      .groups = "drop"
    )

  # Aggregate
  if (pooling == "mean_of_hh") {

    overall_sar <- mean(hh_overall_sar$sar_hh, na.rm = TRUE)

    age_sar <- hh_role_sar %>%
      dplyr::group_by(role) %>%
      dplyr::summarize(
        n_households = dplyr::n(),
        sar          = mean(sar_role,        na.rm = TRUE),
        sar_sd       = stats::sd(sar_role,   na.rm = TRUE),
        sar_se       = sar_sd / sqrt(sum(!is.na(sar_role))),
        sar_lo       = sar - 1.96 * sar_se,
        sar_hi       = sar + 1.96 * sar_se,
        n_secondary  = sum(n_secondary),
        n_suscept    = sum(n_suscept),
        .groups = "drop"
      )

  } else {  # pooled

    overall_sar <- sum(hh_overall_sar$n_secondary, na.rm = TRUE) /
      sum(hh_overall_sar$n_suscept,   na.rm = TRUE)

    age_sar <- hh_role_sar %>%
      dplyr::group_by(role) %>%
      dplyr::summarize(
        n_households = dplyr::n(),
        n_secondary  = sum(n_secondary),
        n_suscept    = sum(n_suscept),
        sar          = ifelse(n_suscept > 0, n_secondary / n_suscept, NA_real_),
        .groups = "drop"
      )
  }

  if (is.na(overall_sar) || is.nan(overall_sar)) overall_sar <- 0
  age_sar$overall_sar <- overall_sar

  # ------------------------------------------------------------------
  # DIAGNOSTICS
  # ------------------------------------------------------------------
  if (verbose) {
    n_hh_total    <- dplyr::n_distinct(hh_df$hh_id)
    n_hh_with_idx <- sum(hh_summary$has_index)
    n_hh_no_idx   <- n_hh_total - n_hh_with_idx

    rlang::inform(sprintf(
      "[SAR] %d households | %d with index | %d no infection | %d co-index tie(s)",
      n_hh_total, n_hh_with_idx, n_hh_no_idx, n_ties
    ))
    rlang::inform(sprintf(
      "[SAR] generation='%s' | pooling='%s' | overall SAR=%.3f",
      generation, pooling, overall_sar
    ))
  }

  list(overall = overall_sar, by_age = age_sar)
}


# ============================================================================
# calculate_sar_quick() -- convenience wrapper
# ============================================================================

#' Quick SAR Calculator (Convenience Wrapper)
#'
#' A drop-in replacement for a legacy \code{calculate_sar_quick()} function.
#' Calls \code{\link{calculate_sar_robust}} with recommended defaults:
#' \code{generation = "strict_secondary"} and \code{pooling = "mean_of_hh"}.
#'
#' @param sim_result A named list containing a \code{hh_df} data frame.
#'   See \code{\link{calculate_sar_robust}} for details.
#'
#' @return See \code{\link{calculate_sar_robust}}.
#'
#' @seealso \code{\link{calculate_sar_robust}}
#' @export
calculate_sar_quick <- function(sim_result) {
  calculate_sar_robust(
    sim_result = sim_result,
    generation = "strict_secondary",
    pooling    = "mean_of_hh",
    alpha_flag = TRUE,
    verbose    = TRUE
  )
}
