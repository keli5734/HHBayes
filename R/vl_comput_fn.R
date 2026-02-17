#' Compute Theoretical Ct Trajectory (Piecewise Linear)
#' @keywords internal
get_theoretical_ct <- function(t, p) {
  # t = days since infection
  # p = list(Cpeak, r, d, t_peak)
  ifelse(t <= p$t_peak,
         p$Cpeak + p$r * (p$t_peak - t),
         p$Cpeak + p$d * (t - p$t_peak))
}

#' Compute Theoretical Log10 Viral Load Trajectory (Double Exponential)
#' @keywords internal
get_theoretical_vl <- function(t, p) {
  # t = days since infection
  # p = list(v_p, t_p, lambda_g, lambda_d)
  numerator <- 2 * 10^(p$v_p)
  denominator <- exp(-p$lambda_g * (t - p$t_p)) + exp(p$lambda_d * (t - p$t_p))
  vt <- numerator / denominator
  return(log10(pmax(1e-9, vt)))
}

#' Fill Missing Viral Data (Ct or Log10) based on Episode Start
#'
#' This function imputes missing viral data during confirmed episodes using
#' theoretical trajectories defined by the parameters.
#'
#' @param df Dataframe containing 'person_id', 'episode_id', 'date', 'role_name', and the viral column.
#' @param viral_col_name String. Name of the column containing viral data (e.g. "ct_value").
#' @param type String. Either "ct_value" or "log10".
#' @param params_list List of parameters for the curves (role-specific).
#' @param detection_limit Numeric. Value to assign for non-infected days (e.g. 45 for Ct, 0 for Log10).
#' @return The original dataframe with NAs in the viral column filled.
#' @export
fill_missing_viral_data <- function(df, viral_col_name, type = c("ct_value", "log10"),
                                    params_list, detection_limit) {

  type <- match.arg(type)

  # Ensure the viral column exists
  if(!viral_col_name %in% names(df)) stop(paste("Column", viral_col_name, "not found in dataframe."))

  # Prepare for non-standard evaluation
  viral_sym <- rlang::sym(viral_col_name)

  df_imputed <- df %>%
    dplyr::group_by(person_id, episode_id) %>%
    dplyr::mutate(
      # 1. Identify start of episode
      episode_start = min(date, na.rm = TRUE),

      # 2. Calculate days since infection (0-indexed)
      days_since_inf = as.numeric(date - episode_start),

      # 3. Impute Logic
      imputed_val = dplyr::case_when(
        !is.na(!!viral_sym) ~ as.numeric(!!viral_sym),  # Keep existing data
        episode_id == 0     ~ detection_limit,          # No infection
        TRUE ~ {
          # Look up params
          role <- unique(role_name)[1]
          if(is.na(role)) role <- "adult"
          p <- params_list[[role]]
          if(is.null(p)) p <- params_list[["adult"]] # Fallback

          # Calculate theoretical value
          if(type == "ct_value") {
            theo <- get_theoretical_ct(days_since_inf, p)
            # Clamp Ct to valid range (e.g., 15 to detection limit)
            pmin(detection_limit, pmax(10, theo))
          } else {
            theo <- get_theoretical_vl(days_since_inf, p)
            # Clamp Log10 to valid range (e.g., detection limit to 12)
            pmax(detection_limit, theo)
          }
        }
      )
    ) %>%
    dplyr::ungroup() %>%
    # Overwrite the original column
    dplyr::mutate(!!viral_col_name := imputed_val) %>%
    dplyr::select(-episode_start, -days_since_inf, -imputed_val)

  return(df_imputed)
}
