#' Summarize Attack Rates & Reinfections
#'
#' Calculates the Primary Attack Rate (proportion of people infected at least once)
#' and a separate summary of Reinfections (secondary episodes).
#'
#' @param sim_result A list object returned by \code{simulate_multiple_households_comm}.
#'
#' @return A list containing four dataframes:
#' \describe{
#'   \item{primary_overall}{Overall primary attack rate statistics.}
#'   \item{primary_by_role}{Primary attack rate stratified by age group.}
#'   \item{reinf_overall}{Overall summary of reinfection counts and rates.}
#'   \item{reinf_by_role}{Reinfection counts stratified by age group.}
#' }
#' @export
summarize_attack_rates <- function(sim_result) {

  if (!"hh_df" %in% names(sim_result)) stop("Input must contain 'hh_df'.")
  df <- sim_result$hh_df

  # ---------------------------------------------------------
  # 1. PREPARE DATA
  # ---------------------------------------------------------

  # A. Total Population (Denominator)
  pop_stats <- df %>%
    dplyr::select(hh_id, person_id, role) %>%
    dplyr::distinct() %>%
    dplyr::group_by(role) %>%
    dplyr::summarise(N_pop = dplyr::n(), .groups="drop")

  total_pop <- sum(pop_stats$N_pop)

  # B. Infection Episodes (Numerator)
  # Filter only rows representing an actual infection
  infected_df <- df %>%
    dplyr::filter(!is.na(infection_time)) %>%
    dplyr::select(hh_id, person_id, role, infection_time)

  # ---------------------------------------------------------
  # 2. PRIMARY ATTACK RATE (First Infection Only)
  # ---------------------------------------------------------

  # Count UNIQUE people infected (at least once)
  primary_counts_by_role <- infected_df %>%
    dplyr::group_by(role) %>%
    dplyr::summarise(n_primaries = dplyr::n_distinct(paste(hh_id, person_id)), .groups="drop")

  # Merge with Population
  primary_by_role <- pop_stats %>%
    dplyr::left_join(primary_counts_by_role, by = "role") %>%
    dplyr::mutate(
      n_primaries = ifelse(is.na(n_primaries), 0, n_primaries),
      Primary_AR  = n_primaries / N_pop,
      Primary_AR_Pct = scales::percent(Primary_AR, accuracy = 0.1)
    )

  # Overall Primary Stats
  total_primaries <- sum(primary_by_role$n_primaries)
  primary_overall <- data.frame(
    Total_Pop = total_pop,
    Total_Infected_People = total_primaries,
    Primary_Attack_Rate = total_primaries / total_pop,
    Primary_Attack_Rate_Pct = scales::percent(total_primaries / total_pop, accuracy = 0.1)
  )

  # ---------------------------------------------------------
  # 3. REINFECTION SUMMARY (Subsequent Infections)
  # ---------------------------------------------------------

  # Calculate Total Episodes vs Unique People
  reinf_counts_by_role <- infected_df %>%
    dplyr::group_by(role) %>%
    dplyr::summarise(
      n_total_episodes = dplyr::n(),
      n_unique_people  = dplyr::n_distinct(paste(hh_id, person_id)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_reinfections = n_total_episodes - n_unique_people,
      # Reinfection Rate: % of infected people who got it AGAIN
      Reinf_Rate_Among_Infected = ifelse(n_unique_people > 0, n_reinfections / n_unique_people, 0)
    )

  # Merge to ensure all roles appear even if 0 reinfections
  reinf_by_role <- pop_stats %>%
    dplyr::select(role) %>%
    dplyr::left_join(reinf_counts_by_role, by = "role") %>%
    dplyr::mutate(
      n_reinfections = ifelse(is.na(n_reinfections), 0, n_reinfections),
      Reinf_Rate_Pct = scales::percent(ifelse(is.na(Reinf_Rate_Among_Infected), 0, Reinf_Rate_Among_Infected), accuracy = 0.1)
    ) %>%
    dplyr::select(role, n_unique_people, n_total_episodes, n_reinfections, Reinf_Rate_Pct)

  # Overall Reinfection Stats
  total_episodes_all <- nrow(infected_df)
  total_reinf_all    <- total_episodes_all - total_primaries

  reinf_overall <- data.frame(
    Total_Primary_Infections = total_primaries,
    Total_Reinfection_Events = total_reinf_all,
    Reinfection_Rate_Among_Infected = ifelse(total_primaries > 0, total_reinf_all / total_primaries, 0)
  )

  # ---------------------------------------------------------
  # 4. RETURN
  # ---------------------------------------------------------

  list(
    primary_overall = primary_overall,
    primary_by_role = primary_by_role,
    reinf_overall   = reinf_overall,
    reinf_by_role   = reinf_by_role
  )
}

#' Plot Dual-Axis Epidemic Curve (Binned)
#'
#' Overlays a stacked bar chart of simulated infections with a line chart of surveillance data.
#' Both datasets are aggregated by 'bin_width' days.
#'
#' @param sim_result Output object from `simulate_multiple_households_comm`.
#' @param surveillance_df Dataframe with `date` and `cases` columns.
#' @param start_date_str String. Start date of the simulation.
#' @param bin_width Integer. Aggregation window in days (default 7).
#' @return A ggplot object.
#' @import ggplot2 dplyr
#' @export
plot_epidemic_curve <- function(sim_result, surveillance_df, start_date_str = "2024-07-01", bin_width = 7) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  start_date <- as.Date(start_date_str)

  # --- 1. PREPARE SIMULATION DATA (BINNED) ---
  df_sim_raw <- sim_result$hh_df %>%
    dplyr::filter(!is.na(infection_time))

  if(nrow(df_sim_raw) == 0) return(NULL)

  df_sim_binned <- df_sim_raw %>%
    dplyr::mutate(
      # Calculate Bin Start Date
      day_idx = infection_time - 1,
      bin_start_idx = floor(day_idx / bin_width) * bin_width,
      date_formatted = start_date + bin_start_idx,

      # Role Mapping
      role_group = dplyr::case_when(
        role == "infant" ~ "Infant",
        role == "toddler" ~ "Toddler",
        role == "adult" ~ "Adult",
        role == "elderly" ~ "Elderly",
        TRUE ~ "Other"
      ),
      role_group = factor(role_group, levels = c("Infant", "Toddler", "Adult", "Elderly"))
    ) %>%
    # Aggregate counts per bin per role
    dplyr::count(date_formatted, role_group, name = "n_infections")

  # --- 2. PREPARE SURVEILLANCE DATA (BINNED) ---
  # We must aggregate surveillance to the same bins for a fair comparison
  surv_binned <- surveillance_df %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::filter(date >= start_date) %>%
    dplyr::mutate(
      day_idx = as.numeric(date - start_date),
      # Only keep positive indices (within study period)
      valid = day_idx >= 0
    ) %>%
    dplyr::filter(valid) %>%
    dplyr::mutate(
      bin_start_idx = floor(day_idx / bin_width) * bin_width,
      date_bin = start_date + bin_start_idx
    ) %>%
    dplyr::group_by(date_bin) %>%
    dplyr::summarise(cases = sum(cases, na.rm = TRUE)) %>%
    dplyr::ungroup()

  # Trim surveillance to simulation range for cleaner plotting
  max_sim_date <- max(df_sim_binned$date_formatted)
  surv_binned <- surv_binned %>% dplyr::filter(date_bin <= max_sim_date + bin_width)

  # --- 3. CALCULATE SCALING ---
  # Max height of stacked bars
  max_bar_val <- df_sim_binned %>%
    dplyr::group_by(date_formatted) %>%
    dplyr::summarise(total = sum(n_infections)) %>%
    dplyr::pull(total) %>%
    max(na.rm = TRUE)

  if(is.infinite(max_bar_val) || max_bar_val == 0) max_bar_val <- 1

  # Max height of line
  max_line_val <- if(nrow(surv_binned) > 0) max(surv_binned$cases, na.rm = TRUE) else 1

  coeff <- max_line_val / max_bar_val
  coeff <- coeff * 1.1 # Buffer

  # --- 4. PLOT ---
  custom_colors <- c(
    "Adult"   = "#00A1D5FF",
    "Infant"  = "#79AF97FF",
    "Toddler" = "#DF8F44FF",
    "Elderly" = "#B24745FF"
  )

  ggplot2::ggplot() +
    # Layer 1: Surveillance Line (Scaled Down)
    # We use geom_step or geom_line depending on preference. Line is usually safer for trends.
    ggplot2::geom_line(
      data = surv_binned,
      ggplot2::aes(x = date_bin, y = cases / coeff, color = "Surveillance"),
      linewidth = 1, alpha = 1
    ) +

    # Layer 2: Simulated Bars (Primary Axis)
    ggplot2::geom_bar(
      data = df_sim_binned,
      ggplot2::aes(x = date_formatted, y = n_infections, fill = role_group),
      stat = "identity",
      width = bin_width - 0.5, # Slightly smaller than bin for gap
      color = "transparent",
      alpha = .6
    ) +

    # Scales
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::scale_color_manual(name = "", values = c("Surveillance" = "black")) +

    ggplot2::scale_y_continuous(
      name = "# of Positive Samples",
      breaks = scales::pretty_breaks(n = 6),
      sec.axis = ggplot2::sec_axis(~ . * coeff, name = "Surveillance Cases")
    ) +

    ggplot2::scale_x_date(date_labels = "%b %d", date_breaks = "4 weeks") +

    ggplot2::labs(x = "Date", fill = "Age Group") +

    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.ticks = ggplot2::element_line(linewidth = 1),
      axis.text.x = ggplot2::element_text(angle = 0, color = "black", size = 18),
      axis.text.y = ggplot2::element_text(color = "black", size = 18),
      axis.title.y.left = ggplot2::element_text(face = "bold", margin = ggplot2::margin(r = 10), size = 18),
      axis.title.y.right = ggplot2::element_text(face = "bold", angle = 90, margin = ggplot2::margin(l = 10), size = 18),
      axis.title.x = ggplot2::element_text(face = "bold", angle = 0, size = 18),
      panel.grid.major = ggplot2::element_blank(),
      legend.position = "top"
    )
}

