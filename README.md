# HHBayes <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/keli5734/HHBayes)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

> **Bayesian Household Transmission Modeling in R**

**HHBayes** is an R package for simulating and estimating infectious disease transmission within households. It couples a stochastic, age-structured household simulator with Bayesian inference via [Stan](https://mc-stan.org/), enabling researchers to study within-household spread, evaluate intervention strategies, and reconstruct transmission chains from longitudinal diagnostic data.

---

## Why HHBayes?

Household transmission studies generate rich but complex data — repeated tests on every member, overlapping infections, reinfections, and age-dependent viral kinetics. HHBayes provides an integrated pipeline:

1. **Simulate** realistic outbreaks with configurable household structures, contact patterns, and interventions.
2. **Prepare** observed or simulated data for Bayesian analysis, with automatic imputation of infection/recovery windows and viral load curves.
3. **Fit** a Stan model that jointly estimates transmission rates, role-specific susceptibility/infectivity, covariate effects, and viral load scaling.
4. **Analyze** results with built-in tools for attack rates, transmission chain reconstruction, and publication-ready plots.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("keli5734/HHBayes")
```

### Prerequisites

HHBayes depends on [RStan](https://mc-stan.org/rstan/), which requires a working C++ toolchain. See the [RStan Getting Started Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for platform-specific setup.

---

## Quick Start

```r
library(HHBayes)
library(rstan)
library(ggpubr)

# 0. Generate surveillance dataset
study_start <- "2024-07-01"
study_end   <- "2025-06-30"
dates_weekly <- seq(from = as.Date(study_start), to = as.Date(study_end), by = "week")
surveillance_data <- data.frame(
  date = dates_weekly,
  # Random epidemic curve (low start, peak middle, low end)
  cases = 0.1 + 100 * exp(-0.0002 * (as.numeric(dates_weekly - mean(dates_weekly)))^2) + abs(rnorm( length(dates_weekly),mean = 0, sd = 10))
)

# 1. Simulate 50 households with ODE-based viral dynamics
sim <- simulate_multiple_households_comm(
  n_households = 50,
  viral_testing = "viral load",
  start_date = "2024-07-01",
  end_date   = "2025-06-30",
  seed = 123,
  surveillance_df = surveillance_data
)

# 2. Summarize attack rates
rates <- summarize_attack_rates(sim)
rates$primary_by_role

# 3. Prepare data and fit the Bayesian model
stan_input <- prepare_stan_data(sim$diagnostic_df,
                                use_vl_data      = 1,
                                surveillance_df  = surveillance_data,
                                study_start_date = as.Date("2024-07-01"),
                                study_end_date   = as.Date("2025-06-30"),
                                seed = 123)
  
options(mc.cores = parallel::detectCores())

fit <- fit_household_model(stan_input, iter = 2000, chains = 4)

# 4. Visualize
plot_posterior_distributions(fit)
```

---

## Full Pipeline Walkthrough

The sections below walk through a complete analysis: from defining a study scenario, through simulation, data preparation, Bayesian model fitting, and visualization.

### Step 1: Define the Study Timeline and Surveillance Data

The simulation is anchored to a calendar period. You can optionally supply external surveillance data (a dataframe with `date` and `cases` columns) to drive time-varying community importation. The case counts are automatically normalized and interpolated to daily resolution.

```r
study_start <- "2024-07-01"
study_end   <- "2025-06-30"

# Weekly surveillance data with a mid-study epidemic peak
dates_weekly <- seq(as.Date(study_start), as.Date(study_end), by = "week")
surveillance_data <- data.frame(
  date  = dates_weekly,
  cases = 0.1 + 100 * exp(-0.0002 * (as.numeric(dates_weekly - mean(dates_weekly)))^2) +
    abs(rnorm(length(dates_weekly), 0, 10))
)
```

If no surveillance data is provided, you can supply custom `seasonal_forcing_list` vectors or leave it flat (constant community force).

### Step 2: Define the Contact Matrix

By default, all household members contact each other equally. A `role_mixing_matrix` lets you specify **differential contact weights** between roles. Each cell represents the relative contact intensity between a source (row) and target (column):

```r
role_mixing_weights <- matrix(c(
# Target: Infant Toddler Adult Elderly
          0.0,   0.5,    1.0,  0.5,    # Source: Infant  (high contact with adults/caregivers)
          0.5,   0.9,    0.7,  0.5,    # Source: Toddler (high peer contact)
          1.0,   0.7,    0.6,  0.7,    # Source: Adult   (caregiver role)
          0.5,   0.5,    0.7,  0.0     # Source: Elderly (limited infant contact)
), nrow = 4, byrow = TRUE,
   dimnames = list(
     c("infant", "toddler", "adult", "elderly"),
     c("infant", "toddler", "adult", "elderly")))
```

This matrix is expanded internally to an N x N individual-level contact matrix for each household, based on the role of each member. For example, in a household with roles `c("adult", "adult", "infant")`, the adult-to-infant weight (1.0) is applied to those pairs.

### Step 3: Define Household Structure

Each simulated household is assembled randomly according to a demographic profile. The `household_profile_list` controls the probability of different compositions:

```r
household_profile <- list(
  prob_adults   = c(0, 0, 1),        # 0% chance of 0 or one adult, 100% chance of two
  prob_infant   = 1.0,               # 100% chance of one infant
  prob_siblings = c(0, 0.8, 0.2),    # 80% chance of one toddler, 20% chance of two
  prob_elderly  = c(0.7, 0.1, 0.2)   # 70% no elderly, 10% one, 20% two
)
```

| Field | Values | Meaning |
|---|---|---|
| `prob_adults` | `c(P0, P1, P2)` | Probability of 0, 1, or 2 adults |
| `prob_infant` | `0` to `1` | Probability of an infant being present |
| `prob_siblings` | `c(P0, P1, P2)` | Probability of 0, 1, or 2 toddler siblings |
| `prob_elderly` | `c(P0, P1, P2)` | Probability of 0, 1, or 2 elderly members |

### Step 4: Define Intervention Strategies

Interventions (vaccination, masking, prophylaxis, etc.) are specified as a list via `covariates_config`. Each entry defines an intervention with a name, what it affects, its efficacy, and role-specific coverage probabilities:

```r
sim_config <- list(
  list(
    name      = "vacc_status",
    efficacy  = 0.5,                # 50% reduction
    effect_on = "both",             # "susceptibility", "infectivity", or "both"
    coverage  = list(
      infant  = 0.00,
      toddler = 0.30,
      adult   = 0.80,
      elderly = 0.90
    )
  )
)
```

**How efficacy works:** For each person, a Bernoulli draw (based on `coverage`) determines if they receive the intervention. If they do, their susceptibility and/or infectivity is multiplied by `(1 - efficacy)`. Multiple interventions stack multiplicatively.

Setting `efficacy = 0` (as in the baseline example) creates the covariate column without any actual effect — useful for testing the pipeline.

### Step 5: Run the Simulation

```r
sim_res <- simulate_multiple_households_comm(
  n_households    = 50,
  viral_testing   = "viral load",    # or "Ct"
  infectious_shape = 10,
  infectious_scale = 1,              # Mean infectious period = 10 days
  waning_shape    = 6,
  waning_scale    = 10,              # Mean immunity = 60 days
  surveillance_interval = 4,         # Test every 4 days
  start_date      = study_start,
  end_date        = study_end,
  surveillance_df = surveillance_data,
  covariates_config      = sim_config,
  household_profile_list = household_profile,
  role_mixing_matrix     = role_mixing_weights,
  seed = 123
)
```

The output is a list with two dataframes:

| Output | Description |
|---|---|
| `sim_res$hh_df` | One row per person per infection episode: household ID, person ID, role, infection time, recovery times, and covariate assignments. |
| `sim_res$diagnostic_df` | Simulated test results at each surveillance time point: viral load or Ct value, binary test result, and episode ID. |

### Step 6: Summarize Attack Rates

```r
rates <- summarize_attack_rates(sim_res)

rates$primary_overall      # Overall primary attack rate
rates$primary_by_role      # Attack rate by age group
rates$reinf_overall        # Overall reinfection summary
rates$reinf_by_role        # Reinfections by age group
```

### Step 7: Plot the Epidemic Curve

Overlays simulated infections (stacked bars by age group) with the surveillance line:

```r
plot_epidemic_curve(sim_res, surveillance_data,
                    start_date_str = study_start, bin_width = 7)
```

### Step 8: Prepare Data for Stan

This step transforms the diagnostic dataframe into the structured list that Stan expects. It handles column renaming, infection window imputation, viral load gap-filling, covariate matrix construction, and prior specification.

**Important:** Covariates from the simulation live in `hh_df`, but Stan needs them joined to `diagnostic_df`:

```r
# Extract 1-row-per-person covariate table
person_covariates <- sim_res$hh_df %>%
  dplyr::select(hh_id, person_id, vacc_status) %>%
  dplyr::distinct()

# Join to diagnostic data
df_for_stan <- sim_res$diagnostic_df %>%
  dplyr::left_join(person_covariates, by = c("hh_id", "person_id"))
```

**Define priors** for the Bayesian model. Supported distributions are `"normal"`, `"uniform"`, and `"lognormal"`:

```r
my_priors <- list(
  beta1     = list(dist = "normal",    params = c(-5, 1)),    # log-transmission rate 1
  beta2     = list(dist = "normal",    params = c(-5, 1)),    # log-transmission rate 2
  alpha     = list(dist = "normal",    params = c(-4, 1)),    # log-community rate
  covariates = list(dist = "normal",   params = c(0, 2)),     # covariate coefficients
  gen_shape = list(dist = "lognormal", params = c(1.5, 0.5)), # generation interval shape
  gen_rate  = list(dist = "lognormal", params = c(0.0, 0.5)), # generation interval rate
  ct50      = list(dist = "normal",    params = c(3, 1)),     # VL reference point
  slope     = list(dist = "lognormal", params = c(0.4, 0.5))  # VL dose-response slope
)
```

**Supply viral load imputation parameters** (used to fill gaps in observed viral data):

```r
VL_params_list <- list(
  adult   = list(v_p = 4.14, t_p = 5.09, lambda_g = 2.31, lambda_d = 2.71),
  infant  = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  toddler = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  elderly = list(v_p = 2.95, t_p = 5.10, lambda_g = 3.15, lambda_d = 0.87)
)
```

**Build the Stan input:**

```r
stan_input <- prepare_stan_data(
  df_clean         = df_for_stan,
  surveillance_df  = surveillance_data,
  study_start_date = as.Date(study_start),
  study_end_date   = as.Date(study_end),
  use_vl_data      = 1,
  imputation_params = VL_params_list,
  priors           = my_priors,
  role_mixing_matrix = role_mixing_weights
)
```

### Step 9: Fit the Bayesian Model

`fit_household_model()` runs HMC sampling via `rstan::sampling()`. You can exclude internal bookkeeping parameters to keep the output clean:

```r
options(mc.cores = parallel::detectCores())

fit <- fit_household_model(
  stan_input,
  pars    = c("log_phi_by_role_raw", "log_kappa_by_role_raw",
              "log_beta1", "log_beta2", "log_alpha_comm",
              "g_curve_est", "V_term_calc"),
  include = FALSE,       # Exclude these parameters from output
  iter    = 2000,
  warmup  = 1000,
  chains  = 4
)

print(fit, probs = c(0.025, 0.5, 0.975))
```

### Step 10: Visualize Results

```r
# Posterior distributions of susceptibility (phi) and infectivity (kappa) by role
plot_posterior_distributions(fit)

# Reconstruct who-infected-whom transmission chains
chains <- reconstruct_transmission_chains(fit, stan_input, min_prob_threshold = 0.001)

# Plot a specific household's infection timeline with transmission arrows
plot_household_timeline(chains, stan_input, target_hh_id = 11)
```

---

## Simulation Parameter Reference

### Transmission

| Parameter | Default | Description |
|---|---|---|
| `beta1` | `0.008` | Baseline within-household transmission rate (contact-driven). |
| `beta2` | `0.008` | Viral-load-dependent transmission rate. Total per-infector force = `beta1 + beta2 * f(VL)`. |
| `delta` | `0` | Household size scaling. Force scaled by `(1/max(n_h,1))^delta` in larger households. |
| `alpha_comm_by_role` | `5e-4` | Daily community acquisition rate. |
| `phi_by_role` | `c(adult=1, infant=3, toddler=4, elderly=1)` | Susceptibility multipliers by role. |
| `kappa_by_role` | `c(adult=1, infant=1, toddler=1.2, elderly=1)` | Infectivity multipliers by role. |

### Infectious Period and Immunity

| Parameter | Default | Description |
|---|---|---|
| `infectious_shape` | `3` | Gamma shape for infectious duration. |
| `infectious_scale` | `1` | Gamma scale for infectious duration. Mean = shape x scale. |
| `waning_shape` | `16` | Gamma shape for immunity waning. |
| `waning_scale` | `10` | Gamma scale for immunity waning. Mean = 160 days. |
| `max_infections` | `Inf` | Max infections per person. Set to `1` to disable reinfection. |

### Viral Dynamics

| Parameter | Default | Description |
|---|---|---|
| `viral_testing` | `"viral load"` | `"viral load"` (log10 scale) or `"Ct"` (cycle threshold). |
| `model_type` | `"empirical"` | `"empirical"` (parametric curves) or `"ODE"` (within-host ODE). |
| `V_ref` | `3.0` | Reference log10 VL for infectiousness scaling. |
| `V_rho` | `2.5` | Power-law exponent for VL scaling. |
| `Ct_50` | `40` | Ct at 50% infectiousness (sigmoid). |
| `Ct_delta` | `2` | Steepness of Ct-infectiousness sigmoid. |

### Testing and Surveillance

| Parameter | Default | Description |
|---|---|---|
| `surveillance_interval` | `1` | Days between routine tests. |
| `test_daily` | `FALSE` | Switch to daily testing after first household detection. |
| `perfect_detection` | `TRUE` | If `FALSE`, detection depends on viral load thresholds. |
| `detect_threshold_log10` | `1e-6` | Min log10 VL for positive test. |
| `detect_threshold_Ct` | `99` | Max Ct for positive test (set `35`-`40` for realistic PCR). |

---

## Key Functions

| Function | Description |
|---|---|
| `simulate_multiple_households_comm()` | Simulate multi-household outbreaks |
| `prepare_stan_data()` | Format data for Stan model |
| `fit_household_model()` | Fit the Bayesian model via HMC |
| `summarize_attack_rates()` | Attack rates and reinfection summaries |
| `reconstruct_transmission_chains()` | Posterior who-infected-whom probabilities |
| `fill_missing_viral_data()` | Impute missing Ct / viral load |
| `plot_posterior_distributions()` | Violin plots of phi and kappa |
| `plot_covariate_effects()` | Forest plot of intervention effects |
| `plot_epidemic_curve()` | Epidemic curve overlaid with surveillance |
| `plot_household_timeline()` | Per-household timeline with transmission arrows |

---

## Model Overview

The daily probability of infection for susceptible individual $i$ in household $h$ at time $t$ is:

$$P(\text{infection}_i(t)) = 1 - \exp\left(-\lambda_i(t)\right)$$

where the total force of infection combines community and household sources:

$$\lambda_{ih}(t) = \phi_{r_i} \cdot m_i \cdot \left( \alpha_{\text{comm}} \cdot S(t) + (\frac{1}{\max(1,n_h)})^{\delta} \sum_{j \in h, j \ne i} C_{ij} \cdot \kappa_{r_j} \cdot m_j \cdot \left(\beta_1 + \beta_2 \cdot f(VL_j(t))\right) \right)$$

| Parameter | Description |
|---|---|
| $\phi_{r_i}$ | Role-specific susceptibility multiplier |
| $\kappa_{r_j}$ | Role-specific infectivity multiplier |
| $\alpha_{\text{comm}}$ | Baseline Community risk |
| $S(t)$ | Seasonal forcing term |
| $\delta$ | Household size scaling |
| $C_{ij}$ | Contact weight between persons $i$ and $j$ |
| $f(VL)$ | Viral load scaling function (power-law or sigmoid) |
| $m_i$ | Intervention-driven modifier (e.g., $1 - \text{efficacy}$) |


Individuals progress through **S -> I -> R -> S** with Gamma-distributed durations. Inference is performed via HMC in Stan (`adapt_delta = 0.95`, `max_treedepth = 15`).

---

## License

MIT

## Citation

If you use HHBayes in your work, please cite:

> Ke Li & Yiren Hou (2026). HHBayes: Bayesian Household Transmission Modeling in R. GitHub: https://github.com/yourusername/HHBayes
