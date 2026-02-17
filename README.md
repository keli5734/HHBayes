# HHBayes <img src="man/figures/logo.png" align="right" height="139" />

> **Bayesian Household Transmission Modeling in R**

**HHBayes** is an R package for simulating and estimating infectious disease transmission within households. It couples a stochastic, age-structured household simulator with Bayesian inference via [Stan](https://mc-stan.org/), enabling researchers to study within-household spread, evaluate intervention strategies, and reconstruct transmission chains from longitudinal diagnostic data.

---

## Why HHBayes?

Household transmission studies generate rich but complex data — repeated tests on every member, overlapping infections, reinfections, and age-dependent viral kinetics. HHBayes provides an integrated pipeline:

1. **Simulate** realistic outbreaks with configurable household structures, contact patterns, and interventions.
2. **Prepare** observed or simulated data for Bayesian analysis.
3. **Fit** a Stan model that jointly estimates transmission rates, role-specific susceptibility/infectivity, covariate effects, and viral load scaling.
4. **Analyze** results with built-in tools for attack rates, transmission chain reconstruction, and publication-ready plots.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("keli5734/HHBayes")
```

### Prerequisites

HHBayes depends on [RStan](https://mc-stan.org/rstan/), which requires a working C++ toolchain. See the [RStan Getting Started Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for platform-specific setup instructions.

---

## Quick Start

```r
library(HHBayes)

# Simulate 200 households over one year
sim <- simulate_multiple_households_comm(
  n_households = 200,
  start_date   = "2024-07-01",
  end_date     = "2025-06-30",
  beta1 = 0.008, beta2 = 0.008,
  seed  = 42
)

# Summarize attack rates by age group
rates <- summarize_attack_rates(sim)
rates$primary_by_role

# Prepare data and fit the Bayesian model
stan_data <- prepare_stan_data(sim$diagnostic_df)
fit <- fit_household_model(stan_data, iter = 2000, chains = 4)

# Visualize posteriors
plot_posterior_distributions(fit)
```

---

## Simulation Guide

The core simulation function is `simulate_multiple_households_comm()`. Below is a detailed guide to its inputs, organized by topic.

### Basic Setup

| Parameter | Default | Description |
|---|---|---|
| `n_households` | `50` | Number of households to simulate. |
| `start_date` | `"2024-07-01"` | Simulation start date (YYYY-MM-DD). |
| `end_date` | `"2025-06-30"` | Simulation end date (YYYY-MM-DD). |
| `seed` | `NULL` | Random seed for reproducibility. Set an integer for reproducible runs. |
| `verbose` | `FALSE` | If `TRUE`, prints progress for each household. |

---

### Household Structure

Households are randomly generated from a demographic profile. Each household contains a mix of four **roles**: `adult`, `infant`, `toddler`, and `elderly`. You control the composition via `household_profile_list`:

```r
# Default profile
household_profile_list = list(
  prob_single_parent = 0,                  # Probability of single-parent household
  prob_siblings      = c(0.10, 0.50, 0.40),# P(0, 1, or 2 toddler siblings)
  prob_elderly       = c(0.9, 0.08, 0.02)  # P(0, 1, or 2 elderly members)
)

# Example: Multi-generational households common in East/South Asia
asian_profile <- list(
  prob_single_parent = 0.05,
  prob_siblings      = c(0.05, 0.30, 0.65),  # More siblings
  prob_elderly       = c(0.20, 0.50, 0.30)   # Grandparents often present
)

sim <- simulate_multiple_households_comm(
  n_households = 200,
  household_profile_list = asian_profile
)
```

---

### Transmission Parameters

| Parameter | Default | Description |
|---|---|---|
| `beta1` | `0.008` | Baseline within-household transmission rate (contact-driven). |
| `beta2` | `0.008` | Viral-load-dependent transmission rate. Total per-infector force = `beta1 + beta2 * f(viral_load)`. |
| `delta` | `0` | Household size scaling exponent. If `delta > 0`, force is scaled by `(1/N)^delta`, reducing per-contact transmission in larger households. |
| `alpha_comm_by_role` | `5e-4` | Daily community acquisition rate. A single value applies equally to all roles. |
| `max_infections` | `Inf` | Maximum number of times a single person can be infected. Set to `1` to disable reinfection. |

---

### Role-Specific Susceptibility and Infectivity

Each role has a **susceptibility multiplier** (phi) and an **infectivity multiplier** (kappa):

```r
sim <- simulate_multiple_households_comm(
  phi_by_role   = c(adult = 1, infant = 4, toddler = 5, elderly = 1),
  kappa_by_role = c(adult = 1, infant = 1, toddler = 1.2, elderly = 1)
)
```

- **`phi_by_role`**: How susceptible each role is relative to the baseline. An infant with `phi = 4` is 4x more likely to become infected per unit of force.
- **`kappa_by_role`**: How infectious each role is. A toddler with `kappa = 1.2` contributes 20% more force to susceptible contacts.

---

### Infectious Period and Immunity Waning

Both durations are drawn from Gamma distributions:

| Parameter | Default | Description |
|---|---|---|
| `infectious_shape` | `3` | Shape of the Gamma distribution for the infectious period. |
| `infectious_scale` | `1` | Scale of the Gamma distribution for the infectious period. Mean duration = shape x scale = 3 days. |
| `waning_shape` | `16` | Shape of the Gamma distribution for immunity duration after recovery. |
| `waning_scale` | `10` | Scale of the Gamma distribution for immunity duration. Mean = 160 days with defaults. |
| `peak_day` | `1` | Day of peak infectiousness (relative to infection onset). |
| `width` | `4` | Width of the infectiousness peak window. |

---

### Contact Matrix

By default, all household members contact each other equally (fully connected). You can specify **differential contact patterns** between roles using `role_mixing_matrix`:

```r
# Define a 4x4 contact weight matrix between roles
# Values are relative contact weights (1.0 = baseline)
mixing <- matrix(c(
# infant toddler adult elderly
  0.5,   1.0,    1.5,  0.3,   # infant's contacts with others
  1.0,   1.0,    0.8,  0.2,   # toddler's contacts
  1.5,   0.8,    1.0,  0.7,   # adult's contacts
  0.3,   0.2,    0.7,  0.5    # elderly's contacts
), nrow = 4, byrow = TRUE,
   dimnames = list(
     c("infant", "toddler", "adult", "elderly"),
     c("infant", "toddler", "adult", "elderly")
   ))

sim <- simulate_multiple_households_comm(
  role_mixing_matrix = mixing
)
```

This matrix is expanded to an N x N individual-level contact matrix for each household based on the roles of its members. For example, in a household with roles `c("adult", "adult", "infant")`, the `adult-infant` weight (1.5) is applied to the corresponding pairs.

Alternatively, you can pass `contact_mat` — a pre-built N x N matrix — for a fixed household size, but `role_mixing_matrix` is recommended as it generalizes across varying household sizes.

---

### Community Forcing and Seasonality

Community importation can be driven by external surveillance data or a custom seasonal curve.

**Option A: Real surveillance data** — provide a dataframe with `date` and `cases` columns. The case counts are automatically normalized to [0, 1] and interpolated to daily resolution:

```r
surv_df <- data.frame(
  date  = seq(as.Date("2024-07-01"), as.Date("2025-06-30"), by = "week"),
  cases = c(10, 15, 30, 80, 150, ...)  # Weekly case counts
)

sim <- simulate_multiple_households_comm(
  surveillance_df = surv_df
)
```

**Option B: Custom forcing curves** — provide a named list of daily forcing vectors (one per role, length = number of simulation days):

```r
n_days <- 366
# Sinusoidal winter peak
winter_curve <- 0.5 + 0.5 * sin(2 * pi * (1:n_days - 90) / 365)

sim <- simulate_multiple_households_comm(
  seasonal_forcing_list = list(
    adult   = winter_curve,
    infant  = winter_curve * 1.2,  # Infants have higher community exposure
    toddler = winter_curve,
    elderly = winter_curve * 0.8   # Elderly have less community contact
  )
)
```

If neither is provided, community forcing is constant (flat = 1.0) across the entire simulation period.

---

### Intervention Strategies (Covariates)

Interventions (e.g., vaccination, masking, prophylaxis) are specified via `covariates_config` — a list where each element defines one intervention. Each intervention has a name, a target effect, an efficacy, and role-specific coverage probabilities:

```r
interventions <- list(
  # Intervention 1: Vaccination
  list(
    name      = "vaccinated",
    effect_on = "susceptibility",   # "susceptibility", "infectivity", or "both"
    efficacy  = 0.6,                # 60% reduction in susceptibility
    coverage  = list(
      adult   = 0.80,   # 80% of adults vaccinated
      elderly = 0.90,   # 90% of elderly
      toddler = 0.30,   # 30% of toddlers
      infant  = 0.00    # Infants not eligible
    )
  ),
  # Intervention 2: Masking
  list(
    name      = "masked",
    effect_on = "both",             # Reduces both susceptibility and infectivity
    efficacy  = 0.3,                # 30% reduction
    coverage  = list(
      adult   = 0.70,
      elderly = 0.60,
      toddler = 0.10,
      infant  = 0.00
    )
  )
)

sim <- simulate_multiple_households_comm(
  covariates_config = interventions
)
```

**How efficacy is applied**: For each person, a Bernoulli draw determines whether they receive the intervention (based on `coverage`). If they do, their susceptibility and/or infectivity is multiplied by `(1 - efficacy)`. Multiple interventions stack multiplicatively. For example, a vaccinated (`efficacy = 0.6`) and masked (`efficacy = 0.3`) adult has their susceptibility multiplied by `(1 - 0.6) * (1 - 0.3) = 0.28`.

---

### Viral Load and Ct Value Dynamics

HHBayes supports two types of viral testing and two modeling approaches for within-host viral trajectories.

#### Testing type (`viral_testing`)

| Value | Description |
|---|---|
| `"viral load"` | Log10 viral load scale (higher = more virus). Positive test when value `>= detect_threshold_log10`. |
| `"Ct"` | Cycle threshold scale (lower = more virus). Positive test when value `<= detect_threshold_Ct`. |

#### Modeling approach (`model_type`)

| Value | Description |
|---|---|
| `"empirical"` | Parametric curves: double-exponential for viral load, piecewise-linear for Ct. Fast and flexible. |
| `"ODE"` | Within-host ODE system (target cell-limited model with innate immunity via deSolve). More mechanistic but slower. |

#### Viral load trajectory parameters

Role-specific parameters define the shape of each person's viral trajectory:

```r
# Empirical viral load (double-exponential) parameters
VL_params_list = list(
  adult   = list(v_p = 4.14, t_p = 5.09, lambda_g = 2.31, lambda_d = 2.71),
  infant  = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  toddler = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  elderly = list(v_p = 2.95, t_p = 5.10, lambda_g = 3.15, lambda_d = 0.87)
)
# v_p:       peak log10 viral load
# t_p:       time to peak (days post-infection)
# lambda_g:  growth rate
# lambda_d:  decay rate

# Empirical Ct trajectory (piecewise-linear) parameters
Ct_params_list = list(
  adult   = list(Cpeak = 33, r = 1.49, d = 1.22, t_peak = 5.14),
  infant  = list(Cpeak = 33.3, r = 2.11, d = 1.38, t_peak = 5.06),
  toddler = list(Cpeak = 34, r = 1.26, d = 1.27, t_peak = 4.75),
  elderly = list(Cpeak = 33, r = 1.49, d = 1.22, t_peak = 5.14)
)
# Cpeak:  minimum Ct at peak (lower Ct = more virus)
# r:      Ct rise rate before peak (Ct/day)
# d:      Ct decline rate after peak (Ct/day)
# t_peak: time to peak (days post-infection)
```

#### Transmission scaling from viral load

| Parameter | Default | Description |
|---|---|---|
| `V_ref` | `3.0` | Reference log10 viral load. Infectiousness scales as `(VL / V_ref)^V_rho`. |
| `V_rho` | `2.5` | Power-law exponent for viral load scaling. Higher values = steeper dose-response. |
| `Ct_50` | `40` | Ct value at 50% of maximum infectiousness (sigmoid function). |
| `Ct_delta` | `2` | Steepness of the Ct-infectiousness sigmoid. Smaller = sharper transition. |

---

### Testing and Surveillance Strategy

These parameters control the diagnostic testing protocol:

| Parameter | Default | Description |
|---|---|---|
| `surveillance_interval` | `1` | Days between routine surveillance tests. Set to `7` for weekly, `14` for biweekly. |
| `test_daily` | `FALSE` | If `TRUE`, switches to daily testing after the **first positive** in a household. Mimics reactive enhanced surveillance. |
| `perfect_detection` | `TRUE` | If `TRUE`, any shedding individual is always detected. Set to `FALSE` for realistic imperfect testing. |
| `detect_threshold_log10` | `1e-6` | Minimum log10 viral load for a positive result (viral load mode). |
| `detect_threshold_Ct` | `99` | Maximum Ct for a positive result (Ct mode). Set to `35`-`40` for realistic PCR thresholds. |

**Example: Weekly screening with reactive daily follow-up**

```r
sim <- simulate_multiple_households_comm(
  viral_testing         = "Ct",
  surveillance_interval = 7,          # Routine weekly swabs
  test_daily            = TRUE,       # Daily testing after first detection
  perfect_detection     = FALSE,      # Imperfect PCR detection
  detect_threshold_Ct   = 35          # Only detect if Ct <= 35
)
```

**Example: Comparing testing frequencies**

```r
# Weekly testing
sim_weekly <- simulate_multiple_households_comm(surveillance_interval = 7, seed = 1)

# Daily testing
sim_daily <- simulate_multiple_households_comm(surveillance_interval = 1, seed = 1)

# Compare detection rates
nrow(sim_weekly$diagnostic_df[sim_weekly$diagnostic_df$test_result == 1, ])
nrow(sim_daily$diagnostic_df[sim_daily$diagnostic_df$test_result == 1, ])
```

---

## Bayesian Model Fitting

### Prepare data for Stan

```r
stan_data <- prepare_stan_data(
  df_clean           = cleaned_observations,
  study_start_date   = as.Date("2024-07-01"),
  study_end_date     = as.Date("2025-06-30"),
  role_levels        = c("adult", "infant", "toddler", "elderly"),
  use_vl_data        = TRUE,
  role_mixing_matrix = mixing,
  covariates_susceptibility = c("vaccinated"),
  covariates_infectivity    = c("masked"),
  seed = 123
)
```

### Fit the model

```r
fit <- fit_household_model(
  stan_data,
  iter   = 2000,
  chains = 4,
  warmup = 1000
)
```

### Analyze results

```r
# Attack rates and reinfection summaries
rates <- summarize_attack_rates(sim)

# Reconstruct who-infected-whom
trans_df <- reconstruct_transmission_chains(fit, stan_data, min_prob_threshold = 0.01)

# Impute missing viral data
filled_df <- fill_missing_viral_data(df, "ct_value", type = "ct_value",
                                      params_list = Ct_params_list,
                                      detection_limit = 45)
```

---

## Visualization

```r
# Posterior distributions of phi (susceptibility) and kappa (infectivity)
plot_posterior_distributions(fit)

# Forest plot of covariate effects (vaccination, masking, etc.)
plot_covariate_effects(fit, stan_data)

# Epidemic curve stratified by age group
plot_epidemic_curve(sim, surveillance_df, start_date_str = "2024-07-01")

# Detailed household timeline with transmission arrows
plot_household_timeline(trans_df, stan_data, target_hh_id = 1)
```

---

## Key Functions

| Function | Description |
|---|---|
| `simulate_multiple_households_comm()` | Simulate multi-household outbreaks with community forcing |
| `prepare_stan_data()` | Format observed/simulated data for the Stan model |
| `fit_household_model()` | Fit the Bayesian household transmission model via HMC |
| `summarize_attack_rates()` | Compute primary attack rates and reinfection summaries by role |
| `reconstruct_transmission_chains()` | Extract posterior who-infected-whom probabilities |
| `fill_missing_viral_data()` | Impute missing Ct or viral load values from theoretical curves |
| `plot_posterior_distributions()` | Violin plots of role-specific posterior estimates |
| `plot_covariate_effects()` | Forest plot of intervention/covariate coefficients |
| `plot_epidemic_curve()` | Weekly/daily epidemic curve stratified by age group |
| `plot_household_timeline()` | Per-household infection timeline with transmission arrows |

---

## Model Overview

The transmission model estimates the daily probability of infection for each susceptible individual *i* in household *h* at time *t*:

```
P(infection_i(t)) = 1 - exp( -lambda_i(t) )
```

where the total force of infection combines community and household sources:

```
lambda_i(t) = phi[role_i] * modifier_i * (
    alpha_comm[role_i](t)
  + SUM_j { C[i,j] * kappa[role_j] * modifier_j * (beta1 + beta2 * f(VL_j(t))) }
)
```

| Symbol | Meaning |
|---|---|
| phi | Role-specific susceptibility multiplier |
| kappa | Role-specific infectivity multiplier |
| alpha_comm(t) | Time-varying community importation rate |
| C[i,j] | Contact weight between persons *i* and *j* |
| f(VL) | Viral load scaling function (power-law or sigmoid) |
| modifier | Intervention-driven reduction (e.g., vaccine efficacy) |

After infection, individuals progress through **S -> I -> R -> S** with Gamma-distributed durations for the infectious period and immunity waning, allowing reinfection.

Inference is performed using Hamiltonian Monte Carlo (HMC) via Stan with adaptive control parameters (`adapt_delta = 0.95`, `max_treedepth = 15`).

---

## License

MIT

## Citation

If you use HHBayes in your work, please cite:

> [Your Name] (2025). HHBayes: Bayesian Household Transmission Modeling in R. GitHub: https://github.com/yourusername/HHBayes
