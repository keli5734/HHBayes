# HHBayes <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/keli5734/HHBayes)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CRAN DOI](https://img.shields.io/badge/DOI-10.32614%2FCRAN.package.HHBayes-blue.svg)](https://doi.org/10.32614/CRAN.package.HHBayes)
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

## Authors

- **Ke Li** — Package creator and maintainer
- **Yiren Hou** — Package creator and maintainer

## Citation

If you use HHBayes in your work, please cite:
> Li K, Hou Y, Mukherjee B, Pitzer VE, Weinberger DM (2026). HHBayes: A Flexible Bayesian Framework for Simulating and Analyzing Household Transmission Dynamics. medRxiv. https://www.medrxiv.org/content/10.64898/2026.04.01.26349903v1

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("keli5734/HHBayes")
```

### Prerequisites

HHBayes depends on [RStan](https://mc-stan.org/rstan/), which requires a working C++ toolchain. See the [RStan Getting Started Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for platform-specific setup.

---

## A Note on Prior Scales

Several parameters are **sampled on the log scale** and then exponentiated to the
biological (natural) scale inside the Stan model. This matters when you specify
priors:

- **Log-scale parameters:** `phi_role`, `kappa_role` (role susceptibility /
  infectivity multipliers) and `beta1`, `beta2`, `alpha` (transmission and
  community-acquisition rates). For these, a `dist = "normal"` prior on the log
  scale **is** a LogNormal prior on the biological quantity, and a
  `dist = "uniform"` prior is a log-uniform prior on the biological quantity.
  A `"lognormal"` prior is **not allowed** for these parameters: on a log-scale
  parameter it becomes one-sided (it forces the biological multiplier to be > 1),
  which is not a lognormal on the biological scale. `prepare_stan_data()` will
  stop with an informative error if you request it.
- **Natural-scale parameters:** `vl_midpoint`, `vl_slope`, `gen_shape`,
  `gen_rate`. These are sampled directly on the biological scale and accept
  `"normal"`, `"uniform"`, or `"lognormal"` priors.
- **Covariate effects** (`beta_susc`, `beta_inf`) are log-linear coefficients and
  accept `"normal"` or `"uniform"` priors only.

> **Rule of thumb:** to place a LogNormal(μ, σ) prior on a role multiplier or a
> transmission rate, use `dist = "normal", params = c(μ, σ)`.

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

# 3. Plot simulated data
my_plot <- plot_epidemic_curve(sim, surveillance_data, start_date_str = study_start, bin_width = 7)
print(my_plot)

# 4. Prepare data and fit the Bayesian model
df_for_stan <- sim$diagnostic_df

# NOTE: phi_role, kappa_role, beta1, beta2, alpha are sampled on the LOG scale,
# so dist = "normal" here == LogNormal on the biological quantity. See
# "A Note on Prior Scales" above.
my_priors <- list(
    beta1      = list(dist = "normal",  params = c(-5, 1)),
    beta2      = list(dist = "normal",  params = c(-5, 1)),
    alpha      = list(dist = "normal",  params = c(-7, 1)),
    phi_role   = list(dist = "normal",  params = c(0, 0.5)),
    kappa_role = list(dist = "normal",  params = c(0, 0.5)),
    vl_midpoint = list(dist = "normal", params = c(3, 0.5)),  # adjust if using Ct values
    vl_slope    = list(dist = "normal", params = c(2.5,  0.5))   # adjust if using Ct values
)

VL_params_list <- list(
    adult   = list(v_p=4.14, t_p=5.09, lambda_g=2.31, lambda_d=2.71),
    infant  = list(v_p=5.84, t_p=4.09, lambda_g=2.82, lambda_d=1.01),
    toddler = list(v_p=5.84, t_p=4.09, lambda_g=2.82, lambda_d=1.01),
    elderly = list(v_p=2.95, t_p=5.1,  lambda_g=3.15, lambda_d=0.87)
)

stan_input <- prepare_stan_data(
    df_clean          = df_for_stan,
    surveillance_df   = surveillance_data,
    study_start_date  = as.Date(study_start),
    study_end_date    = as.Date(study_end),
    use_vl_data       = TRUE,
    use_curve_logic   = FALSE,
    delta             = 0,
    imputation_params = VL_params_list,
    priors            = my_priors
)

options(mc.cores = parallel::detectCores())

# fit_household_model() automatically drops internal bookkeeping quantities AND
# any viral-component parameters that are inert under this configuration, so the
# printed fit shows only parameters informed by the likelihood.
fit <- fit_household_model(stan_input, iter = 2000, warmup = 1000, chains = 4, seed = 123, cores = 4)

# 5. Visualize
print(fit, probs = c(0.025, 0.5, 0.975))
p_post <- plot_posterior_distributions(fit)
chains <- reconstruct_transmission_chains(fit = fit, stan_data = stan_input, min_prob_threshold =  .01)
selected_hh = 1
p_hh <- plot_household_timeline(chains, stan_input, target_hh_id = selected_hh) # define selected_hh.
print(p_hh)
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

By default, all household members contact each other equally. A `role_mixing_matrix` lets you specify **differential contact weights** between roles.

**Orientation convention (important for asymmetric matrices):** each cell `[target, source]` is the contact weight *from* the source role (column) *to* the target role (row). This matches how the weight enters the force of infection — force on a susceptible target accumulates the infectivity of each source scaled by this weight. For symmetric matrices the orientation is irrelevant; for asymmetric ones, make sure rows are targets and columns are sources.

```r
role_mixing_weights <- matrix(c(
# Source:  Infant Toddler Adult Elderly
          0.0,   0.5,    1.0,  0.5,    # Target: Infant  (receives high weight from adults)
          0.5,   0.9,    0.7,  0.5,    # Target: Toddler (high peer contact)
          1.0,   0.7,    0.6,  0.7,    # Target: Adult
          0.5,   0.5,    0.7,  0.0     # Target: Elderly (limited infant contact)
), nrow = 4, byrow = TRUE,
   dimnames = list(
     c("infant", "toddler", "adult", "elderly"),   # rows = targets
     c("infant", "toddler", "adult", "elderly")))  # cols = sources
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

Interventions (e.g., vaccination) are specified as a list via `covariates_config`. Each entry defines an intervention with a name, what it affects, its efficacy, and role-specific coverage probabilities:

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

**Define priors** for the Bayesian model.

Supported distributions depend on the parameter's scale (see "A Note on Prior Scales" above):

- `phi_role`, `kappa_role`, `beta1`, `beta2`, `alpha` — sampled on the **log
  scale**; use `"normal"` (== LogNormal on the biological scale) or `"uniform"`
  (== log-uniform). `"lognormal"` is rejected for these.
- `vl_midpoint`, `vl_slope`, `gen_shape`, `gen_rate` — sampled on the **natural
  scale**; accept `"normal"`, `"uniform"`, or `"lognormal"`.
- covariate effects — accept `"normal"` or `"uniform"`.

```r
my_priors <- list(
  beta1      = list(dist = "normal",  params = c(-5, 1)),
  beta2      = list(dist = "normal",  params = c(-5, 1)),
  alpha      = list(dist = "normal",  params = c(-7, 1)),
  phi_role   = list(dist = "normal",  params = c(0, 1)),   # LogNormal(0,1) on the multiplier
  kappa_role = list(dist = "normal",  params = c(0, 1)),   # LogNormal(0,1) on the multiplier
  vl_midpoint = list(dist = "normal", params = c(4, 1.0)),  # adjust if Ct values
  vl_slope    = list(dist = "normal", params = c(1.0,  0.5))   # adjust if Ct values
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
  use_vl_data      = TRUE,
  use_curve_logic  = FALSE,
  delta = 0,
  imputation_params = VL_params_list,
  priors           = my_priors,
  role_mixing_matrix = role_mixing_weights)
```

### Step 9: Fit the Bayesian Model

`fit_household_model()` runs HMC sampling via `rstan::sampling()`. By default
(`drop_inert = TRUE`) it automatically excludes internal bookkeeping quantities
**and** any viral-component parameters that are inert under your data
configuration, so the printed fit shows only the parameters that were actually
informed by the likelihood.

```r
options(mc.cores = parallel::detectCores())

fit <- fit_household_model(
  stan_input,
  iter   = 2000,
  warmup = 1000,
  chains = 4,
  seed   = 123,
  cores  = 4
)

print(fit, probs = c(0.025, 0.5, 0.975))
```

**Which parameters are dropped?** `prepare_stan_data()` records the inert set in
`stan_input$inert_pars` based on your `use_vl_data` / `use_curve_logic` choice,
and `fit_household_model()` excludes them from the returned object:

| Configuration | Viral component used | Inert (dropped from output) |
|---|---|---|
| `use_vl_data = TRUE` | Viral-load / Ct scaling (`vl_midpoint`, `vl_slope`) | `gen_shape`, `gen_rate` |
| `use_vl_data = FALSE`, `use_curve_logic = TRUE` | Gamma generation-interval curve (`gen_shape`, `gen_rate`) | `vl_midpoint`, `vl_slope` |
| `use_vl_data = FALSE`, `use_curve_logic = FALSE` | Constant infectiousness | `gen_shape`, `gen_rate`, `vl_midpoint`, `vl_slope` |

Inert parameters are still *sampled* internally (from a neutral pinning prior)
but never mistaken for estimates, because they are removed from `print(fit)` and
`rstan::extract(fit)`.

> **Do not pass your own `pars`/`include`** unless you want to take over
> parameter selection entirely — doing so disables the automatic dropping and the
> inert parameters will reappear in the output. To keep *everything* for
> debugging, use `fit_household_model(stan_input, drop_inert = FALSE)`.

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
| `vl_midpoint` | `3.0 or 40` | Reference log10 VL for infectiousness scaling or Ct at 50% infectiousness (sigmoid) |
| `vl_slope` | `2.5 or 2` | Power-law exponent for log10 VL scaling. or Steepness of Ct-infectiousness sigmoid. |

### Testing and Surveillance

| Parameter | Default | Description |
|---|---|---|
| `surveillance_interval` | `1` | Days between routine tests. |
| `test_daily` | `FALSE` | Switch to daily testing after first household detection. |
| `perfect_detection` | `TRUE` | If `FALSE`, detection depends on viral load thresholds. |
| `detect_threshold_log10` | `1e-6` | Min log10 VL for positive test. |
| `detect_threshold_Ct` | `40` | Max Ct for positive test (set `35`-`40` for realistic PCR). |

### Data Preparation & Model Control

| Parameter | Default | Description |
|---|---|---|
| `use_vl_data` | `TRUE` | Use actual viral load data vs. generation interval curve. |
| `use_curve_logic` | `TRUE` | When `use_vl_data = FALSE`, estimate generation interval curve vs. constant infectivity. |
| `delta` | `0` | Household size scaling. Force scaled by `(1/max(n_h,1))^delta`. |

### Viral Load Parameters (when `vl_type = 1`)

| Parameter | Default | Description |
|---|---|---|
| `vl_midpoint` | `3.0` | Reference log10 VL for infectiousness scaling. |
| `vl_slope` | `2.5` | Power-law exponent for VL scaling: `(VL/vl_midpoint)^vl_slope`. |

### Ct Parameters (when `vl_type = 0`)

| Parameter | Default | Description |
|---|---|---|
| `vl_midpoint` | `35.0` | Ct value at 50% infectiousness. |
| `vl_slope` | `2.0` | Steepness of Ct-infectiousness sigmoid. |

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

$$P(\text{infection}_{ih}(t)) = 1 - \exp\left(-\lambda_{ih}(t)\right)$$

where the total force of infection combines community and household sources:

$$\lambda_{ih}(t) = \phi_{r_i} \cdot \exp(\boldsymbol{\beta}_{\text{susc}} \cdot \mathbf{X}_{\text{susc},i}) \cdot \left( \alpha_{\text{comm}} \cdot S(t) + (\frac{1}{\max(1,n_h)})^{\delta} \cdot \sum_{j \in h, j \ne i} C_{ij} \cdot \kappa_{r_j} \cdot \exp(\boldsymbol{\beta}_{\text{inf}} \cdot \mathbf{X}_{\text{inf},j}) \cdot \left(\beta_1 + \beta_2 \cdot f(VL_j(t))\right) \right)$$

| Parameter | Description |
|---|---|
| $\phi_{r_i}$ | Role-specific susceptibility multiplier |
| $\kappa_{r_j}$ | Role-specific infectivity multiplier |
| $\alpha_{\text{comm}}$ | Baseline Community risk |
| $S(t)$ | Seasonal forcing term |
| $\delta$ | Household size scaling |
| $C_{ij}$ | Contact weight between persons $i$ and $j$ |
| $f(VL)$ | Viral load scaling function (power-law or sigmoid) |
| $\boldsymbol{\beta}_{\text{susc}}$ | Log-linear coefficients for susceptibility covariates |
| $\boldsymbol{\beta}_{\text{inf}}$ | Log-linear coefficients for infectivity covariates |
| $\mathbf{X}_{\text{susc},i}$ | Covariate vector for person $i$ (susceptibility) |
| $\mathbf{X}_{\text{inf},j}$ | Covariate vector for person $j$ (infectivity) |

Individuals progress through **S -> I -> R -> S** with Gamma-distributed duration. Inference is performed via HMC in Stan (`adapt_delta = 0.95`, `max_treedepth = 15`).

### Parameterization and Prior Scale

The role multipliers and rates are parameterized on the **log scale** and mapped to the biological scale inside the model:

$$\phi_{r} = \phi_{\text{ref}} \cdot \exp(\theta^{\phi}_{r}), \qquad
\kappa_{r} = \kappa_{\text{ref}} \cdot \exp(\theta^{\kappa}_{r}), \qquad
\beta_1 = \exp(\theta^{\beta_1}), \ \text{etc.}$$

Because the sampled quantity is $\theta^{\phi}_{r} = \log \phi_r$, a $\text{Normal}(\mu,\sigma)$ prior placed on it corresponds to a $\text{LogNormal}(\mu,\sigma)$ prior on the biological multiplier $\phi_r$, and a $\text{Uniform}$ prior corresponds to a log-uniform prior. Priors on the viral-load scaling parameters (`vl_midpoint`, `vl_slope`) and generation-interval parameters (`gen_shape`, `gen_rate`) are placed directly on the natural (biological) scale. This is why the software disallows a `"lognormal"` option for the log-scale parameters: applied to a log-scale parameter it would be one-sided and would not represent a lognormal prior on the biological quantity.

---

#### **How Covariate Effects Work**

The simulation applies **fixed efficacies** during data generation, but the Stan model estimates **log-linear coefficients** to recover those effects:

**Simulation time:**

```r
# Vaccination reduces susceptibility by 50%
if (vaccinated) susceptibility *= 0.5
```

**Stan estimation:**

```r
# Stan estimates beta_susc ≈ -0.693, because exp(-0.693) ≈ 0.5
susceptibility = phi_role * exp(beta_susc * vacc_status)
```

**Interpretation:**

- `β = 0`: No effect (`exp(0) = 1`)
- `β = -0.693`: 50% reduction (`exp(-0.693) ≈ 0.5`)
- `β = +0.693`: 2× increase (`exp(0.693) ≈ 2.0`)
- Multiple covariates: Effects multiply on natural scale, add on log scale

---

### Viral Component Scenarios

The model supports **four scenarios** for modeling infectiousness, controlled by `use_vl_data`, `vl_type`, and `use_curve_logic`:

**Scenario 1: Log10 Viral Load Data** (`use_vl_data = TRUE, vl_type = 1`)

```r
v_component = (max(0, VL) / vl_midpoint)^vl_slope
```

**Estimated parameters:** `vl_midpoint`, `vl_slope`

**Scenario 2: Ct Value Data** (`use_vl_data = TRUE, vl_type = 0`)

```r
v_component = 1 / (1 + exp((VL - vl_midpoint) / vl_slope))
```

**Estimated parameters:** `vl_midpoint`, `vl_slope`

**Scenario 3: Generation Interval Curve** (`use_vl_data = FALSE, use_curve_logic = TRUE`)

```r
v_component = g_curve[days_since_infection]
```

Where `g_curve` is a normalized Gamma distribution.

**Estimated parameters:** `gen_shape`, `gen_rate`

**Scenario 4: Constant Infectiousness** (`use_vl_data = FALSE, use_curve_logic = FALSE`)

```r
v_component = 1.0  # No viral dynamics
```

**No viral-component parameters estimated.**

The final infectiousness is always: `beta1 + beta2 * v_component`

> In every scenario, the viral-component parameters that do **not** apply are
> inert (sampled from a neutral prior, not informed by the data) and are dropped
> from the fitted object by `fit_household_model()`. See Step 9.

---

## Estimated Parameters

The Stan model uses **conditional parameter estimation** — which viral-component parameters are informed by the likelihood depends on your data configuration. Parameters that are not informed under a given configuration are inert (prior-only) and are dropped from the fitted object.

### Always Estimated

| Parameter | Description |
|---|---|
| `beta1`, `beta2` | Baseline and viral-dependent transmission rates (log-scale priors) |
| `alpha_comm` | Community importation rate (log-scale prior) |
| `phi_role` | Role-specific susceptibility multipliers (log-scale priors) |
| `kappa_role` | Role-specific infectivity multipliers (log-scale priors) |

### Conditionally Estimated

| Parameter | Description | When Estimated |
|---|---|---|
| `vl_midpoint`, `vl_slope` | Viral-load / Ct infectiousness scaling | `use_vl_data = TRUE` |
| `gen_shape`, `gen_rate` | Generation-interval curve | `use_vl_data = FALSE` AND `use_curve_logic = TRUE` |
| `beta_susc` | Susceptibility covariate effects | `covariates_susceptibility` provided |
| `beta_inf` | Infectivity covariate effects | `covariates_infectivity` provided |

> Under Scenario 4 (`use_vl_data = FALSE`, `use_curve_logic = FALSE`) none of the
> viral-component parameters (`vl_midpoint`, `vl_slope`, `gen_shape`, `gen_rate`)
> are estimated.

### Flexible Prior Support

Priors are placed on the scale at which each parameter is sampled (see "A Note on Prior Scales"):

- **Log-scale** (`phi_role`, `kappa_role`, `beta1`, `beta2`, `alpha`): `"normal"`
  or `"uniform"`. A `"normal"` prior here is a LogNormal prior on the biological
  quantity. `"lognormal"` is **not permitted** and raises an error.
- **Natural-scale** (`vl_midpoint`, `vl_slope`, `gen_shape`, `gen_rate`):
  `"normal"`, `"uniform"`, or `"lognormal"`.
- **Covariate effects**: `"normal"` or `"uniform"`.

```r
my_priors <- list(
  phi_role   = list(dist = "normal", params = c(0, 0.5)),  # LogNormal(0,0.5) on the multiplier
  kappa_role = list(dist = "normal", params = c(0, 0.3)),  # LogNormal(0,0.3) on the multiplier
  vl_slope   = list(dist = "lognormal", params = c(0, 0.3)) # natural-scale param: lognormal OK
)
```

### Example

```r
# Baseline: no covariates → estimates core params
stan_input <- prepare_stan_data(
  ...,
  use_vl_data = TRUE,
  covariates_susceptibility = NULL,
  covariates_infectivity    = NULL
)

# With covariates → also estimates vaccine effects
stan_input <- prepare_stan_data(
  ...,
  use_vl_data = TRUE,
  covariates_susceptibility = c("vacc_status"),
  covariates_infectivity    = c("vacc_status")
)
```

When `K_susc = 0` or `K_inf = 0`, the corresponding coefficient vectors do not enter the model. This means the same Stan model adapts automatically — no need to switch between different model files.

---

## Acknowledgments

This work was supported by a grant from the National Institutes of Health (R01AI137093). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
