# HouseholdSim

**HouseholdSim** is an R package for simulating and estimating infectious disease transmission within households. It combines a stochastic household transmission simulator with Bayesian inference via Stan, supporting age-structured susceptibility/infectivity, viral load dynamics, seasonal community forcing, and reinfection.

## Features

- **Household simulation**: Simulate outbreaks across multiple households with configurable composition (adults, infants, toddlers, elderly), role-specific susceptibility and infectivity multipliers, and community-to-household importation.
- **Viral load modeling**: Empirical (double-exponential / piecewise-linear Ct) or mechanistic (within-host ODE) viral trajectories, with imputation of missing viral data.
- **Bayesian model fitting**: Prepare data and fit a pre-compiled Stan model estimating transmission rates, role-specific effects, covariate coefficients, and viral load scaling parameters.
- **Transmission chain reconstruction**: Recover who-infected-whom probabilities from posterior samples.
- **Visualization**: Posterior distributions, covariate forest plots, epidemic curves, and per-household infection timelines with transmission arrows.

## Installation

```r
# install.packages("devtools")
devtools::install_github("yourusername/HouseholdSim")
```

### Prerequisites

HouseholdSim depends on [RStan](https://mc-stan.org/rstan/), which requires a working C++ toolchain. See the [RStan Getting Started Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for setup instructions.

## Quick Start

### 1. Simulate households

```r
library(HouseholdSim)

sim <- simulate_multiple_households_comm(
  n_households = 100,
  beta1 = 0.008,
  beta2 = 0.008,
  phi_by_role = c(adult = 1, infant = 4, toddler = 5, elderly = 1),
  kappa_by_role = c(adult = 1, infant = 1, toddler = 1.2, elderly = 1),
  start_date = "2024-07-01",
  end_date = "2025-06-30",
  seed = 42
)
```

### 2. Summarize attack rates

```r
rates <- summarize_attack_rates(sim)
rates$primary_by_role
rates$reinf_by_role
```

### 3. Prepare data and fit the Stan model

```r
stan_data <- prepare_stan_data(
  df_clean = your_data,
  study_start_date = as.Date("2024-07-01"),
  study_end_date = as.Date("2025-06-30")
)

fit <- fit_household_model(stan_data, iter = 2000, chains = 4)
```

### 4. Visualize results

```r
# Posterior distributions of role-specific parameters
plot_posterior_distributions(fit)

# Covariate effects (forest plot)
plot_covariate_effects(fit, stan_data)

# Epidemic curve
plot_epidemic_curve(sim, surveillance_df, start_date_str = "2024-07-01")

# Household-level transmission timeline
trans_df <- reconstruct_transmission_chains(fit, stan_data)
plot_household_timeline(trans_df, stan_data, target_hh_id = 1)
```

## Key Functions

| Function | Description |
|---|---|
| `simulate_multiple_households_comm()` | Simulate multi-household outbreaks with community forcing |
| `prepare_stan_data()` | Format data for Stan model input |
| `fit_household_model()` | Fit the Bayesian household transmission model |
| `summarize_attack_rates()` | Compute primary attack rates and reinfection summaries |
| `reconstruct_transmission_chains()` | Extract who-infected-whom probabilities |
| `fill_missing_viral_data()` | Impute missing Ct / viral load values from theoretical curves |
| `plot_posterior_distributions()` | Violin plots of posterior estimates |
| `plot_covariate_effects()` | Forest plot of covariate coefficients |
| `plot_epidemic_curve()` | Epidemic curve stratified by age group |
| `plot_household_timeline()` | Timeline of infections within a single household |

## Model Overview

The transmission model estimates the daily probability of infection for each susceptible individual from two sources:

- **Community force of infection**: driven by external surveillance data or seasonal forcing curves.
- **Within-household force of infection**: a function of the number and infectiousness of currently infected household members, scaled by role-specific susceptibility (φ) and infectivity (κ) multipliers.

Viral load can optionally modulate infectiousness via an empirical trajectory or a mechanistic within-host ODE model. Immunity waning after recovery allows for reinfection.

Inference is performed using Hamiltonian Monte Carlo (HMC) in Stan.

## License

MIT

## Citation

If you use HouseholdSim in your work, please cite:

> [Your Name] (2025). HouseholdSim: Household Transmission Simulation and Inference in R. GitHub: https://github.com/yourusername/HouseholdSim

