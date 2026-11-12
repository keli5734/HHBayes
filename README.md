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
2. **Estimate viral kinetics and infection timing (Stage 1)** from the viral load or Ct measurements, including the negative tests before and after each infection episode.
3. **Prepare** observed or simulated data for Bayesian analysis, using the Stage 1 infection times, infectious windows, and viral load curves (or a simple random imputation when no quantitative viral data are available).
4. **Fit** a Stan model of household transmission (Stage 2) that estimates transmission rates, role-specific susceptibility/infectivity, covariate effects, and viral load scaling, and pool fits across Stage 1 draws so that the uncertainty in infection timing is carried into the estimates.
5. **Analyze** results with built-in tools for attack rates, transmission chain reconstruction, and publication-ready plots.

---

## What's New

This release adds a two-stage estimation framework and makes the simulator and the estimator consistent with each other.

- **Stage 1 viral kinetics model** (`prepare_kinetics_data()`, `fit_viral_kinetics()`): a hierarchical Stan model (`viral_kinetics.stan`) that estimates viral load (log10) or Ct trajectories by age group and each episode's time from infection to the first positive test. Negative tests before and after an episode enter as censored observations.
- **Uncertainty in infection timing:** `kinetics_draw_to_timelines()` turns one Stage 1 posterior draw into infection dates, infectious windows, and viral load curves for every episode. The transmission model is fitted once per draw and the fits are pooled, so the posterior reflects the uncertainty in when each person was infected.
- **New arguments in `prepare_stan_data()`:** `episode_timelines` / `episode_curves` (use Stage 1 timing), `viral_scale` (`"log10"` or `"Ct"`, chosen explicitly), and `risk_window` (see below). Without `episode_timelines`, the original random interval imputation is used, so existing scripts run unchanged.
- **At-risk windows:** each person is now at risk only while under observation, from their first to their last sample (`risk_window = "observed"`, the default), which handles late enrolment, early drop-out, and infections already present at enrolment. A reinfection can no longer be dated before the previous episode has resolved. Use `risk_window = "study"` to reproduce the previous behaviour (at risk for the whole study).
- **Simulator:**
  - `infectious_period = "viral"`: a person is infectious exactly while their viral load is detectable, matching the assumption of the two-stage estimator. The default `"gamma"` keeps the previous behaviour.
  - `hh_df` gains `infectious_start` and the **true infector** of every infection (`infected_by`, `infected_by_episode`), so reconstructed transmission chains can be validated.
  - Recorded test values now come from the same viral trajectory that drives transmission (previously they were evaluated one day ahead).
- **Fixes:** imputed log10 viral loads are floored at 0 (no negative values); `reconstruct_transmission_chains()` works when inert parameters have been dropped from the fit; the Stan models use the current array syntax (RStan ≥ 2.26).

---

## Authors

- **Ke Li** — Package creator and maintainer
- **Yiren Hou** — Package creator and maintainer

## Citation

If you use HHBayes in your work, please cite:
> Li K, Hou Y, Mukherjee B, Pitzer VE, Weinberger DM (2026). HHBayes: A Flexible Bayesian Framework for Simulating and Analyzing Household Transmission Dynamics. medRxiv. https://www.medrxiv.org/content/10.64898/2026.04.01.26349903v2

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("keli5734/HHBayes")
```

### Prerequisites

HHBayes depends on [RStan](https://mc-stan.org/rstan/) (version 2.26 or later), which requires a working C++ toolchain. See the [RStan Getting Started Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for platform-specific setup. Installing from source compiles two Stan models (the household transmission model and the Stage 1 viral kinetics model), which takes a few minutes.

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

**Stage 1 priors** (`prepare_kinetics_data(priors = ...)`) are set separately. Each element (`mu`, `mu_sd`, `tau_grp`, `tau_ep`) is a vector of length 5 in the order (level, log rise rate, log decay rate, peak day relative to the first positive, log time from infection to first positive), where level is `v_p` (log10) or `Cpeak` (Ct). Sensible defaults are provided for both scales.

---

## Quick Start

The recommended workflow has two stages: Stage 1 estimates viral kinetics and infection timing once; Stage 2 fits the household transmission model once per Stage 1 posterior draw, and the fits are pooled.

```r
library(HHBayes)
library(rstan)

# 0. Generate surveillance dataset
study_start <- "2024-07-01"
study_end   <- "2025-06-30"
dates_weekly <- seq(from = as.Date(study_start), to = as.Date(study_end), by = "week")
surveillance_data <- data.frame(
  date = dates_weekly,
  # Random epidemic curve (low start, peak middle, low end)
  cases = 0.1 + 100 * exp(-0.0002 * (as.numeric(dates_weekly - mean(dates_weekly)))^2) + abs(rnorm( length(dates_weekly),mean = 0, sd = 10))
)

# 1. Simulate 50 households; each person is infectious while their viral load is detectable
sim <- simulate_multiple_households_comm(
  n_households = 50,
  viral_testing = "viral load",
  infectious_period = "viral",
  surveillance_interval = 3,
  start_date = study_start,
  end_date   = study_end,
  seed = 123,
  surveillance_df = surveillance_data
)

# 2. Summarize attack rates and plot simulated data
rates <- summarize_attack_rates(sim)
rates$primary_by_role
my_plot <- plot_epidemic_curve(sim, surveillance_data, start_date_str = study_start, bin_width = 7)
print(my_plot)

df_for_stan <- sim$diagnostic_df
options(mc.cores = parallel::detectCores())

# 3. Stage 1: viral kinetics and infection timing (fitted once)
kin  <- prepare_kinetics_data(df_for_stan, study_start_date = as.Date(study_start),
                              viral_scale = "log10")
kfit <- fit_viral_kinetics(kin, iter = 2000, warmup = 1000, chains = 4, seed = 1)
summarize_kinetics(kfit)

# 4. Stage 2: household transmission, one fit per Stage 1 posterior draw
# NOTE: phi_role, kappa_role, beta1, beta2, alpha are sampled on the LOG scale,
# so dist = "normal" here == LogNormal on the biological quantity. See
# "A Note on Prior Scales" above.
my_priors <- list(
    beta1      = list(dist = "normal",  params = c(-5, 1)),
    beta2      = list(dist = "normal",  params = c(-5, 1)),
    alpha      = list(dist = "normal",  params = c(-7, 1)),
    phi_role   = list(dist = "normal",  params = c(0, 1)),
    kappa_role = list(dist = "normal",  params = c(0, 1)),
    vl_midpoint = list(dist = "normal", params = c(3, 0.5)),  # adjust if using Ct values
    vl_slope    = list(dist = "normal", params = c(2.5,  0.5))   # adjust if using Ct values
)

ex <- kinetics_extract(kfit)
M  <- 10                                                   # number of Stage 1 draws
draw_ids <- round(seq(1, dim(ex$tinc)[1], length.out = M))
fits <- vector("list", M); stan_inputs <- vector("list", M)
for (m in seq_len(M)) {
  tl <- kinetics_draw_to_timelines(kfit, draw_ids[m], ex)
  stan_inputs[[m]] <- prepare_stan_data(
      df_clean          = df_for_stan,
      surveillance_df   = surveillance_data,
      study_start_date  = as.Date(study_start),
      study_end_date    = as.Date(study_end),
      viral_scale       = "log10",
      use_vl_data       = TRUE,
      use_curve_logic   = FALSE,
      delta             = 0,
      priors            = my_priors,
      episode_timelines = tl$episode_timelines,
      episode_curves    = tl$episode_curves
  )
  fits[[m]] <- fit_household_model(stan_inputs[[m]], iter = 2000, warmup = 1000, chains = 4, seed = m)
}

# 5. Pool the M fits (the pooled posterior carries the infection-timing uncertainty)
fit <- rstan::sflist2stanfit(fits)
print(fit, pars = c("beta1", "beta2", "alpha_comm", "phi_by_role", "kappa_by_role"),
      probs = c(0.025, 0.5, 0.975))

# 6. Visualize
p_post <- plot_posterior_distributions(fit)
chains <- reconstruct_transmission_chains(fit = fit, stan_data = stan_inputs[[1]], min_prob_threshold =  .01)
selected_hh = 1
p_hh <- plot_household_timeline(chains, stan_inputs[[1]], target_hh_id = selected_hh) # define selected_hh.
print(p_hh)
```

> **Runtime:** Stage 2 is fitted `M` times. The fits are independent, so on a cluster they can run in parallel (one job per draw). See Step 10.

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
  infectious_period = "viral",       # infectious while viral load is detectable
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

**Infectious period.** With `infectious_period = "viral"`, each infected person is infectious exactly while their viral load is detectable (log10 VL ≥ `detect_threshold_log10`, or Ct ≤ `detect_threshold_Ct`), starting at least one day after infection, and both `beta1` and `beta2` act only on those days. This is the assumption made by the two-stage estimator, so use it to generate data for validating the estimation. With the default `infectious_period = "gamma"`, the infectious period is Gamma(`infectious_shape`, `infectious_scale`) days, as in previous versions.

The output is a list with two dataframes:

| Output | Description |
|---|---|
| `sim_res$hh_df` | One row per person per infection episode: household ID, person ID, role, infection time, start and end of infectiousness, recovery time, the true infector (`infected_by`: person ID in the same household, or `NA` for community infection; `infected_by_episode`), and covariate assignments. Times are integer study days. |
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

### Step 8: Stage 1 — Estimate Viral Kinetics and Infection Timing

The viral load (or Ct) measurements tell us more than which tests were positive: their shape over time shows where each infection was in its course. Stage 1 fits a hierarchical model of viral trajectories by age group and, for every infection episode, the time from infection to the first positive test. It is fitted **once**.

```r
kin <- prepare_kinetics_data(
  df_for_stan,
  study_start_date = as.Date(study_start),
  viral_scale      = "log10"            # or "Ct"
)
kfit <- fit_viral_kinetics(kin, iter = 2000, warmup = 1000, chains = 4, seed = 1)

summarize_kinetics(kfit)                # group-level kinetics and time to first positive
```

`prepare_kinetics_data()` collects, for each episode, its positive samples plus the negative samples up to 21 days before and 28 days after (never crossing into a neighbouring episode of the same person). Negatives enter as censored observations ("below the detection limit"), so they constrain when the viral load rose and fell. Episodes with fewer than 3 samples or 2 positives are not fitted individually; they receive their age group's typical trajectory. Optional symptom-onset dates (`symptom_onset`) add information on timing.

**From Stage 1 to infection timelines.** Each Stage 1 posterior draw gives, for every episode:

- **Infection day:** the first positive test minus the estimated time to first positive (at least `min_incub_days = 1` day).
- **Infectious window:** the days on which the fitted viral load is detectable, starting at least one day after infection, never before the day after the last negative test or after the first positive, and never ending before the last positive or on or after the next negative.
- **Viral load curve** on the infectious days (observed values replace the curve on tested days).

```r
ex <- kinetics_extract(kfit)
tl <- kinetics_draw_to_timelines(kfit, draw = 1, ex = ex)
head(tl$episode_timelines)   # participantid, episode_id, date_infection,
                             # date_infectious_start, date_infectious_end
head(tl$episode_curves)      # participantid, episode_id, date, value
```

> **What Stage 1 can and cannot identify.** The viral load data pin down the shape of each trajectory and its timing relative to the first positive test well. The time from infection to the first positive is only weakly identified by viral load data alone and leans on its prior (default median about 3 days for log10 data, 4 days for Ct). For real data, set this prior from the literature (`priors$mu[5]`, on the log scale) and supply symptom-onset dates if available.

### Step 9: Prepare Data for Stan

This step transforms the diagnostic dataframe into the structured list that Stan expects. It handles column renaming, infection and infectious windows (from Stage 1, or by random imputation), viral load filling, at-risk windows, covariate matrix construction, and prior specification.

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

Supported distributions depend on the parameter's scale (see "A Note on Prior
Scales" above):

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
  vl_midpoint = list(dist = "normal", params = c(3, 0.5)),  # adjust if Ct values
  vl_slope    = list(dist = "normal", params = c(2.5,  0.5))   # adjust if Ct values
)
```

**Build the Stan input from one Stage 1 draw:**

```r
tl <- kinetics_draw_to_timelines(kfit, draw = 1, ex = ex)

stan_input <- prepare_stan_data(
  df_clean          = df_for_stan,
  surveillance_df   = surveillance_data,
  study_start_date  = as.Date(study_start),
  study_end_date    = as.Date(study_end),
  viral_scale       = "log10",            # or "Ct"; set explicitly
  use_vl_data       = TRUE,
  use_curve_logic   = FALSE,
  delta = 0,
  priors            = my_priors,
  role_mixing_matrix = role_mixing_weights,
  covariates_susceptibility = "vacc_status",
  covariates_infectivity    = "vacc_status",
  episode_timelines = tl$episode_timelines,
  episode_curves    = tl$episode_curves)
```

**At-risk windows.** By default (`risk_window = "observed"`) each person is at risk from their first to their last sample, so late enrolment and early drop-out are handled, and an infection already present at the first sample is not counted as an event. After an infection, a person is at risk again only once the previous episode has resolved (end of infectiousness plus an immunity period drawn from `recovery_params`). Use `risk_window = "study"` to keep everyone at risk for the whole study, as in previous versions.

**Viral scale.** Set `viral_scale = "log10"` or `"Ct"` explicitly. With `"auto"`, the scale is guessed from the data (values above 15 are taken as Ct) and a message is printed.

### Step 10: Fit the Bayesian Model

`fit_household_model()` runs HMC sampling via `rstan::sampling()`. By default
(`drop_inert = TRUE`) it automatically excludes internal bookkeeping quantities
**and** any viral-component parameters that are inert under your data
configuration, so the printed fit shows only the parameters that were actually
informed by the likelihood.

Fit the model once per Stage 1 draw, then pool the fits:

```r
options(mc.cores = parallel::detectCores())

ex <- kinetics_extract(kfit)
M  <- 10
draw_ids <- round(seq(1, dim(ex$tinc)[1], length.out = M))
fits <- vector("list", M); stan_inputs <- vector("list", M)

for (m in seq_len(M)) {
  tl <- kinetics_draw_to_timelines(kfit, draw_ids[m], ex)
  stan_inputs[[m]] <- prepare_stan_data(
    df_clean = df_for_stan, surveillance_df = surveillance_data,
    study_start_date = as.Date(study_start), study_end_date = as.Date(study_end),
    viral_scale = "log10", use_vl_data = TRUE, use_curve_logic = FALSE, delta = 0,
    priors = my_priors, role_mixing_matrix = role_mixing_weights,
    covariates_susceptibility = "vacc_status", covariates_infectivity = "vacc_status",
    episode_timelines = tl$episode_timelines, episode_curves = tl$episode_curves)

  fits[[m]] <- fit_household_model(
    stan_inputs[[m]],
    iter   = 2000,
    warmup = 1000,
    chains = 4,
    seed   = m,
    cores  = 4
  )
}

# Pooled posterior: all draws from the M fits together
fit <- rstan::sflist2stanfit(fits)
print(fit, probs = c(0.025, 0.5, 0.975))
```

- **How many draws?** `M = 10`–`20` is usually enough; more draws give a smoother pooled posterior.
- **Convergence:** check Rhat and effective sample size on each of the `M` fits. Do not read Rhat from the pooled object: it mixes fits to different timelines, and the spread between them inflates Rhat by design.
- **Running in parallel:** the `M` fits are independent. On a cluster, run one job per draw (save `tl` from each Stage 1 draw with `saveRDS()`), then pool.

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

**Without Stage 1.** If your data have only positive/negative results, or you want a quick fit, call `prepare_stan_data()` without `episode_timelines`. Infection days are then drawn uniformly between the last negative and the first positive test, the end of infectiousness between the last positive and the next negative test, and viral loads on untested days come from the role-specific curves in `imputation_params`. This is a single imputation, so the uncertainty in infection timing is not propagated.

### Step 11: Visualize Results

```r
# Posterior distributions of susceptibility (phi) and infectivity (kappa) by role
plot_posterior_distributions(fit)

# Reconstruct who-infected-whom transmission chains
chains <- reconstruct_transmission_chains(fit, stan_inputs[[1]], min_prob_threshold = 0.001)

# Plot a specific household's infection timeline with transmission arrows
plot_household_timeline(chains, stan_inputs[[1]], target_hh_id = 11)
```

A household member is a candidate source only on days inside their infectious window. With simulated data, compare the reconstructed sources with the true infectors in `sim_res$hh_df$infected_by`.

> **Interpreting role multipliers.** Susceptibility (`phi`) of one group and infectivity (`kappa`) of another are separated only through variation in who lives with whom. When most households have the same composition, the data mainly determine their products (how likely an infected person of one group is to infect a person of another), and the individual multipliers and `beta1` are less well identified. Households of varied composition, and transmission that is mostly viral-load driven (`beta1` small relative to `beta2`), make them easier to separate.

---

## Simulation Parameter Reference

### Transmission

| Parameter | Default | Description |
|---|---|---|
| `beta1` | `0.008` | Baseline within-household transmission rate (contact-driven). |
| `beta2` | `0.008` | Viral-load-dependent transmission rate. Total per-infector force = `beta1 + beta2 * f(VL)`. |
| `delta` | `0` | Household size scaling. Force scaled by `(1/max(n_h,1))^delta` in larger households. |
| `alpha_comm_by_role` | `5e-4` | Daily community acquisition rate. |
| `phi_by_role` | `c(adult=1, infant=4, toddler=5, elderly=1)` | Susceptibility multipliers by role. |
| `kappa_by_role` | `c(adult=1, infant=1, toddler=1.2, elderly=1)` | Infectivity multipliers by role. |

### Infectious Period and Immunity

| Parameter | Default | Description |
|---|---|---|
| `infectious_period` | `"gamma"` | `"gamma"`: Gamma-distributed infectious period. `"viral"`: infectious exactly while the viral load is detectable (matches the two-stage estimator). |
| `infectious_shape` | `3` | Gamma shape for infectious duration (`infectious_period = "gamma"` only). |
| `infectious_scale` | `1` | Gamma scale for infectious duration. Mean = shape x scale. |
| `max_shedding_days` | `60` | Days over which the viral trajectory is followed when `infectious_period = "viral"`. |
| `waning_shape` | `16` | Gamma shape for immunity waning. |
| `waning_scale` | `10` | Gamma scale for immunity waning. Mean = 160 days. |
| `max_infections` | `Inf` | Max infections per person. Set to `1` to disable reinfection. |

### Viral Dynamics

| Parameter | Default | Description |
|---|---|---|
| `viral_testing` | `"viral load"` | `"viral load"` (log10 scale) or `"Ct"` (cycle threshold). |
| `model_type` | `"empirical"` | `"empirical"` (parametric curves) or `"ODE"` (within-host ODE). |
| `VL_params_list` / `Ct_params_list` | built-in | Role-specific trajectory parameters (`v_p`, `t_p`, `lambda_g`, `lambda_d`, or `Cpeak`, `r`, `d`, `t_peak`). |
| `V_ref` / `Ct_50` | `3.0` / `40` | Reference log10 VL for infectiousness scaling, or Ct at 50% infectiousness (sigmoid). Estimated as `vl_midpoint`. |
| `V_rho` / `Ct_delta` | `2.5` / `2` | Power-law exponent for log10 VL scaling, or steepness of the Ct-infectiousness sigmoid. Estimated as `vl_slope`. |

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
| `viral_scale` | `"auto"` | `"log10"` or `"Ct"`. With `"auto"` the scale is guessed from the data (values above 15 are treated as Ct). Set it explicitly. |
| `episode_timelines` | `NULL` | Infection dates and infectious windows per episode (from `kinetics_draw_to_timelines()`). `NULL` = random interval imputation. |
| `episode_curves` | `NULL` | Daily viral load / Ct per episode (from `kinetics_draw_to_timelines()`). `NULL` = role curves in `imputation_params`. |
| `risk_window` | `"observed"` | `"observed"`: at risk from first to last sample. `"study"`: at risk for the whole study (previous behaviour). |
| `recovery_params` | Gamma(2, 3) | Immunity period after an episode, per role (`list(role = list(shape, scale))`); a reinfection is at risk only after it. |
| `imputation_params` | built-in | Role-specific viral curves used to fill viral loads when `episode_curves` is not supplied. |
| `use_vl_data` | `TRUE` | Use actual viral load data vs. generation interval curve. |
| `use_curve_logic` | `TRUE` | When `use_vl_data = FALSE`, estimate generation interval curve vs. constant infectivity. |
| `delta` | `0` | Household size scaling. Force scaled by `(1/max(n_h,1))^delta`. |

### Stage 1 Viral Kinetics (`prepare_kinetics_data()`)

| Parameter | Default | Description |
|---|---|---|
| `viral_scale` | `"auto"` | `"log10"` or `"Ct"`. |
| `lod` | `1e-6` / `40` | Detection limit (log10 VL / Ct); negative tests are censored at this value. |
| `groups` | role | Grouping for the hierarchy (default: the person's role). |
| `pre_window`, `post_window` | `21`, `28` | Days of negative tests used before the first and after the last positive. |
| `min_samples`, `min_positive` | `3`, `2` | Episodes with fewer samples are not fitted individually and use their group's trajectory. |
| `symptom_onset` | `NULL` | Optional symptom-onset dates (adds information on infection timing). |
| `nu` | `4` | Degrees of freedom of the Student-t measurement error (0 = normal). |
| `tinc_bounds` | `c(1, 21)` | Plausible range (days) from infection to the first positive test. |
| `peak_bounds` | `c(1, 14)` | Plausible range (days) from infection to the viral peak. |
| `priors` | scale-specific | Population means and SDs; see "A Note on Prior Scales". |

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
| `prepare_kinetics_data()` | Format viral load / Ct data for the Stage 1 kinetics model |
| `fit_viral_kinetics()` | Fit the Stage 1 viral kinetics and infection-timing model |
| `summarize_kinetics()` | Group-level kinetics and time to first positive |
| `kinetics_extract()` | Extract Stage 1 posterior draws |
| `kinetics_draw_to_timelines()` | Infection dates, infectious windows and viral curves from one Stage 1 draw |
| `kinetics_curve()` | Evaluate a viral load / Ct trajectory |
| `infectious_window_from_curve()` | Infection day and infectious window for one episode |
| `prepare_stan_data()` | Format data for the transmission model |
| `fit_household_model()` | Fit the Bayesian transmission model via HMC |
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


Individuals progress through **S -> I -> R -> S**. In the simulator the infectious period is Gamma-distributed (`infectious_period = "gamma"`) or lasts while the viral load is detectable (`"viral"`), and immunity wanes after a Gamma-distributed period. In the transmission model, person $i$ contributes to the force of infection only on the days of their infectious window, and is at risk from their first to their last sample (or until infection). Inference is performed via HMC in Stan (`adapt_delta = 0.95`, `max_treedepth = 15`).

### Stage 1: Viral Kinetics and Infection Timing

For infection episode $e$ in age group $g$, the expected viral load at $s$ days after infection follows a double exponential (log10 viral load)

$$V_e(s) = \log_{10}\left(\frac{2 \cdot 10^{v_{p,e}}}{\exp(-\lambda_{g,e}(s - t_{p,e})) + \exp(\lambda_{d,e}(s - t_{p,e}))}\right)$$

or a piecewise-linear "tent" (Ct): $C_e(s) = C_{\text{peak},e} + r_e (t_{\text{peak},e} - s)$ before the peak and $C_{\text{peak},e} + d_e (s - t_{\text{peak},e})$ after. Days since infection are $s = t - t^{\text{first pos}}_e + t^{\text{inc}}_e$, where $t^{\text{inc}}_e$ is the time from infection to the first positive test.

The episode parameters (level, log rise rate, log decay rate, peak day relative to the first positive, $\log t^{\text{inc}}$) are hierarchical: population mean, plus a group deviation, plus an episode deviation (non-centred). Peak timing is sampled relative to the first positive test, $t_{\text{peak}} = c + t^{\text{inc}}$, because the data identify $c$ directly. Positive samples follow a Student-t measurement model; negative samples are censored at the detection limit. Soft bounds keep $t^{\text{inc}}$, the time to peak, and the end of detectability in plausible ranges.

**Combining the stages.** For each of $M$ Stage 1 posterior draws $\mathcal{T}^{(m)}$ (infection days, infectious windows, viral curves), the transmission model is fitted, and the draws are pooled:

$$p(\theta \mid y) \approx \frac{1}{M} \sum_{m=1}^{M} p\left(\theta \mid y, \mathcal{T}^{(m)}\right)$$

Information flows one way: the viral data inform the timelines, but the household transmission data do not feed back into Stage 1.

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

### Stage 1 (Viral Kinetics) Parameters

| Parameter | Description |
|---|---|
| `v_p`, `t_p`, `lambda_g`, `lambda_d` | log10 trajectory by group: peak level, peak time, rise and decay rates |
| `Cpeak`, `r`, `d`, `t_peak` | Ct trajectory by group: peak Ct, rise and decay slopes, peak time |
| `t_incub` | Time from infection to the first positive test, by group and episode |
| `sigma_y` | Measurement error of viral load / Ct |

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

