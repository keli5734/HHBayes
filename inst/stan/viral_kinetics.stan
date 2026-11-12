// =============================================================================
// viral_kinetics.stan  --  "Stage 1": within-host viral kinetics + infection timing
//
// Estimates, for every episode, (i) the viral trajectory parameters and (ii) the
// infection time, jointly. The posterior draws feed Stage 2 (household_transmission)
// through episode timelines (see R/viral_kinetics.R).
//
// SCALE (user choice, vl_type):
//   vl_type = 0  Ct values. Mean curve is the piecewise-linear "tent"
//                  Ct(s) = Cpeak + r (t_peak - s)   for s <= t_peak   (rise)
//                        = Cpeak + d (s - t_peak)   for s >  t_peak   (decay)
//                parameters cp = (Cpeak, r, d, t_peak)
//                a sample is NEGATIVE when the latent Ct is >= lod  (right-censored)
//   vl_type = 1  log10 viral load. Mean curve is the double exponential
//                  log10 VL(s) = log10( 2 * 10^v_p / (exp(-lambda_g (s - t_p)) + exp(lambda_d (s - t_p))) )
//                parameters cp = (v_p, t_p, lambda_g, lambda_d)
//                a sample is NEGATIVE when the latent value is < lod (left-censored)
//   s = days since infection = t_rel + tinc, where t_rel = sample day - first positive
//   day of the episode, and tinc >= 0 is the time from infection to the first
//   positive sample. tinc is the quantity that makes onset timing explicit.
//
// HOW THE NEGATIVES INFORM ONSET: samples taken before the first positive and after
// the last positive enter as CENSORED observations. A negative before the first
// positive says the curve was still undetectable then, which bounds how early the
// curve can rise relative to the first positive, and hence tinc. Negatives after the
// last positive bound the decay (end of the detectable period).
//
// HIERARCHY (non-centred): population -> group (e.g. age group) -> episode, for
//   theta = (level, log rise rate, log decay rate, c, log tinc)
// where level = Cpeak (Ct) or v_p (log10 VL), rise/decay are r,d (Ct) or lambda_g,lambda_d
// (log10 VL), and c is the peak day RELATIVE TO THE FIRST POSITIVE sample. The time from
// infection to peak is derived, t_peak = c + tinc. Sampling c (which the data identify)
// instead of t_peak removes the ridge between t_peak and tinc that otherwise makes NUTS
// take very long trajectories: the data fix t_peak - tinc, and only the priors/bounds
// separate the two. tinc is log-normal with soft truncation.
// =============================================================================

functions {
  // Gaussian-tail soft bound: 0 inside [lo, hi], quadratic penalty outside.
  real soft_bound_lp(real x, real lo, real hi, real s) {
    real pen = 0;
    if (x < lo) {
      pen += -0.5 * square((lo - x) / s);
    }
    if (x > hi) {
      pen += -0.5 * square((x - hi) / s);
    }
    return pen;
  }
  
  // Natural-scale curve parameters from the sampled quantities.
  //   Ct:        cp = (Cpeak, r, d, t_peak)
  //   log10 VL:  cp = (v_p, t_p, lambda_g, lambda_d)
  // with t_peak (t_p) = c + tinc, days since infection.
  row_vector natural_params(int vl_type, real level, real log_rise,
                            real log_decay, real c, real tinc) {
    row_vector[4] cp;
    if (vl_type == 0) {
      cp[1] = level;
      cp[2] = exp(log_rise);
      cp[3] = exp(log_decay);
      cp[4] = c + tinc;
    } else {
      cp[1] = level;
      cp[2] = c + tinc;
      cp[3] = exp(log_rise);
      cp[4] = exp(log_decay);
    }
    return cp;
  }
  
  // Expected value on the native scale at s days since infection.
  real mean_curve(real s, int vl_type, row_vector cp) {
    if (vl_type == 0) {
      if (s <= cp[4]) {
        return cp[1] + cp[2] * (cp[4] - s);
      }
      return cp[1] + cp[3] * (s - cp[4]);
    }
    return log10(2) + cp[1]
           - log_sum_exp(-cp[3] * (s - cp[2]), cp[4] * (s - cp[2])) / log(10);
  }
}
data {
  int<lower=1> N_ep; // episodes
  int<lower=1> N_obs; // samples across all episodes
  int<lower=1> G; // groups (e.g. infant / toddler / adult)
  
  int<lower=0, upper=1> vl_type; // 1 = log10 VL, 0 = Ct
  real lod; // detection threshold on the native scale
  
  array[N_obs] int<lower=1, upper=N_ep> ep; // episode of each sample
  array[N_obs] real t_rel; // day of sample minus day of the episode's first positive
  array[N_obs] real y; // measured value (ignored when cens = 1)
  array[N_obs] int<lower=0, upper=1> cens; // 1 = negative sample (censored), 0 = quantified positive
  
  array[N_ep] int<lower=1, upper=G> grp; // group of each episode
  array[N_ep] real<lower=0> tinc_sd_mult; // 1 normally; >1 when no earlier negative exists (wider prior)
  
  // optional symptom-onset information: (tinc - sx_lag) ~ N(sx_mean, sx_sd),
  // sx_lag = first positive sample day minus symptom-onset day
  array[N_ep] int<lower=0, upper=1> has_sx;
  array[N_ep] real sx_lag;
  real sx_mean;
  real<lower=0> sx_sd;
  
  // observation model
  real nu; // Student-t df; nu <= 0 selects a normal likelihood
  real<lower=0> sigma_floor; // lower bound on the residual SD
  real<lower=0> prior_sigma_y; // half-normal scale for the residual SD
  
  // priors on the population means of theta = (level, log rise, log decay, c, log tinc)
  vector[5] prior_mu;
  vector<lower=0>[5] prior_mu_sd;
  vector<lower=0>[5] prior_tau_grp; // half-normal scales, group-level SDs (per parameter)
  vector<lower=0>[5] prior_tau_ep; // half-normal scales, episode-level SDs (per parameter)
  
  // soft biological bounds (days)
  vector[2] tinc_bounds; // infection -> first positive sample
  vector[2] peak_bounds; // infection -> peak
  real tend_max; // infection -> end of detectable period (upper)
  real<lower=0> tend_soft; // softness of the tend_max bound
  real<lower=0> bound_soft; // softness of the tinc / peak bounds
}
parameters {
  vector[5] mu_pop;
  vector<lower=0>[5] tau_grp;
  vector<lower=0>[5] tau_ep;
  matrix[G, 5] z_grp;
  matrix[N_ep, 5] z_ep;
  real<lower=sigma_floor> sigma_y;
}
transformed parameters {
  matrix[N_ep, 5] theta;
  matrix[N_ep, 4] cp; // natural-scale curve parameters of each episode (see natural_params)
  vector[N_ep] tinc; // days from infection to first positive sample
  for (i in 1 : N_ep) {
    for (k in 1 : 5) {
      real sd_mult = (k == 5) ? tinc_sd_mult[i] : 1.0;
      theta[i, k] = mu_pop[k] + tau_grp[k] * z_grp[grp[i], k]
                    + tau_ep[k] * sd_mult * z_ep[i, k];
    }
    tinc[i] = exp(theta[i, 5]);
    cp[i] = natural_params(vl_type, theta[i, 1], theta[i, 2], theta[i, 3],
                           theta[i, 4], tinc[i]);
  }
}
model {
  int pk = (vl_type == 1) ? 2 : 4; // column of cp holding the time-to-peak parameter
  
  // ---- hyperpriors ----
  mu_pop ~ normal(prior_mu, prior_mu_sd);
  tau_grp ~ normal(0, prior_tau_grp); // tau_* are declared <lower=0>: half-normal
  tau_ep ~ normal(0, prior_tau_ep);
  to_vector(z_grp) ~ std_normal();
  to_vector(z_ep) ~ std_normal();
  sigma_y ~ normal(0, prior_sigma_y);
  
  // ---- biological constraints, per episode ----
  for (i in 1 : N_ep) {
    real s_end; // days from infection until the curve returns to the LOD
    target += soft_bound_lp(tinc[i], tinc_bounds[1], tinc_bounds[2],
                            bound_soft);
    target += soft_bound_lp(cp[i, pk], peak_bounds[1], peak_bounds[2],
                            bound_soft);
    
    if (vl_type == 0) {
      s_end = cp[i, 4] + (lod - cp[i, 1]) / cp[i, 3];
    } else {
      // asymptotic crossing time of the decay branch (good enough for a soft bound)
      s_end = cp[i, 2] + (cp[i, 1] + log10(2) - lod) * log(10) / cp[i, 4];
    }
    if (s_end > tend_max) {
      target += -0.5 * square((s_end - tend_max) / tend_soft);
    }
    
    if (has_sx[i] == 1) {
      target += normal_lpdf(tinc[i] - sx_lag[i] | sx_mean, sx_sd);
    }
  }
  
  // ---- likelihood: quantified positives + censored negatives ----
  for (k in 1 : N_obs) {
    int i = ep[k];
    real mu = mean_curve(t_rel[k] + tinc[i], vl_type, cp[i]);
    if (cens[k] == 0) {
      if (nu > 0) {
        y[k] ~ student_t(nu, mu, sigma_y);
      } else {
        y[k] ~ normal(mu, sigma_y);
      }
    } else if (vl_type == 0) {
      // negative <=> latent Ct >= lod
      if (nu > 0) {
        target += student_t_lccdf(lod | nu, mu, sigma_y);
      } else {
        target += normal_lccdf(lod | mu, sigma_y);
      }
    } else // negative <=> latent log10 VL < lod
    if (nu > 0) {
      target += student_t_lcdf(lod | nu, mu, sigma_y);
    } else {
      target += normal_lcdf(lod | mu, sigma_y);
    }
  }
}
generated quantities {
  // Group-level typical trajectory and incubation (used for episodes that have no
  // Stage 1 fit of their own, e.g. too few samples).
  matrix[G, 4] grp_curve;
  vector[G] grp_tinc;
  for (g in 1 : G) {
    grp_tinc[g] = exp(mu_pop[5] + tau_grp[5] * z_grp[g, 5]);
    grp_curve[g] = natural_params(vl_type,
                                  mu_pop[1] + tau_grp[1] * z_grp[g, 1],
                                  mu_pop[2] + tau_grp[2] * z_grp[g, 2],
                                  mu_pop[3] + tau_grp[3] * z_grp[g, 3],
                                  mu_pop[4] + tau_grp[4] * z_grp[g, 4],
                                  grp_tinc[g]);
  }
}

