data {
  int<lower=1> N; // Total number of people
  int<lower=1> T; // Total days
  int<lower=1> H; // Total households
  int<lower=1> R; // Number of Roles (4: Adult, Infant, Toddler, Elderly)
  real delta;     // Density scaling parameter (usually fixed or passed as data)

  // --- FLAGS ---
  int<lower=0, upper=1> use_vl_data;     // 1 = Use Viral Load, 0 = Use simple curve
  int<lower=0, upper=1> vl_type;         // 1 = Log10, 0 = Ct
  int<lower=0, upper=1> use_curve_logic; // 1 = Use Gamma Curve Fallback

  // --- INDEXING ---
  int<lower=1, upper=H> hh_id[N];   // Household ID for each person
  int<lower=1, upper=R> role_id[N]; // Role ID for each person

  // --- OUTCOMES ---
  int<lower=0, upper=1> I[N, T];        // 1 if infected on day t (Event)
  int<lower=0, upper=1> Y[N, T];        // 1 if infectious on day t (State)
  real V[N, T];                         // Viral Load value on day t
  int<lower=1, upper=T> start_risk[N];  // Day when risk starts for person n

  // --- HOUSEHOLD STRUCTURE ---
  int<lower=1> hh_size_people[H];

  // --- SPARSE CONTACT MATRIX (The Speed Fix) ---
  int<lower=0> n_contacts;                       // Total number of directed edges
  int<lower=1, upper=N> contact_src[n_contacts]; // Infector ID (Source)
  int<lower=1, upper=N> contact_tgt[n_contacts]; // Susceptible ID (Target)
  real<lower=0> contact_w[n_contacts];           // Contact Weight (0 to 1)

  // --- FORCING & REFERENCE ---
  matrix[T, R] seasonal_forcing_mat;
  real<lower=0> reference_phi;   // Reference susceptibility (Adult=1.0)
  real<lower=0> reference_kappa; // Reference infectivity (Adult=1.0)

  // --- COVARIATES ---
  int<lower=0> K_susc;
  matrix[N, K_susc] X_susc;
  int<lower=0> K_inf;
  matrix[N, K_inf] X_inf;

  // --- PRIORS (Flexible Distribution Types) ---
  // 1=Normal, 2=Uniform, 3=LogNormal
  int<lower=0> prior_beta1_type; vector[2] prior_beta1_params;
  int<lower=0> prior_beta2_type; vector[2] prior_beta2_params;
  int<lower=0> prior_alpha_type; vector[2] prior_alpha_params;
  int<lower=0> prior_cov_type;   vector[2] prior_cov_params;

  // Role-specific priors
  int<lower=0> prior_phi_type;   vector[2] prior_phi_params;
  int<lower=0> prior_kappa_type; vector[2] prior_kappa_params;

  // Viral Dynamics Priors
  int<lower=0> prior_shape_type; vector[2] prior_shape_params;
  int<lower=0> prior_rate_type;  vector[2] prior_rate_params;
  int<lower=0> prior_ct50_type;  vector[2] prior_ct50_params;
  int<lower=0> prior_slope_type; vector[2] prior_slope_params;
  int<lower=0> prior_vref_type;  vector[2] prior_vref_params;
  int<lower=0> prior_rho_type;   vector[2] prior_rho_params;
}

transformed data {
  // Pre-calculate the exact day of infection for efficiency
  // 0 means "Never infected"
  int infection_day[N];
  for (n in 1:N) {
    infection_day[n] = 0;
    for (t in start_risk[n]:T) {
      if (I[n, t] == 1) {
        infection_day[n] = t;
        break;
      }
    }
  }

  // Calculate conditional parameter counts for array sizing
  int n_gen_params = use_curve_logic;
  int n_ct_params = (use_vl_data == 1 && vl_type == 0) ? 1 : 0;
  int n_vl_params = (use_vl_data == 1 && vl_type == 1) ? 1 : 0;
}

parameters {
  // Role Effects (Relative to Reference) - ALWAYS ESTIMATED
  vector[R-1] log_phi_by_role_raw;
  vector[R-1] log_kappa_by_role_raw;

  // Transmission Parameters - ALWAYS ESTIMATED
  real log_beta1;
  real log_beta2;
  real log_alpha_comm;

  // Covariate Effects - ALWAYS DECLARED (size based on data)
  vector[K_susc] beta_susc;
  vector[K_inf] beta_inf;

  // CONDITIONAL VIRAL DYNAMICS PARAMETERS

  // Generation interval parameters (only when use_curve_logic == 1)
  vector<lower=1.0, upper=20.0>[n_gen_params] gen_shape_raw;
  vector<lower=0.1, upper=5.0>[n_gen_params] gen_rate_raw;

  // Ct parameters (only when use_vl_data == 1 AND vl_type == 0)
  vector<lower=0>[n_ct_params] Ct50_raw;
  vector<lower=0>[n_ct_params] slope_ct_raw;

  // Log10 VL parameters (only when use_vl_data == 1 AND vl_type == 1)
  vector<lower=0>[n_vl_params] V_ref_raw;
  vector<lower=0>[n_vl_params] rho_raw;
}

transformed parameters {
  // Expand parameters to natural scale
  vector<lower=0>[R] phi_by_role;
  vector<lower=0>[R] kappa_by_role;
  vector[T] g_curve_est;
  real<lower=0> alpha_comm = exp(log_alpha_comm);
  real beta1 = exp(log_beta1);
  real beta2 = exp(log_beta2);

  // Extract scalar viral parameters (with defaults for unused cases)
  real gen_shape = n_gen_params > 0 ? gen_shape_raw[1] : 3.0;
  real gen_rate = n_gen_params > 0 ? gen_rate_raw[1] : 1.0;
  real Ct50 = n_ct_params > 0 ? Ct50_raw[1] : 35.0;
  real slope_ct = n_ct_params > 0 ? slope_ct_raw[1] : 2.0;
  real V_ref = n_vl_params > 0 ? V_ref_raw[1] : 3.0;
  real rho = n_vl_params > 0 ? rho_raw[1] : 2.5;

  // 1. Construct Age-Specific Susceptibility (Phi)
  phi_by_role[1] = reference_phi;
  for (r in 2:R) phi_by_role[r] = reference_phi * exp(log_phi_by_role_raw[r-1]);

  // 2. Construct Age-Specific Infectivity (Kappa)
  kappa_by_role[1] = reference_kappa;
  for (r in 2:R) kappa_by_role[r] = reference_kappa * exp(log_kappa_by_role_raw[r-1]);

  // 3. Pre-calculate Gamma Curve (only when use_curve_logic == 1)
  if (use_curve_logic == 1) {
    vector[T] raw_curve;
    for(d in 1:T) raw_curve[d] = exp(gamma_lpdf(d | gen_shape, gen_rate));
    g_curve_est = raw_curve / sum(raw_curve);
  } else {
    g_curve_est = rep_vector(0.0, T); // Dummy values when not used
  }

  // 4. Pre-calculate Viral Term (V_term) for the whole matrix
  matrix[N, T] V_term_calc = rep_matrix(0.0, N, T);
  if (use_vl_data == 1) {
    for (n in 1:N) {
      for (t in 1:T) {
        if (Y[n, t] == 1) {
          real val = V[n, t];

          // Log10 viral load case (vl_type == 1)
          if (vl_type == 1) {
            V_term_calc[n, t] = pow(fmax(0.0, val) / V_ref, rho);
          }
          // Ct case (vl_type == 0)
          else {
            V_term_calc[n, t] = inv_logit((Ct50 - val) / slope_ct);
          }
        }
      }
    }
  }
}

model {
  // =========================================================
  // 1. PRIORS
  // =========================================================

  // Role-specific priors (flexible)
  if (prior_phi_type == 1) log_phi_by_role_raw ~ normal(prior_phi_params[1], prior_phi_params[2]);
  else if (prior_phi_type == 2) log_phi_by_role_raw ~ uniform(prior_phi_params[1], prior_phi_params[2]);
  else if (prior_phi_type == 3) log_phi_by_role_raw ~ lognormal(prior_phi_params[1], prior_phi_params[2]);

  if (prior_kappa_type == 1) log_kappa_by_role_raw ~ normal(prior_kappa_params[1], prior_kappa_params[2]);
  else if (prior_kappa_type == 2) log_kappa_by_role_raw ~ uniform(prior_kappa_params[1], prior_kappa_params[2]);
  else if (prior_kappa_type == 3) log_kappa_by_role_raw ~ lognormal(prior_kappa_params[1], prior_kappa_params[2]);

  // Transmission Priors
  if (prior_beta1_type == 1) log_beta1 ~ normal(prior_beta1_params[1], prior_beta1_params[2]);
  else if (prior_beta1_type == 2) log_beta1 ~ uniform(prior_beta1_params[1], prior_beta1_params[2]);

  if (prior_beta2_type == 1) log_beta2 ~ normal(prior_beta2_params[1], prior_beta2_params[2]);
  else if (prior_beta2_type == 2) log_beta2 ~ uniform(prior_beta2_params[1], prior_beta2_params[2]);

  if (prior_alpha_type == 1) log_alpha_comm ~ normal(prior_alpha_params[1], prior_alpha_params[2]);
  else if (prior_alpha_type == 2) log_alpha_comm ~ uniform(prior_alpha_params[1], prior_alpha_params[2]);

  // Covariate Priors
  if (K_susc > 0) {
    if (prior_cov_type == 1) beta_susc ~ normal(prior_cov_params[1], prior_cov_params[2]);
    else beta_susc ~ uniform(prior_cov_params[1], prior_cov_params[2]);
  }
  if (K_inf > 0) {
    if (prior_cov_type == 1) beta_inf ~ normal(prior_cov_params[1], prior_cov_params[2]);
    else beta_inf ~ uniform(prior_cov_params[1], prior_cov_params[2]);
  }

  // CONDITIONAL VIRAL DYNAMICS PRIORS (only applied when parameters exist)

  // Generation interval priors
  if (n_gen_params > 0) {
    if (prior_shape_type == 1) gen_shape_raw ~ normal(prior_shape_params[1], prior_shape_params[2]);
    else if (prior_shape_type == 2) gen_shape_raw ~ uniform(prior_shape_params[1], prior_shape_params[2]);
    else if (prior_shape_type == 3) gen_shape_raw ~ lognormal(prior_shape_params[1], prior_shape_params[2]);

    if (prior_rate_type == 1) gen_rate_raw ~ normal(prior_rate_params[1], prior_rate_params[2]);
    else if (prior_rate_type == 2) gen_rate_raw ~ uniform(prior_rate_params[1], prior_rate_params[2]);
    else if (prior_rate_type == 3) gen_rate_raw ~ lognormal(prior_rate_params[1], prior_rate_params[2]);
  }

  // Ct parameters priors
  if (n_ct_params > 0) {
    if (prior_ct50_type == 1) Ct50_raw ~ normal(prior_ct50_params[1], prior_ct50_params[2]);
    else if (prior_ct50_type == 2) Ct50_raw ~ uniform(prior_ct50_params[1], prior_ct50_params[2]);
    else if (prior_ct50_type == 3) Ct50_raw ~ lognormal(prior_ct50_params[1], prior_ct50_params[2]);

    if (prior_slope_type == 1) slope_ct_raw ~ normal(prior_slope_params[1], prior_slope_params[2]);
    else if (prior_slope_type == 2) slope_ct_raw ~ uniform(prior_slope_params[1], prior_slope_params[2]);
    else if (prior_slope_type == 3) slope_ct_raw ~ lognormal(prior_slope_params[1], prior_slope_params[2]);
  }

  // Log10 VL parameters priors
  if (n_vl_params > 0) {
    if (prior_vref_type == 1) V_ref_raw ~ normal(prior_vref_params[1], prior_vref_params[2]);
    else if (prior_vref_type == 2) V_ref_raw ~ uniform(prior_vref_params[1], prior_vref_params[2]);
    else if (prior_vref_type == 3) V_ref_raw ~ lognormal(prior_vref_params[1], prior_vref_params[2]);

    if (prior_rho_type == 1) rho_raw ~ normal(prior_rho_params[1], prior_rho_params[2]);
    else if (prior_rho_type == 2) rho_raw ~ uniform(prior_rho_params[1], prior_rho_params[2]);
    else if (prior_rho_type == 3) rho_raw ~ lognormal(prior_rho_params[1], prior_rho_params[2]);
  }

  // =========================================================
  // 2. OPTIMIZED LIKELIHOOD (Vectorized Time Loop)
  // =========================================================
  // A. Pre-calculate Static Multipliers for every person
  vector[N] susc_eff;
  vector[N] inf_eff;
  vector[N] density_eff;

  for (n in 1:N) {
    real log_susc = (K_susc > 0) ? dot_product(X_susc[n], beta_susc) : 0.0;
    real log_inf  = (K_inf > 0)  ? dot_product(X_inf[n], beta_inf)  : 0.0;

    susc_eff[n] = phi_by_role[role_id[n]] * exp(log_susc);
    inf_eff[n]  = kappa_by_role[role_id[n]] * exp(log_inf);
    density_eff[n] = pow(1.0 / max(hh_size_people[hh_id[n]], 1), delta);
  }

  // B. Iterate TIME
  for (t in 1:T) {
    // --- Step 1: Calculate "Current Infectivity" for everyone ---
    vector[N] current_infectivity = rep_vector(0.0, N);
    for (n in 1:N) {
      if (Y[n, t] == 1) { // If infectious today
        real v_comp = 1.0;
        if (use_vl_data == 1) {
          v_comp = V_term_calc[n, t]; // Lookup pre-calculated value
        } else {
          // Empirical Curve Logic Fallback
          int inf_day = infection_day[n];
          if (inf_day != 0 && t >= inf_day) {
            int dt = t - inf_day + 1;
            if (use_curve_logic == 1 && dt <= T) v_comp = g_curve_est[dt];
          }
        }
        // Formula: Density * Kappa * (Beta1 + Beta2 * V)
        current_infectivity[n] = density_eff[n] * inf_eff[n] * (beta1 + beta2 * v_comp);
      }
    }

    // --- Step 2: Accumulate Force via Sparse Matrix ---
    vector[N] hh_force = rep_vector(0.0, N);
    for (k in 1:n_contacts) {
      hh_force[contact_tgt[k]] += contact_w[k] * current_infectivity[contact_src[k]];
    }

    // --- Step 3: Calculate Likelihood ---
    for (n in 1:N) {
      if (t >= start_risk[n]) {
        int t_stop = (infection_day[n] == 0) ? T : infection_day[n];
        if (t <= t_stop) {
          // Total Force = Susceptibility * (Community + Household)
          real lambda = susc_eff[n] * (alpha_comm * seasonal_forcing_mat[t, role_id[n]] + hh_force[n]);
          // Numerical stability clamp
          if (lambda > 20) lambda = 20;
          // Outcome: 1 if they got infected today, 0 otherwise
          int outcome = (t == infection_day[n]) ? 1 : 0;
          // Bernoulli Likelihood: Prob of infection = 1 - exp(-lambda)
          target += bernoulli_lpmf( outcome | 1.0 - exp(-lambda) );
        }
      }
    }
  }
}
