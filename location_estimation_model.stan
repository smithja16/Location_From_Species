
data {
  int<lower=1> N;                     // Total number of observations
  int<lower=1> S;                     // Number of species
  int<lower=1> M;                     // Number of months
  int<lower=1> N_known;               // Number of observations with known accurate locations
  int<lower=1> N_unknown;             // Number of observations with unknown true locations
  int<lower=1, upper=N> known_ids[N_known];     // Indices of observations with known locations
  int<lower=1, upper=N> unknown_ids[N_unknown]; // Indices of observations with unknown locations
  matrix[N, 2] loc_obs;               // Observed locations [x, y]
  int<lower=1, upper=M> month[N];     // Month of each observation
  int<lower=0> counts[N, S];          // Species counts
  real<lower=0> loc_sd;               // Prior SD for unknown locations
}

parameters {
  vector[S] alpha;                    // Species-specific intercepts
  vector[S] beta_x;                   // Species-specific x-coordinate effects
  vector[S] beta_y;                   // Species-specific y-coordinate effects
  matrix[S, M] gamma;                 // Species-specific monthly effects
  vector<lower=0>[S] phi;             // Negative binomial dispersion parameters
  matrix[N_unknown, 2] loc_true_unknown; // True locations for unknown observations
}

transformed parameters {
  matrix[N, 2] loc_true;              // True locations for all observations

  // Known locations are exactly as observed
  for (i in 1:N_known)
    loc_true[known_ids[i]] = loc_obs[known_ids[i]];

  // Unknown locations are parameters to estimate
  for (j in 1:N_unknown)
    loc_true[unknown_ids[j]] = loc_true_unknown[j];
}

model {
  // Priors
  alpha ~ normal(1, 1);
  beta_x ~ normal(0, 0.5); // Use normal(0, 0.5) or student_t(3, 0, 0.5)
  beta_y ~ normal(0, 0.5);
  to_vector(gamma) ~ normal(0, 0.5);
  phi ~ exponential(1);  // Prior for negative binomial dispersion parameter; gamma(2, 0.5) also possible

  // Prior for unknown locations - normal centered on observed locations
  for (j in 1:N_unknown)
    loc_true_unknown[j] ~ normal(loc_obs[unknown_ids[j]], loc_sd);

  // Likelihood
  for (i in 1:N)
    for (s in 1:S) {
      // Linear predictor
      real eta = alpha[s]
               + beta_x[s] * loc_true[i, 1]
               + beta_y[s] * loc_true[i, 2]
               + gamma[s, month[i]];

      // Negative binomial likelihood
      counts[i, s] ~ neg_binomial_2_log(eta, phi[s]);
    }
}

generated quantities {
  // Store log-likelihood for model assessment
  matrix[N, S] log_lik;

  for (i in 1:N) {
    for (s in 1:S) {
      real eta = alpha[s]
               + beta_x[s] * loc_true[i, 1]
               + beta_y[s] * loc_true[i, 2]
               + gamma[s, month[i]];

      log_lik[i, s] = neg_binomial_2_log_lpmf(counts[i, s] | eta, phi[s]);
    }
  }
}


