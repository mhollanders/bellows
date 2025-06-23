functions {
  /**
   * Hawkes process with exponential intensity log density. See 
   * Daniel O'Shea. 2020. Modelling with Hawkes Processes, eqs. 4.26--4.27
   * 
   * @param y      Waiting y of detections
   * @param mu     Baseline detection rate
   * @param hp     Hawkes process parameters alpha and beta
   * @param Delta  Interval length
   * @param tau    Time remaining from each detection to end of interval
   * 
   * @return       Log density of waiting y
   */
  real hawkes_exp_lpdf(vector y, real mu, vector hp, real Delta, vector tau) {
    int D = size(y) - 1;
    vector[D] exp_neg_beta_y = exp(-hp[2] * y[1:D]);
    vector[D] A = zeros_vector(D);
    for (d in 2:D) {
      A[d] = exp_neg_beta_y[d] * (1 + A[d - 1]);
    }
    return sum(log(mu + hp[1] * A))
           - mu * Delta
           + hp[1] / hp[2] * sum(expm1(-hp[2] * tau));
  }
}

data {
  array[2] int<lower=0> N;  // number of drone counts and ARU observations
  int<lower=0> I, J, X, B_max;  // number of observations, sites, and covariates, and maximum number of bellows
  array[N[2]] int<lower=1, upper=I> site;  // site indicator
  array[N[2]] int<lower=1, upper=J> survey;  // survey indicator
  array[N[1]] int<lower=0> C;  // drone counts
  array[N[2]] int<lower=0> B;  // number of bellows
  row_vector<lower=0>[N[2]] Delta;  // drone survey lengths (hours)
  array[N[2]] vector<lower=0>[B_max + 1] y;  // waiting times between bellows
  matrix[X, N[2]] x;  // covariates
  array[2] int<lower=0, upper=1> NB;  // indicators for negative binomial
  int<lower=0, upper=1> HP;  // indicator for Hawkes process
}

transformed data {
  array[N[2]] int Bp_1;  // number of bellows + 1
  array[N[2]] vector[B_max + 1] tau;
  for (n in 1:N[2]) {
    Bp_1[n] = B[n] + 1;
    tau[n, :Bp_1[n]] = Delta[n] - cumulative_sum(y[n, :Bp_1[n]]);
  }
  int idx = NB[1] ? 2 : 1;  // index for negative binomial overdispersion for bellows
}

parameters {
  vector<lower=0>[2] gamma0;  // intercepts (original scale)
  vector<lower=0, upper=1>[2] R2;  // R-squared
  row_stochastic_matrix[2, X + 2] zeta;  // variance partitions
  array[2] cholesky_factor_corr[2] Omega_L;  // Cholesky factors of correlation matrices
  matrix[2, X + I + J] z;  // standard normal variates for coefficients and random effects
  vector<lower=0>[sum(NB)] phi;  // negative binomial overdispersion
  row_vector<lower=0>[NB[2] * N[2]] u;  // negative binomial variates
  positive_ordered[HP * 2] hp;  // Hawkes process parameters
}

transformed parameters {
  // R2D2 pseudo-variance and scales
  vector[2] tau2 = R2 ./ (1 - R2), pseudo_var = inv(gamma0);
  pseudo_var[1] += NB[1] ? inv(phi[1]) : 0;
  pseudo_var[2] += NB[2] ? inv(phi[idx]) : 0;
  matrix[2, X + 2] scales = sqrt(diag_pre_multiply(tau2 .* pseudo_var, zeta));
  
  // coefficients and random site and survey effects
  matrix[2, X] gamma = scales[:, 1:X] .* z[:, :X];
  matrix[2, I] theta = diag_pre_multiply(scales[:, X + 1], Omega_L[1])
                       * z[:, (X + 1):(X + I)];
  matrix[2, J] epsilon = diag_pre_multiply(scales[:, X + 2], Omega_L[2])
                         * z[:, X + I + 1:X + I + J];
  
  // log-linear expectations including negative binomial variates
  matrix[2, N[2]] mu;
  for (d in 1:2) {
    mu[d, :N[d]] = exp(log(gamma0[d]) 
                       + gamma[d] * x[:, :N[d]]
                       + theta[d, site[:N[d]]]
                       + epsilon[d, survey[:N[d]]]);
  }
  if (NB[2]) {
    mu[2] .*= u;
  }
  
  // priors (gamma and Omega_L are implicit)
  real lprior = exponential_lpdf(to_vector(gamma0) | 0.1)
                + beta_lpdf(R2 | 1, 1)
                + std_normal_lpdf(to_vector(z));
  if (sum(NB)) {
    lprior += gamma_lpdf(phi | 2, 0.1);
    if (NB[2]) {
      lprior += gamma_lpdf(u | phi[idx], phi[idx]);
    }
  }
  if (HP) {
    lprior += exponential_lpdf(hp | 0.1);
  }
}

model {
  target += lprior;
  target += NB[1] ?
            neg_binomial_2_lupmf(C | mu[1, 1:N[1]], phi[1])
            : poisson_lupmf(C | mu[1, 1:N[1]]);
  for (n in 1:N[2]) {
    if (B[n]) {
      if (HP) {
        target += hawkes_exp_lpdf(y[n, :Bp_1[n]] | mu[2, n], hp, Delta[n], 
                                  tau[n, :Bp_1[n]]);
      } else {
        target += exponential_lupdf(y[n, :B[n]] | mu[2, n])
                  + exponential_lccdf(y[n, Bp_1[n]] | mu[2, n]);
      }
    } else {
      target += exponential_lccdf(Delta[n] | mu[2, n]);
    }
  }
}

generated quantities {
  // correlations
  vector[2] rho;
  for (i in 1:2) {
    rho[i] = multiply_lower_tri_self_transpose(Omega_L[i])[1, 2];
  }
  
  // log likelihoods
  vector[N[1]] log_lik1;
  vector[N[2]] log_lik2, log_lik;
  for (n in 1:N[1]) {
    log_lik1[n] = NB[1] ?
                  neg_binomial_2_lpmf(C[n] | mu[1, n], phi[1])
                  : poisson_lpmf(C[n] | mu[1, n]);
  }
  for (n in 1:N[2]) {
    if (B[n]) {
      if (HP) {
        log_lik2[n] = hawkes_exp_lpdf(y[n, :Bp_1[n]] | mu[2, n], hp, Delta[n], 
                                      tau[n, :Bp_1[n]]);
      } else {
        log_lik2[n] = exponential_lpdf(y[n, :B[n]] | mu[2, n])
                      + exponential_lccdf(y[n, Bp_1[n]] | mu[2, n]);
      }
    } else {
      log_lik2[n] = exponential_lccdf(Delta[n] | mu[2, n]);
    }
  }
  log_lik = log_lik2;
  log_lik[:N[1]] += log_lik1;
  
  // posterior predictions
  array[N[1]] int C_rep;
  array[N[2]] int B_rep;
  {
    C_rep = NB[1] ?
            neg_binomial_2_rng(mu[1, :N[1]], phi[1])
            : poisson_rng(mu[1, :N[1]]);
    row_vector[N[2]] mu_Delta = mu[2] .* Delta 
                               ./ (NB[2] ? u : ones_row_vector(N[2]))
                               ./ (HP ? 1 - hp[1] / hp[2] : 1);
    B_rep = NB[2] ?
            neg_binomial_2_rng(mu_Delta, phi[idx])
            : poisson_rng(mu_Delta);
  }
}
