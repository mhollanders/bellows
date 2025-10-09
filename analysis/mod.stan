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
  int<lower=0> I, J, R, X, B_max;  // number of sites, surveys, regions, and covariates, and maximum number of bellows
  array[N[2]] int<lower=1, upper=I> site;  // site indicator
  array[I] int<lower=1, upper=R> region_i;  // region indicator for sites
  array[J] int<lower=1, upper=R> region_j;  // region indicator for survey
  array[N[2]] int<lower=1, upper=J> survey;  // survey indicator
  array[N[1]] int<lower=0> C;  // drone counts
  array[N[2]] int<lower=0> B;  // number of bellows
  row_vector<lower=0>[N[1]] area;  // areas surveyed with drones
  row_vector<lower=0>[N[2]] Delta;  // ARU survey lengths (hours)
  array[N[2]] vector<lower=0>[B_max + 1] y;  // waiting times between bellows
  matrix[X, N[2]] x;  // covariates
  array[2] int<lower=0, upper=1> NB;  // indicators for negative binomial
  int<lower=0, upper=1> HP;  // indicator for Hawkes process
  int<lower=0> pred;  // number of predictions
}

transformed data {
  array[N[2]] int Bp_1;  // number of bellows + 1
  array[N[2]] vector[B_max + 1] tau;
  for (n in 1:N[2]) {
    Bp_1[n] = B[n] + 1;
    tau[n, :Bp_1[n]] = Delta[n] - cumulative_sum(y[n, :Bp_1[n]]);
  }
  int idx = NB[1] ? 2 : 1;  // index for negative binomial overdispersion for bellows
  real r_scale = inv(sqrt(1 - inv(R)));
}

parameters {
  vector<lower=0>[2] gamma0;  // intercepts (original scale) and scales
  array[2] sum_to_zero_vector[R] gamma0_r;  // region-level offsets
  vector<lower=0, upper=1>[2] R2;  // R-squared
  row_stochastic_matrix[2, X + 3] zeta;  // variance partitions
  vector<lower=-1, upper=1>[2] rho_a;  // site and survey correlations
  vector<lower=0>[2] rho_t;  // scales
  array[2] sum_to_zero_vector[R] rho_r;  // region-level offsets
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
  matrix[2, X + 3] scales = sqrt(diag_pre_multiply(tau2 .* pseudo_var, zeta));
  
  // coefficients and random region, site, and survey effects
  matrix[2, X] gamma = scales[:, :X] .* z[:, :X];
  matrix[R, 2] rho;
  matrix[2, I] theta = rep_matrix(log(gamma0), I);
  matrix[2, J] epsilon;
  {
    matrix[2, 2] O = diag_matrix(ones_vector(2));
    array[2, R] matrix[2, 2] S_L;
    for (d in 1:2) {
      rho[:, d] = tanh(atanh(rho_a[d]) + rho_t[d] * rho_r[d]);
      for (r in 1:R) {
        O[1, 2] = rho[r, d];
        O[2, 1] = rho[r, d];
        S_L[d, r] = diag_pre_multiply(scales[:, X + d], cholesky_decompose(O));
      }
      theta[d] += scales[d, X + 3] * gamma0_r[d, region_i]';
    }
    for (i in 1:I) {
      theta[:, i] += S_L[1, region_i[i]] * z[:, X + i];
    }
    for (j in 1:J) {
      epsilon[:, j] = S_L[2, region_j[j]] * z[:, X + I + j];
    }
  }
  
  // log-linear expectations including negative binomial variates
  matrix[2, N[2]] mu;
  for (d in 1:2) {
    mu[d, :N[d]] = exp(theta[d, site[:N[d]]]
                       + gamma[d] * x[:, :N[d]]
                       + epsilon[d, survey[:N[d]]]);
  }
  mu[1, :N[1]] .*= area;
  if (NB[2]) {
    mu[2] .*= u;
  }
  
  // priors
  real lprior = exponential_lpdf(gamma0 | 0.1)
                + beta_lpdf(R2 | 1, 1)
                + std_normal_lpdf(to_vector(z))
                + exponential_lpdf(rho_t | 2);
  for (d in 1:2) {
    lprior += normal_lpdf(gamma0_r[d] | 0, r_scale)
              + normal_lpdf(rho_r[d] | 0, r_scale);
  }
  if (sum(NB)) {
    lprior += inv_gamma_lpdf(phi | 0.4, 0.3);
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
  
  // // posterior site-level predictions
  // array[2] matrix[pred, 2] theta_pred;
  // for (r in 1:2) {
  //   real log_C0 = log(gamma0[1]) + scales[1, X + 3] .* gamma0_r[1, r],
  //        log_B0 = log(gamma0[2]) + scales[2, X + 3] .* gamma0_r[2, r];
  //   theta_pred[r] = append_col(linspaced_vector(pred, 0, 3),
  //                              rep_vector(log_B0, pred));
  //   theta_pred[r, :, 2] += rho[1, r] * scales[2, X + 1] / scales[1, X + 1]
  //                          * (theta_pred[r, :, 1] - log_C0);
  // }
  
  // posterior site-level predictions
  matrix[pred, 2] theta_pred = append_col(linspaced_vector(pred, 0, 3),
                                          rep_vector(log(gamma0[2]), pred));
  theta_pred[:, 2] += rho_a[1] * scales[2, X + 1] / scales[1, X + 1]
                      * (theta_pred[:, 1] - log(gamma0[1]));
  array[pred] real theta_rep = normal_rng(theta_pred[:, 2],
                                          scales[2, X + 1]
                                          * sqrt(1 - square(rho_a[1])));
}
