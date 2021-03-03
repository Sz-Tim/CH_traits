// based on 
// https://www.r-bloggers.com/bayesian-varying-effects-models-in-r-and-stan/

data {
  int N_clny; // number of colonies
  int N_wkr;  // number of workers
  int S; // number of participants
  int P_mn; // number of predictors + intercept
  int P_sd; // number of predictors + intercept
  int spp_id[N_clny]; // participant id vector
  int clny_id[N_wkr]; // participant id vector
  matrix[N_clny,P_mn] x_mn; // matrix of predictors
  matrix[N_clny,P_sd] x_sd; // matrix of predictors
  real y[N_wkr]; // y vector: observed colony means
}

parameters {
  real<lower=1e-10, upper=1e2> sigma_mn_clny_to_mn;
  vector[N_clny] err_mn_clny_to_mn;
  real<lower=1e-10, upper=1e2> sigma_sd_clny_to_mn;
  vector[N_clny] err_sd_clny_to_mn;
  
  vector[P_sd] alpha;  // intercept and slope for variance
  
  vector[P_mn] beta; // intercept and slope hyper-priors
  matrix[P_mn,S] z_b; // species specific intercept and slope
  vector<lower=0>[P_mn] sigma_b; // sd for intercept and slope
  cholesky_factor_corr[P_mn] L_b; // cholesky correlation matrix
}

transformed parameters {
  matrix[P_mn, S] b; // non-centered version of z_b
  vector<lower=0>[N_clny] sigma_wkr_to_clny;
  vector[N_clny] mu;
  vector[N_clny] delta;
  vector[N_clny] clny_mn;
  
  b = diag_pre_multiply(sigma_b, L_b) * z_b;

  for(i in 1:N_clny) {
    mu[i] = x_mn[i,] * (beta + b[,spp_id[i]]);
    delta[i] = x_sd[i,] * alpha;
    clny_mn[i] = mu[i] + err_mn_clny_to_mn[i];
    sigma_wkr_to_clny[i] = exp(delta[i] + err_sd_clny_to_mn[i]);
  }
}

model {
  // priors
  err_mn_clny_to_mn ~ normal(0, sigma_mn_clny_to_mn);
  sigma_mn_clny_to_mn ~ exponential(1);
  
  err_sd_clny_to_mn ~ normal(0, sigma_sd_clny_to_mn);
  sigma_sd_clny_to_mn ~ exponential(5);
  
  alpha ~ normal(0, 1);
  
  beta ~ normal(0, 1);
  sigma_b ~ exponential(1);
  L_b ~ lkj_corr_cholesky(2);
  to_vector(z_b) ~ normal(0, 1);

  // likelihood
  y ~ normal(clny_mn[clny_id], sigma_wkr_to_clny[clny_id]);
}

generated quantities {
  vector[N_wkr] log_lik;
  matrix[P_mn, P_mn] Omega_b;
  matrix[P_mn,S] b_mn; // non-centered version of z_b
  vector[N_clny] delta_exp = exp(delta);
  
  for(p in 1:P_mn) {
    b_mn[p,] = beta[p] + b[p,];
  }
  
  for(i in 1:N_wkr) {
    log_lik[i] = normal_lpdf(y[i] | clny_mn[clny_id[i]], sigma_wkr_to_clny[clny_id[i]]);
  }
  Omega_b = multiply_lower_tri_self_transpose(L_b);
}
