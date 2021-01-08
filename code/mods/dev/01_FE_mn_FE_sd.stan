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
  real<lower=0> sigma_mn_clny_to_mn; // should be species-specific...?
  vector[N_clny] err_mn_clny_to_mn;
  real<lower=0> sigma_sd_clny_to_mn;  // should be species-specific...?
  vector[N_clny] err_sd_clny_to_mn;
  
  vector[P_sd] alpha;  // intercept and slope for variance
  
  vector[P_mn] beta; // intercept and slope hyper-priors
}

transformed parameters {
  vector<lower=0>[N_clny] sigma_wkr_to_clny;
  vector[N_clny] mu;
  vector[N_clny] delta;
  vector[N_clny] clny_mn;

  for(i in 1:N_clny) {
    mu[i] = x_mn[i,] * beta;
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

  // likelihood
  y ~ normal(clny_mn[clny_id], sigma_wkr_to_clny[clny_id]);
}

generated quantities {
  vector[N_wkr] log_lik;
  vector[N_clny] delta_exp = exp(delta);
  
  for(i in 1:N_wkr) {
    log_lik[i] = normal_lpdf(y[i] | clny_mn[clny_id[i]], sigma_wkr_to_clny[clny_id[i]]);
  }
}
