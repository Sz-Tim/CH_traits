// based on 
// https://www.r-bloggers.com/bayesian-varying-effects-models-in-r-and-stan/

data {
  
  int<lower=0> N_clny; // number of colonies
  int<lower=0> N_wkr;  // number of workers
  int<lower=0> S; // number of participants
  int<lower=0> P_mn; // number of predictors + intercept
  int<lower=0> P_sd; // number of predictors + intercept
  int<lower=0, upper=S> spp_id[N_clny]; // species id vector
  int<lower=0, upper=N_clny> clny_id[N_wkr]; // colony id vector
  matrix[N_clny,P_mn] x_mn; // matrix of predictors
  matrix[N_clny,P_sd] x_sd; // matrix of predictors
  real y[N_wkr]; // y vector: observed worker traits
  
}



parameters {
  
  real<lower=0> sigma_clny_1;  // scale: y_bar ~ Norm(mu, sigma_clny_1)
  real<lower=0> sigma_clny_2;  // scale: log(d) ~ Norm(delta, sigma_clny_2)
  vector[N_clny] y_bar;  // within-colony mean
  vector<lower=0>[N_clny] d;  // within-colony sd
  
  vector[P_sd] alpha;  // intercept and slope for variance
  
  vector[P_mn] beta; // intercept and slope hyper-priors
  
  real<lower=1> nu_w;
  
}



transformed parameters {
  
  vector[N_clny] mu;  // predicted mean
  vector[N_clny] delta;  // predicted sd

  for(j in 1:N_clny) {
    mu[j] = x_mn[j,] * beta;
    delta[j] = x_sd[j,] * alpha;
  }
  
}



model {
  
  // priors
  y_bar ~ normal(mu, sigma_clny_1);
  sigma_clny_1 ~ exponential(1);
  
  d ~ lognormal(delta, sigma_clny_2);
  sigma_clny_2 ~ exponential(5);
  
  alpha ~ normal(0, 1);
  
  beta ~ normal(0, 1);
  
  nu_w ~ gamma(2, 0.5);

  // likelihood
  y ~ student_t(nu_w, y_bar[clny_id], d[clny_id]);
  
}



generated quantities {
  
  vector[N_wkr] log_lik;
  vector[N_wkr] y_pred;
  vector[N_wkr] y_pred_;
  vector[N_clny] delta_exp = exp(delta);
  vector[N_clny] d_log = log(d);
  vector[N_clny] y_bar_pred;
  vector[N_clny] d_pred;
  
  for(j in 1:N_clny) {
    y_bar_pred[j] = normal_rng(mu[j], sigma_clny_1);
    d_pred[j] = lognormal_rng(delta[j], sigma_clny_2);
  }
  for(i in 1:N_wkr) {
    log_lik[i] = student_t_lpdf(y[i] | nu_w, y_bar[clny_id[i]], d[clny_id[i]]);
    y_pred[i] = normal_rng(y_bar[clny_id[i]], d[clny_id[i]]);
    y_pred_[i] = normal_rng(y_bar_pred[clny_id[i]], d_pred[clny_id[i]]);
  }
  
}
