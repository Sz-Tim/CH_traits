// based on 
// https://www.r-bloggers.com/bayesian-varying-effects-models-in-r-and-stan/

data {
  
  int<lower=0> N_clny; // number of colonies
  int<lower=0> N_wkr;  // number of workers
  int<lower=0> S; // number of participants
  int<lower=0> G;  // number of genera
  int<lower=0, upper=S> tax_i[S,2];  // species-genus lookup
  int<lower=0> P_mn; // number of predictors + intercept
  int<lower=0, upper=S> spp_id[N_clny]; // species id vector
  int<lower=0, upper=N_clny> clny_id[N_wkr]; // colony id vector
  matrix[N_clny,P_mn] x_mn; // matrix of predictors
  real y[N_wkr]; // y vector: observed worker traits
  
}



parameters {
  
  real<lower=1e-10, upper=1e2> sigma_clny_1_global;
  real<lower=1e-10, upper=1e2> sigma_clny_2_global;  // scale: log(d) ~ Norm(delta, sigma_clny_2)
  vector[N_clny] y_bar;  // within-colony mean
  vector<lower=0>[N_clny] d;  // within-colony sd
  
  real alpha;  // intercept and slope for variance
  
  // trait_mn slopes: beta, B, a
  // B ~ mvNorm(beta, sigma_B)
  // b ~ Norm(B, sigma_b)
  vector[P_mn] beta; // intercept and slope hyper-priors
  matrix[1,G] z_B; // species specific intercept and slope
  vector<lower=0>[1] sigma_B; // sd for intercept and slope
  cholesky_factor_corr[1] L_B; // cholesky correlation matrix
  matrix[1,S] z_b; // species specific intercept and slope
  vector<lower=0>[1] sigma_b; // sd for intercept and slope
  
}



transformed parameters {
  
  vector[N_clny] mu;  // predicted mean
  matrix[1,G] B; // non-centered version of z_B
  matrix[1,S] b; // non-centered version of z_b
  
  B = diag_pre_multiply(sigma_B, L_B) * z_B;
  b[1,] = B[1,tax_i[,2]] + z_b[1,] * sigma_b[1];
  
  for(j in 1:N_clny) {
    mu[j] = beta[1] + b[1,spp_id[j]] + x_mn[j,2:P_mn] * beta[2:P_mn];
  }
  
}



model {
  
  // priors
  y_bar ~ normal(mu, sigma_clny_1_global);
  sigma_clny_1_global ~ normal(5, 3);
  
  d ~ lognormal(alpha, sigma_clny_2_global);
  sigma_clny_2_global ~ normal(5, 3);
  
  alpha ~ normal(0, 1);
  
  beta ~ normal(0, 1);
  sigma_B ~ exponential(1);
  sigma_b ~ exponential(1);
  L_B ~ lkj_corr_cholesky(2);
  to_vector(z_B) ~ normal(0, 1);
  to_vector(z_b) ~ normal(0, 1);
  
  // likelihood
  y ~ normal(y_bar[clny_id], d[clny_id]);
  
}



generated quantities {
  
  vector[N_wkr] log_lik;
  vector[N_wkr] y_pred;
  vector[N_wkr] y_pred_;
  vector[N_clny] d_log = log(d);
  vector[N_clny] y_bar_pred;
  vector[N_clny] d_pred;
  //matrix[1,1] Omega_B;
  matrix[1,S] b_mn = beta[1] + b; 
  matrix[1,G] B_mn = beta[1] + B; 
  
  for(j in 1:N_clny) {
    y_bar_pred[j] = normal_rng(mu[j], sigma_clny_1_global);
    d_pred[j] = lognormal_rng(alpha, sigma_clny_2_global);
  }
  for(i in 1:N_wkr) {
    log_lik[i] = normal_lpdf(y[i] | y_bar[clny_id[i]], d[clny_id[i]]);
    y_pred[i] = normal_rng(y_bar[clny_id[i]], d[clny_id[i]]);
    y_pred_[i] = normal_rng(y_bar_pred[clny_id[i]], d_pred[clny_id[i]]);
  }
  
  //Omega_B = multiply_lower_tri_self_transpose(L_B);
  
}
