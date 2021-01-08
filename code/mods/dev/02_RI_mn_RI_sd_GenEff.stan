// based on 
// https://www.r-bloggers.com/bayesian-varying-effects-models-in-r-and-stan/

data {
  
  int<lower=0> N_clny; // number of colonies
  int<lower=0> N_wkr;  // number of workers
  int<lower=0> S; // number of species
  int<lower=0> G;  // number of genera
  int<lower=0> tax_i[S,2];  // species-genus lookup
  int<lower=0> P_mn; // number of predictors + intercept for trait_mean
  int<lower=0> P_sd; // number of predictors + intercept for trait_sd
  int spp_id[N_clny]; // species id vector
  int clny_id[N_wkr]; // colony id vector
  matrix[N_clny,P_mn] x_mn; // matrix of predictors for trait_mean
  matrix[N_clny,P_sd] x_sd; // matrix of predictors for trait_sd
  real y[N_wkr]; // y vector: observed worker values
  
}



parameters {
  
  // error terms
  vector[N_clny] err_mn_clny_to_pred;  // error for modeled_clny_mn ~ predicted_mn
  vector[N_clny] err_sd_clny_to_pred;  // error for modeled_clny_sd ~ predicted_sd
  real<lower=0> sigma_mn_clny_to_pred;  // sd for distribution of err_mn_clny_to_pred
  real<lower=0> sigma_sd_clny_to_pred;  // sd for distribution of err_sd_clny_to_pred
  
  // trait_sd slopes: alpha, A, a
  // A ~ mvNorm(alpha, sigma_A)
  // a ~ Norm(A, sigma_a)
  vector[P_sd] alpha;  // intercept and slope for trait_sd
  matrix[1,G] z_A; // species specific intercept and slope
  vector<lower=0>[1] sigma_A; // sd for intercept and slope
  cholesky_factor_corr[1] L_A; // cholesky correlation matrix
  matrix[1,S] z_a; // species specific intercept and slope
  vector<lower=0>[1] sigma_a; // sd for intercept and slope
  
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
  
  matrix[1,G] A; // non-centered version of z_A
  matrix[1,S] a; // non-centered version of z_a
  matrix[1,G] B; // non-centered version of z_B
  matrix[1,S] b; // non-centered version of z_b
  vector<lower=0>[N_clny] sigma_wkr_to_clny;
  vector[N_clny] mu;
  vector[N_clny] delta;
  vector[N_clny] clny_mn;
  
  A = diag_pre_multiply(sigma_A, L_A) * z_A;
  a[1,] = A[1,tax_i[,2]] + z_a[1,] * sigma_a[1];
  B = diag_pre_multiply(sigma_B, L_B) * z_B;
  b[1,] = B[1,tax_i[,2]] + z_b[1,] * sigma_b[1];

  for(i in 1:N_clny) {
    mu[i] = beta[1] + b[1,spp_id[i]] + x_mn[i,2:P_mn] * beta[2:P_mn];
    delta[i] = alpha[1] + a[1,spp_id[i]] + x_sd[i,2:P_sd] * alpha[2:P_sd];
    clny_mn[i] = mu[i] + err_mn_clny_to_pred[i];
    sigma_wkr_to_clny[i] = exp(delta[i] + err_sd_clny_to_pred[i]);
  }
  
}



model {
  
  // priors
  err_mn_clny_to_pred ~ normal(0, sigma_mn_clny_to_pred);
  sigma_mn_clny_to_pred ~ exponential(1);
  
  err_sd_clny_to_pred ~ normal(0, sigma_sd_clny_to_pred);
  sigma_sd_clny_to_pred ~ exponential(5);
  
  alpha ~ normal(0, 1);
  sigma_A ~ exponential(1);
  sigma_a ~ exponential(1);
  L_A ~ lkj_corr_cholesky(2);
  to_vector(z_A) ~ normal(0, 1);
  to_vector(z_a) ~ normal(0, 1);
  
  beta ~ normal(0, 1);
  sigma_B ~ exponential(1);
  sigma_b ~ exponential(1);
  L_B ~ lkj_corr_cholesky(2);
  to_vector(z_B) ~ normal(0, 1);
  to_vector(z_b) ~ normal(0, 1);

  // likelihood
  y ~ normal(clny_mn[clny_id], sigma_wkr_to_clny[clny_id]);
  
}



generated quantities {
  
  vector[N_wkr] log_lik;
  matrix[1,1] Omega_A;
  matrix[1,1] Omega_B;
  matrix[1,S] a_mn = alpha[1] + a; 
  matrix[1,S] b_mn = beta[1] + b; 
  vector[N_clny] delta_exp = exp(delta);
  vector[N_wkr] y_pred;
  
  for(i in 1:N_wkr) {
    log_lik[i] = normal_lpdf(y[i] | clny_mn[clny_id[i]], sigma_wkr_to_clny[clny_id[i]]);
    y_pred[i] = normal_rng(clny_mn[clny_id[i]], sigma_wkr_to_clny[clny_id[i]]);
  }
  Omega_A = multiply_lower_tri_self_transpose(L_A);
  Omega_B = multiply_lower_tri_self_transpose(L_B);
  
}
