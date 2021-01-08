// based on 
// https://www.r-bloggers.com/bayesian-varying-effects-models-in-r-and-stan/

data {
  
  int<lower=0> N_clny; // number of colonies
  int<lower=0> N_wkr;  // number of workers
  int<lower=0> S; // number of participants
  int<lower=0> G;  // number of genera
  int<lower=0, upper=S> tax_i[S,2];  // species-genus lookup
  int<lower=0> P_mn; // number of predictors + intercept
  int<lower=0> P_sd; // number of predictors + intercept
  int<lower=0, upper=S> spp_id[N_clny]; // species id vector
  int<lower=0, upper=N_clny> clny_id[N_wkr]; // colony id vector
  matrix[N_clny,P_mn] x_mn; // matrix of predictors
  matrix[N_clny,P_sd] x_sd; // matrix of predictors
  real y[N_wkr]; // y vector: observed worker traits
  
}



parameters {
  
  real<lower=0> sigma_clny_1_global;
  real<lower=0> sigma_clny_1[S];  // scale: y_bar ~ Norm(mu, sigma_clny_1)
  real<lower=0> sigma_clny_2_global;
  real<lower=0> sigma_clny_2[S];  // scale: log(d) ~ Norm(delta, sigma_clny_2)
  vector[N_clny] y_bar;  // within-colony mean
  vector<lower=0>[N_clny] d;  // within-colony sd
  
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
  matrix[P_mn,G] z_B; // species specific intercept and slope
  vector<lower=0>[P_mn] sigma_B; // sd for intercept and slope
  cholesky_factor_corr[P_mn] L_B; // cholesky correlation matrix
  matrix[P_mn,S] z_b; // species specific intercept and slope
  vector<lower=0>[P_mn] sigma_b; // sd for intercept and slope
  
}



transformed parameters {
  
  vector[N_clny] mu;  // predicted mean
  vector[N_clny] delta;  // predicted sd
  matrix[1,G] A; // non-centered version of z_A
  matrix[1,S] a; // non-centered version of z_a
  matrix[P_mn,G] B; // non-centered version of z_B
  matrix[P_mn,S] b; // non-centered version of z_b
  
  A = diag_pre_multiply(sigma_A, L_A) * z_A;
  a[1,] = A[1,tax_i[,2]] + z_a[1,] * sigma_a[1];
  B = diag_pre_multiply(sigma_B, L_B) * z_B;
  for(p in 1:P_mn) {
    b[p,] = B[p,tax_i[,2]] + z_b[p,] * sigma_b[p];
  }

  for(j in 1:N_clny) {
    mu[j] = x_mn[j,] * (beta + b[,spp_id[j]]);
    delta[j] = alpha[1] + a[1,spp_id[j]] + x_sd[j,2:P_sd] * alpha[2:P_sd];
  }
  
}



model {
  
  // priors
  y_bar ~ normal(mu, sigma_clny_1[spp_id]);
  to_vector(sigma_clny_1) ~ exponential(sigma_clny_1_global);
  sigma_clny_1_global ~ normal(5, 3);
  
  d ~ lognormal(delta, sigma_clny_2[spp_id]);
  to_vector(sigma_clny_2) ~ exponential(sigma_clny_2_global);
  sigma_clny_2_global ~ normal(5, 3);
  
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
  y ~ normal(y_bar[clny_id], d[clny_id]);
  
}



generated quantities {
  
  vector[N_wkr] log_lik;
  vector[N_wkr] y_pred;
  vector[N_wkr] y_pred_;
  vector[N_clny] delta_exp = exp(delta);
  vector[N_clny] d_log = log(d);
  vector[N_clny] y_bar_pred;
  vector[N_clny] d_pred;
  //matrix[1,1] Omega_A;
  matrix[1,S] a_mn = alpha[1] + a; 
  //matrix[P_mn,P_mn] Omega_B;
  matrix[P_mn,S] b_mn; 
  matrix[P_mn,G] B_mn; 
  
  for(j in 1:N_clny) {
    y_bar_pred[j] = normal_rng(mu[j], sigma_clny_1[spp_id[j]]);
    d_pred[j] = lognormal_rng(delta[j], sigma_clny_2[spp_id[j]]);
  }
  for(i in 1:N_wkr) {
    log_lik[i] = normal_lpdf(y[i] | y_bar[clny_id[i]], d[clny_id[i]]);
    y_pred[i] = normal_rng(y_bar[clny_id[i]], d[clny_id[i]]);
    y_pred_[i] = normal_rng(y_bar_pred[clny_id[i]], d_pred[clny_id[i]]);
  }
  
  //Omega_A = multiply_lower_tri_self_transpose(L_A);
  //Omega_B = multiply_lower_tri_self_transpose(L_B);
  for(p in 1:P_mn) {
    b_mn[p,] = beta[p] + b[p,];
    B_mn[p,] = beta[p] + B[p,];
  }
  
}
