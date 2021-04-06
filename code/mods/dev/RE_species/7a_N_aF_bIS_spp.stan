// based on 
// https://www.r-bloggers.com/bayesian-varying-effects-models-in-r-and-stan/

data {
  
  int<lower=0> N_clny; // number of colonies
  int<lower=0> N_wkr;  // number of workers
  int<lower=0> S; // number of participants
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
  
  real<lower=1e-10, upper=1e2> sigma_clny_1_global;
  real<lower=1e-10, upper=1e2> sigma_clny_1[S];  // scale: y_bar ~ Norm(mu, sigma_clny_1)
  real<lower=1e-10, upper=1e2> sigma_clny_2_global;
  real<lower=1e-10, upper=1e2> sigma_clny_2[S];  // scale: log(d) ~ Norm(delta, sigma_clny_2)
  vector[N_clny] y_bar;  // within-colony mean
  vector<lower=0>[N_clny] d;  // within-colony sd
  
  // trait_sd slopes: alpha, A, a
  // A ~ mvNorm(alpha, sigma_A)
  // a ~ Norm(A, sigma_a)
  vector[P_sd] alpha;  // intercept and slope for trait_sd
  
  // trait_mn slopes: beta, B, a
  // B ~ mvNorm(beta, sigma_B)
  // b ~ Norm(B, sigma_b)
  vector[P_mn] beta; // intercept and slope hyper-priors
  cholesky_factor_corr[P_mn] L_b; // cholesky correlation matrix
  matrix[P_mn,S] z_b; // species specific intercept and slope, on N(0,1)
  vector<lower=0>[P_mn] sigma_b; // sd for intercept and slope
  
}



transformed parameters {
  
  vector[N_clny] mu;  // predicted mean
  vector[N_clny] delta;  // predicted sd
  matrix[P_mn,S] b; // non-centered version of z_b
  
  b = diag_pre_multiply(sigma_b, L_b) * z_b;

  for(j in 1:N_clny) {
    mu[j] = x_mn[j,] * (beta + b[,spp_id[j]]);
    delta[j] = x_sd[j,] * alpha;
  }
  
}



model {
  
  // model structure:
  // worker i, colony j, species k
  // y[i,j,k] ~ N(y_bar[j,k], d[j,k])
  // y_bar[j,k] ~ N(mu[j,k], sigma_clny_1[k])
  //    mu[j,k] = X[j] * b[k]
  //    sigma_clny_1[k] ~ Exp(sigma_clny_1_global)
  //    sigma_clny_1_global ~ N(5,3)[0,]
  // d[j,k] ~ logN(delta[j,k], sigma_clny_2[k])
  //    delta[j,k] = X[j] * alpha
  //    sigma_clny_2[k] ~ Exp(sigma_clny_2_global)
  //    sigma_clny_2_global ~ N(5,3)[0,]
  
  // priors
  y_bar ~ normal(mu, sigma_clny_1[spp_id]);
  to_vector(sigma_clny_1) ~ exponential(sigma_clny_1_global);
  sigma_clny_1_global ~ normal(5, 3);
  
  d ~ lognormal(delta, sigma_clny_2[spp_id]);
  to_vector(sigma_clny_2) ~ exponential(sigma_clny_2_global);
  sigma_clny_2_global ~ normal(5, 3);
  
  alpha ~ normal(0, 1);
  
  beta ~ normal(0, 1);
  sigma_b ~ normal(0, 1);
  L_b ~ lkj_corr_cholesky(1);
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
  matrix[P_mn,S] b_mn; 
  
  for(j in 1:N_clny) {
    y_bar_pred[j] = normal_rng(mu[j], sigma_clny_1[spp_id[j]]);
    d_pred[j] = lognormal_rng(delta[j], sigma_clny_2[spp_id[j]]);
  }
  for(i in 1:N_wkr) {
    log_lik[i] = normal_lpdf(y[i] | y_bar[clny_id[i]], d[clny_id[i]]);
    y_pred[i] = normal_rng(y_bar[clny_id[i]], d[clny_id[i]]);
    y_pred_[i] = normal_rng(y_bar_pred[clny_id[i]], d_pred[clny_id[i]]);
  }
  
  for(p in 1:P_mn) {
    b_mn[p,] = beta[p] + b[p,];
  }
  
}
