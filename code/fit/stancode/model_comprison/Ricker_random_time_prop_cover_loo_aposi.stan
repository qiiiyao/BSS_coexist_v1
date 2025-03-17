data {
  int<lower=0> N;   // number of observations
  int<lower=0> P;  // number of predicators
  int<lower=0> ftime;
  vector[N] time_diff;
  matrix[N,P-1] X_t_inter; 
  vector[N] X_t_intra;
  vector<lower=0, upper=1>[N] Y_t1; 
  int<lower=1, upper=ftime> time[N];
}

parameters {
  real<lower=0> r;
  real<upper=0> a_intra; // species interaction coefficients no permission
  vector[P-1] a_inter;
  vector[ftime] ran_t;
  real<lower=0> sigma_t;
  real<lower=0> sigma_e;
}

transformed parameters {
 // Rciker model 
  vector[N] yhat;
  
  for (i in 1:N)
    yhat[i] = X_t_intra[i]*exp((r+X_t_intra[i]*a_intra+X_t_inter[i]*a_inter+ran_t[time[i]])*time_diff[i]);
  
}
model {
  // priors 
  r ~ normal(0, 1);
  for (i in 1:(P-1)){
      a_inter[i] ~ normal(0, 1);
  }  
  a_intra ~ normal(-1, 1);
  sigma_t ~ normal(0, 1);
  ran_t ~ normal(0, sigma_t);
  sigma_e ~ normal(0, 1);
  // likelihood
  Y_t1 ~ lognormal(yhat, sigma_e); 
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) { 
    log_lik[n] = lognormal_lpdf(Y_t1[n] | yhat[n], sigma_e);
  }
}

