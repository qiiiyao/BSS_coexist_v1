data {
  int<lower=0> N;   // number of observations
  int<lower=0> P;  // number of predicators
  int<lower=0> ftime;
  vector[N] time_diff;
  matrix[N,P] X;   
  vector[N] G;     
  int<lower=1, upper=ftime> time[N];
}

//transformed data {
  // Define transformed data variables here
  //vector[N] log_G; // Example transformation: Log-transform of observations

  // Compute the transformations
  //for (n in 1:N) {
    //log_G[n] = log(G[n])/time_diff[n];
  //}
//}

parameters {
  real b;
  real r;
  vector[P] a;
  vector[ftime] ran_t;
  real<lower=0> sigma_t;
  real<lower=0> sigma_e;
}
transformed parameters {
 // Rciker model 
  vector[N] yhat;
  for (i in 1:N)
    yhat[i] = exp((r+ran_t[time[i]]) * time_diff[i])/((1+X[i]*a)^b);
  
}

model {
  // priors 
  r ~ normal(0, 1);
  b ~ normal(0, 1);
  for (i in 1:P){
      a[i] ~ normal(0, 1);
  }  
  sigma_t ~ normal(0, 1);
  ran_t ~ normal(0, sigma_t);
  sigma_e ~ normal(0, 1);
  // likelihood
  G ~ lognormal(yhat, sigma_e);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) { 
    log_lik[n] = lognormal_lpdf(G[n] | yhat[n], sigma_e);
  }
}

