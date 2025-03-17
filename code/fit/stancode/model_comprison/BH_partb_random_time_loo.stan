data {
  int<lower=0> N;   // number of observations
  int<lower=0> P;  // number of predicators
  int<lower=0> ftime;
  vector[N] time_diff;
  matrix[N,P] X;   
  vector[N] G;     
  int<lower=1, upper=ftime> time[N];
}
 
parameters {
  real <lower=0> b;
  real<lower=0> r;
  vector<lower=0>[P] a;
  vector[ftime] ran_t;
  real<lower=0> sigma_t;
  real<lower=0> sigma_e;
}
transformed parameters {
 // Rciker model 
  vector[N] yhat;
  for (i in 1:N)
    yhat[i] = exp((r+ran_t[time[i]])*time_diff[i])/(1+((X[i]*a)^b));
  
}
model {
  // priors 
  r ~ normal(0, 1);
  for (i in 1:P){
      a[i] ~ normal(0, 1);
  }  
  sigma_t ~ normal(0, 1);
  ran_t ~ normal(0, sigma_t);
  sigma_e ~ normal(0, 1);
  // likelihood
  log(G) ~ student_t(2.5,log(yhat), sigma_e);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) { 
    log_lik[n] = student_t_lpdf(log(G[n]) |2.5, log(yhat[n]), sigma_e);
  }
}

