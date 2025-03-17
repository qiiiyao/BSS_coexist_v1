data {
  int<lower=0> N;   // number of observations
  int<lower=0> P;  // number of predicators
  int<lower=0> time_diff[N]; // year difference in our data
  int<lower=0> ftime;
  int<lower=0> ffield;
  matrix[N,P] X;     
  vector[N] G;     
  int<lower=1, upper=ftime> time[N];
  int<lower=1> field[N];
}


//transformed data {
  // Define transformed data variables here
 // vector[N] log_G; // Example transformation: Log-transform of observations

  // Compute the transformations
  //for (n in 1:N) {
   // log_G[n] = log(G[n]);
  //}
//}

parameters {
  
  real<lower=0> r; // intrinsic growth rate
  vector[P] a; // species interaction coefficients no permission
  
  vector[ftime] ran_t_raw;
  real<lower=0> sigma_ran_t;
  vector[ffield] ran_f_raw;
  real<lower=0> sigma_ran_f;
  real<lower=0> sigma_e; // sigma_e is observational error
}

transformed parameters {
  vector[ftime] ran_t;
  vector[ffield] ran_f; 

  //non-centred parameterization to improve efficiency and convergence of hierachical model
  // followed Yu_et_al_EL_2024
  ran_t = sigma_ran_t * ran_t_raw;
  ran_f = sigma_ran_f * ran_f_raw;
  
 // Rciker model 
  vector[N] yhat;
  
  for (i in 1:N){
   real exp_sum = 0;
   for (j in 1:P) {
    exp_sum += (X[i, j]^a[j]);
     } 
    yhat[i] = exp((r+X[i]*a+ran_t[time[i]]+ran_f[field[i]])*time_diff[i])/(1+exp_sum);
  }
}

model {
  // priors 
  r ~ normal(0, 1);
  for (i in 1:P){
      a[i] ~ normal(0, 1);
  }  
  sigma_ran_t ~ normal(0, 0.1);
  ran_t_raw ~ normal(0, 1);
  sigma_ran_f ~ normal(0, 0.1);
  ran_f_raw ~ normal(0, 1);
  sigma_e ~ normal(0, 1);
  
  // likelihood
  G ~ lognormal(yhat, sigma_e); // robust regression using student distribution following McElreath_2020
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = lognormal_lpdf(G[n] | yhat[n], sigma_e);
  }
}

