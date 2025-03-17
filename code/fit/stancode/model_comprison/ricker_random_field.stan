data{
  int<lower=0> n;              // length of time-series
  int<lower=0> nplot;           // number of mixes (i.e. communities)
  array[nplot] int Fieldnum;  //name of Field
  array[n,nplot] real g; // data: dim1=time, dim2=mixcode
  array[n,50,nplot] real N;
  array[n] int time;     // time-series
}


parameters {
  real<lower=0> r;
  vector[50] a;
  real<lower=0> sigma_e;
  real<lower=0> sigma_field;
  real<lower=0> sigma_rand_t;
  vector[n] rand_t_raw;
  real<lower=0> sigma_rand_f;
  vector[10] rand_f_raw;
}

transformed parameters {
 // Rciker model 
  vector[n] gsim;
  matrix[n,50] Nsim;

//non-centred parameterization to improve efficiency and convergence of hierachical model,
// following Yu_et_al_2024_Ecology Letters

  rand_t = rand_t_raw * sigma_rand_t;
  rand_f = rand_f_raw * sigma_rand_f;
  
 // Ricker competition model
  for (x in 1:nplot){ // loop for mixes
    Nsim = to_matrix(N[,,x]);
    // fitting
      for (t in 1:n){
          gsim[t] = exp(r-Nsim[t]*a+ rand_t[t] + rand_f[Fieldnum[x]]);
        } // t
      } // x
}


model {
  
  // priors
  r ~ normal(0,1);
  for (i in 1:50){
      a[i] ~ normal(0,1);
  }
  // random effect structure
  sigma_rand_t ~ normal(0, 0.1);
  sigma_rand_f ~ normal(0, 0.1);
  #sigma_e ~ gamma(0.001, 0.001);
  rand_t_raw ~ std_normal();
  rand_f_raw ~ std_normal();
  
  // likelihood
  for (x in 1:nplot){ // loop for mixes
    // fitting
      for (t in 1:n){
          log(g[t,x]) ~ normal(log(gsim[t]), sigma_e);
        } // t
      } // x
  } 
  

generated quantities {
  matrix[n, nplot] log_lik;
  for (x in 1:nplot){ // loop for mixes
    // fitting
      for (t in 1:n){
    log_lik[t,x] = normal_lpdf(log(g[t,x]) | log(gsim[t]), sigma_e);
  }
 }
}

  
