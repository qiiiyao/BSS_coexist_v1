
data {
  int<lower=0> s; // number of plots
  int<lower=0> timestep; //number of timesteps = ?
  int Y_O[s,timestep]; //species occurence data
  real Y_C[s,timestep]; // species cover data
}


parameters {
  vector[timestep] A_t;
  vector[s] C_f;
  real b1; //bio interaction effects on occurence
  real sigmaplot_o;
  real psi_o;
 
}

transformed parameters{
  real<lower=0> tauplot_o;
  tauplot_o = sqrt(1/(sigmaplot_o*sigmaplot_o));
  array[s,timestep-1] real<lower=0,upper=1> mu;
  
   for(i in 1:s){
  for(t in 2:timestep){
    mu[i,t-1] = inv_logit(A_t[t-1]+b1*Y_C[i,t-1]+C_f[i]);
}}

}

model {
  //Occurence
  psi_o~uniform(0, 10);
  b1~normal(0, 0.1); //bio interaction effects on occurence
  sigmaplot_o~uniform(0, 0.1);
  
  for(i in 1:10){
        C_f[i] ~ normal(0,tauplot_o);
  }
      
  for (i in 1:(timestep-1)){
        A_t[i] ~ normal(0, 0.1);
  }
  // likelihood 
  for(i in 1:s){
  for(t in 2:timestep){
    Y_O[i,t]~bernoulli(mu[i,t-1]);
    //Y_O[i,t]~bernoulli_logit(A_t[t]+b1*Y_C[i,t]+C_f[i]);
    target += log(mu[i,t-1]);
  }}
  
}
