
data {
  int<lower=0> s; // number of plots
  int<lower=0> timestep; //number of timesteps = ?
  real Y_G[s,timestep]; // species focal cover data
  real X_F[s,timestep]; // species focal cover data
  real X_E[s,timestep]; // species exotic cover data
  real X_O[s,timestep]; // species all other cover data
  
}


parameters {
  vector[timestep] A_t;
  vector[s] C_f;
  real<lower=0> r;
  real b_i; //bio interaction effects on occurence
  real b_j; //bio interaction effects on occurenc
  real b_o; //bio interaction effects on occurenc
  real sigmaplot_g;
  real<lower=0> tau;
 
}

transformed parameters{
  real<lower=0> tauplot_g;
  tauplot_g = sqrt(1/(sigmaplot_g*sigmaplot_g));
  array[s,timestep] real mu;
  
   for(i in 1:s){
  for(t in 1:timestep){
    mu[i,t] = r-X_F[i,t]*b_i-X_E[i,t]*b_j-X_O[i,t]*b_o+A_t[t]+C_f[i];
}}

}

model {
  //growth
  b_i~normal(0, 0.001); //bio interaction effects from intra
  b_j~normal(0, 0.001); //bio interaction effects from inter
  b_o~normal(0, 0.001); //bio interaction effects from all other inter
  r~normal(0, 0.001); //intrinsic growth rate
  sigmaplot_g~uniform(0, 0.1);
  tau~cauchy(0, 5);
  
  for(i in 1:10){
        C_f[i] ~ normal(0,tauplot_g);
  }
      
  for (i in 1:(timestep)){
        A_t[i] ~ normal(0, 0.001);
  }
  // likelihood 
  for(i in 1:s){
  for(t in 2:timestep){
    Y_G[i,t]~normal(mu[i,t],tau);

  }}
  
}
