//NDDM (the Ricker model) & RIM joint model (don't permit positive interaction!)///
// modified from Bimler_et_al_2023_MEE
data {
  int<lower=1> S;          // total number of focal species (s in the manuscript) 
  int<lower=1> N;          // total number of observations (rows in model matrix) (n in the manuscript)
  int<lower=0> T;          // total number of interaction partners (columns in model matrix) (t in the manuscript)
  int<lower=0> I;          // total number of identifiable interactions 
  
  array[N] int species_ID;   // index matching observations to focal species
  array[N] real perform;      // response variable: per captia growth rate (in the manuscript)
    
  array[I] int icol;  // indices matching pairwise inferrable to location in interaction matrix
  array[I] int irow;  
  int<lower=0> ftime; // number of years in our data
  int<lower=0> time_diff[N]; // year difference in our data
  matrix[N,T] X;         // neighbour abundances (the model matrix)
  array[N] int time; // years in our data
} 


transformed data {
  // Define transformed data variables here
  vector[N] log_perform; // Example transformation: Log-transform of observations

  // Compute the transformations
  for (n in 1:N) {
    log_perform[n] = log(perform[n]);
  }
}

parameters {
  
  vector<lower=0>[S] gamma_i;    // species-specific intrinsic growth rate (must positive) 
    
  vector<lower=0>[S] sigma; // species-specific error
    
  vector<lower=0>[I] beta_ij;     // vector of interactions which are inferrable by the NDDM 
  vector[ftime] ran_t;     // vector of random time intercepts
  real<lower=0> weight;    // weighting value controlling the average strength of interactions (must be positive)
  vector<lower=0>[S] raw_response; // species-specific effect parameter (res in the manuscript but without weight)
  vector<lower=0>[T] raw_effect;   // species-specific effect parameter (eff in the manuscript)

} 

transformed parameters {
  
  // transformed res and eff to unit vectors
  vector[S] response;
  {
    real norm = sqrt(dot_self(raw_response)); // 
    response = raw_response / norm; // 
  }
  
  vector[T] effect;
  {
    real norm = sqrt(dot_self(raw_effect)); // 
    effect = raw_effect / norm; // 
  }
  
  vector[N] mu;              // the RIM linear predictor for perform (here per captia growth rate)
  vector[N] mu2;             // the joint model linear predictor for perform (here per captia growth rate)
  
  matrix[S, T] ri_betaij;   // interaction matrix for response - impact estimates (res*eff in the manuscript)
  matrix[S, T] ndd_betaij;  // interaction matrix for joint model interaction estimates (alpha in the manuscript)
  
  // get RIM interactions
  ri_betaij = weight * response*effect';  // the apostrophe transposes the effect vector
   
  // Estimate response-impact interactions
  for(n in 1:N) {

       mu[n] = (gamma_i[species_ID[n]] - dot_product(X[n], ri_betaij[species_ID[n], ]) + ran_t[time[n]])*time_diff[n];  
  }
  
  // NDDM estimates identifiable interactions, and uses RIM estimates when non-identifiable:
  ndd_betaij = ri_betaij; // initialise nddm interaction matrix to rim estimates
  // match identifiable interactions parameters to the correct position in the interaction matrix
  for(i in 1:I) {
    ndd_betaij[irow[i], icol[i]] = beta_ij[i];
  }

  // estimate identifiable interactions
  for(n in 1:N) {

        mu2[n] = (gamma_i[species_ID[n]] - dot_product(X[n], ndd_betaij[species_ID[n], ]) + ran_t[time[n]])*time_diff[n];  
  
  }
   
} 

model {

  // priors
  gamma_i ~ cauchy(0, 10);   // prior for the intercept following Gelman 2008
  sigma ~ cauchy(0, 1); 
  beta_ij ~ normal(0, 1);    // prior for interactions inferred by the NDDM
  weight ~ normal(0, 10); // constrained by parameter definition to be positive
  ran_t ~ normal(0, 1);
  // no prior needed for response or effect as we can use the default prior for the unit_vector

  // maximising the likelihood for the joint model 
  for(n in 1:N) {
  
    // maximise the likelihood of the RIM for all observations
    log_perform[n] ~ normal(mu[n], sigma[species_ID[n]]);
    // maximise the likelihood of the NDDM over identifiable interactions
    target += normal_lpdf(log_perform[n] | mu2[n], sigma[species_ID[n]]);
    
  }
  
}

generated quantities {

  vector[N] log_lik_rim;	// log-likelihood for the RIM
  vector[N] log_lik_nddm;	// log-likelihood for the joint model
  
  for(n in 1:N) {	
    log_lik_rim[n] = normal_lpdf(log_perform[n] | mu[n], sigma[species_ID[n]]);
	
    log_lik_nddm[n] = normal_lpdf(log_perform[n] | mu2[n], sigma[species_ID[n]]);
  }
  
}
