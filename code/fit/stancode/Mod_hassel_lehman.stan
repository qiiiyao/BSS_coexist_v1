data{
  int<lower = 1> T; // Number of Timestep
  int<lower = 1> S; // Number of species
  real Fecundity[T]; // Fecundity of the focal species in each plot
  //int Fecundity[T]; // Fecundity of the focal species in each plot
  #int reserve[T];   // Indicator variable for the reserve each plot is located in
  matrix[T,S] SpMatrix; // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  //vector[T] env;   // Environmental values for each plot
  int<lower = 0> Intra[S]; // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
  // The below values define the regularized horseshoe priors used for species-specific parameters
  real tau0; 		// determines the scale of the global shrinkage parameter (tau)
  real slab_scale;	// scale for significant alpha_sp values
  real slab_df;		// effective degrees of freedom for significant alpha_sp values
}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  real lambdas;
  real alpha_generic_tilde;
  real alpha_intra_tilde;
  vector[S] alpha_hat_ij_tilde;
  vector<lower = 0>[S] local_shrinkage_ij;
  real<lower = 0> c2_tilde;
  //real<lower = 0> tau_tilde;
  real<lower = 0> theta;
}

transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  real c2;
  real tau;
  vector[S] alpha_hat_ij;
  vector[S] local_shrinkage_ij_tilde;
  real alpha_generic;
  real alpha_intra;

  tau = tau0; //tau = tau0*tau_tilde; 	 
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017

    for(s in 1:S){
      local_shrinkage_ij_tilde[s] = sqrt( c2 * square(local_shrinkage_ij[s]) / (c2 + square(tau) * square(local_shrinkage_ij[s])) );
      alpha_hat_ij[s] = tau * local_shrinkage_ij_tilde[s] * alpha_hat_ij_tilde[s];
    }
  

  // scale the lambdas and alphas values
  alpha_generic = 1 * alpha_generic_tilde - 0;
  alpha_intra = 3 * alpha_intra_tilde - 6;
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*T values for each species.
  vector[T] F_hat;
  vector[T] interaction_effects;
  row_vector[S] alpha_ij;
  real lambda_i;
  
  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  lambdas ~ normal(0, 1);


  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
  
    alpha_hat_ij_tilde ~ normal(0,1);
    local_shrinkage_ij ~ cauchy(0,1);
    
    theta ~ normal(0, 1);
    //tau_tilde ~ cauchy(0,1);
    c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  // implement the biological model
  for(i in 1:T){
    lambda_i = exp(lambdas);
    for(s in 1:S){
        alpha_ij[s] = exp((1-Intra[s]) * alpha_generic + Intra[s] * alpha_intra + (1-Intra[s]) * alpha_hat_ij[s]);
    }
    
    interaction_effects[i] = sum(alpha_ij.* SpMatrix[i,]);
    F_hat[i] = lambda_i-interaction_effects[i];
  }
  Fecundity ~ normal(F_hat, theta);
  //Fecundity ~ poisson(F_hat);
}

