# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript
### load work dictionary
library('cmdstanr')
library("coda")
library("deSolve")
library("BayesianTools")
library("parallel")
library("posterior")
library("dplyr")
library('tidyr')
library('data.table')

PrelimFit_l = list()
load("code/data preparation/transformed data/fit_field_alltime_same_0.20_raw.RData")
load("code/data preparation/transformed data/fit_fp_alltime_same_0.20_raw.RData")
BSScover.list = fit_fp_alltime_same_0.20_raw[[i]][[5]]
focal_dat = data.frame()

for (i in 1:length(BSScover.list)) {
  i = 1

BSScover.dat = BSScover.list[[i]]

SpMatrix = BSScover.dat$X
Fecundity = BSScover.dat$G
SpNames = colnames(SpMatrix)
Intra = ifelse(SpNames == SpNames[i], 1, 0)
S = ncol(SpMatrix)
DataVec = list(T = nrow(SpMatrix), S = ncol(SpMatrix), Fecundity = Fecundity,
                 SpMatrix = SpMatrix, Intra = Intra, tau0 = 8,
                 slab_scale = sqrt(2), slab_df = 4)
  
  # Now run a perliminary fit of the model to assess parameter shrinkage
  mod_Prelim = cmdstan_model(stan_file = 'code/fit/stancode/Mod_hassel_lehman.stan')
  PrelimFit = mod_Prelim$sample(data = DataVec, seed = 123, 
                                       chains = 4, 
                                       parallel_chains = 4,
                                       iter_warmup = 1000,
                                       iter_sampling = 2500,
                                       refresh = 10, # print update every 10 iters
                                       max_treedepth = 15,
                                       adapt_delta = 0.99)
  print(PrelimFit$summary(), n = 50)
PrelimPosteriors = PrelimFit$draws(format = 'df')
colnames(PrelimPosteriors)

Inclusion_ij = matrix(data = 0, nrow = 1, ncol = S)
IntLevel = 0.5

for(s in 1:S){
  #s = 1
  Ints_ij = HDInterval::hdi(as.numeric(unlist(PrelimPosteriors[,paste('alpha_hat_ij[',s,']',sep = '')])),
                            credMass = IntLevel)
  if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
    Inclusion_ij[s] = 1
  }
}

sum(Inclusion_ij)
Inclusion_ij = as.numeric(Inclusion_ij)
 if (sum(Inclusion_ij) > 0) {
   focal_dat = rbind(focal_dat, data.frame(plot = z, sp = i))
 }
DataVec = list(T = nrow(SpMatrix), S = ncol(SpMatrix), Fecundity = Fecundity,
                SpMatrix = SpMatrix, Intra = Intra, Inclusion_ij = Inclusion_ij)
mod_Final = cmdstan_model(stan_file = 'code/fit/stancode/Mod_hassel_lehman_final.stan')
FinalFit = mod_Final$sample(data = DataVec, seed = 123, 
                              chains = 4, 
                              parallel_chains = 4,
                              iter_warmup = 1000,
                              iter_sampling = 2500,
                              refresh = 10, # print update every 10 iters
                              max_treedepth = 15,
                              adapt_delta = 0.99)
print(FinalFit$summary(), n = 25)
FinalPosteriors <- extract(FinalFit)

}






