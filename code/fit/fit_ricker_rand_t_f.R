rm(list=ls())

# set up R environment
library(cmdstanr)
options(mc.cores = parallel::detectCores(logical = F)) 
library(dplyr)
library(parallel)
library(doParallel)
require(loo)
setwd('D:/R projects/BSS')

# Load data
load("code/data preparation/transformed data/fit_ranfield_top50_ages1_35_sp.RData")

# run model chains separately then combine then to avoid memory issues
stan.seed = 52
file_ricker_ran_t_f = file.path("code/fit/stancode/models_comparison/Ricker_random_t_f.stan") 
mod_ricker_ran_t_f = cmdstan_model(file_ricker_ran_t_f)

numer_cores = parallel::detectCores(logical = F)
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10

for (i in 1:50) {
#i = 1
data = fit_ranfield_top50_ages1_35_sp[[i]]
species = names(fit_ranfield_top50_ages1_35_sp)[i]
fit_ricker_ran_t_f = mod_ricker_ran_t_f$sample(data = data,
                       seed = stan.seed, 
                       chains = chains, 
                       parallel_chains = chains,
                       iter_warmup = warmup,
                       iter_sampling = sampling,
                       refresh = 100, # print update every 100 iters
                       thin = thin,
                       max_treedepth = 20,
                       adapt_delta = 0.99)

fit_ricker_ran_t_f$save_object(file = paste0('fit_results/ran_t_f_ages1_35_top50', '/posterior/fit_',
                                             species, '.RDS'))

}
