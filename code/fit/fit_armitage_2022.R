### Dave Armitage (david.armitage@oist.jp)
### Code adapted from https://doi.org/10.5281/zenodo.7083314 by M. Van Dyke
### Please download the author's original data from the link above prior to running the analysis
### Last edit: 20 Dec 2022

# This script fits and compares Bayesian models of competition to the data provided in Van Dyke et al.
# Where possible, it is an exact replica of the bootstrap analysis performed in the original paper
rm(list = ls())
library(dplyr)
library(cmdstanr)
library(tidyr)
library(loo)

options(scipen = 5)
set.seed(1234)
setwd("D:/R projects/BSS")

### load work dictionary
file1 = file.path("code/model_comprison/exp_r_random_time.stan")
mod1 = cmdstan_model(file1)
file2 = file.path("code/model_comprison/Ricker_random_time_loo.stan") 
mod2 = cmdstan_model(file2)
#file2 = file.path("code/model_comprison/Ricker_loo.stan") 
#mod2 = cmdstan_model(file2)
file3 = file.path("code/model_comprison/Ricker_log(n+1)random_time_loo.stan")
mod3 = cmdstan_model(file3)
file4 = file.path("code/model_comprison/BH_random_time_loo.stan")
mod4 = cmdstan_model(file4)
file5 = file.path("code/model_comprison/BH_log(n+1)_random_time_loo.stan")
mod5 = cmdstan_model(file5)
file6 = file.path("code/model_comprison/BH_b_random_time_loo.stan")
mod6 = cmdstan_model(file6)
file7 = file.path("code/model_comprison/BH_1+b_random_time_loo.stan")
mod7 = cmdstan_model(file7)
file8 = file.path("code/model_comprison/BH_n^a_random_time.stan")
mod8 = cmdstan_model(file8)

# Load data
#load("code/data preparation/transformed data/fit_fp_same_ages_top50_early_suc.RData")
load("code/data preparation/transformed data/fit_fp_same_ages_mergetop50_early_suc.RData")

#Create a data frame with all combinations of treatment and species 
#spp_list = sort(na.omit( unique(seed_data$focal)))
#treat_list = sort( na.omit( unique(seed_data$Tr)))
#spp_combos = expand.grid(species = spp_list, treatment = treat_list)
# loop through the species list, fit the model, return the dataframe 
out_bayes_compare.waic = list() #list for for loop output
out_bayes_compare.loo = list() #list for for loop output
numer_cores = parallel::detectCores(logical = F)
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10

plot_list = sapply(seq(1, 480, 48), function(x){y = c(x:(x+47))})

for(i in 1:480){ 
   i = 1
   dat_fit = fit_fp_same_ages_mergetop50_early_suc[[i]][[3]]
   out_bayes_compare.waic_1 = list()
   out_bayes_compare.loo_1 = list()
   f_p = unique(fit_fp_same_ages_mergetop50_early_suc[[i]]$re_cover_ab$f_p)
  for (j in 1:length(dat_fit)) {
   j = 2
    
  data = list(N = dat_fit[[j]]$N,
                   P = dat_fit[[j]]$P,
                   X = dat_fit[[j]]$X,
                   ftime = dat_fit[[j]]$ftime,
                   time = dat_fit[[j]]$time,
                   G = exp(dat_fit[[j]]$G))

focal_sp = gsub("1", "", names(dat_fit[[j]]$G)[1])
background_sp = colnames(dat_fit[[j]]$X)

# Next, we fit the bayesian models for comparison purposes. 
# This code will fit a bunch of models so it can take some time to run
# It takes around 1 hour on an M1 chip with 8 cores in parallel.

growth_bayes_fit_1= mod1$sample(data = data, seed = 123,
                              chains = chains,
                              #  init = 50,
                              parallel_chains = chains,
                              iter_warmup = warmup,
                              iter_sampling = sampling,
                              refresh = 10, # print update every 10 iters
                              thin = thin,
                              max_treedepth = 20,
                              adapt_delta = 0.99)
growth_bayes_fit_2= mod2$sample(data = data, seed = 123,
                              chains = chains,
                              #  init = 50,
                              parallel_chains = chains,
                              iter_warmup = warmup,
                              iter_sampling = sampling,
                              refresh = 10, # print update every 10 iters
                              thin = thin,
                              max_treedepth = 20,
                              adapt_delta = 0.99)
growth_bayes_fit_3= mod3$sample(data = data, seed = 123,
                              chains = chains,
                              #  init = 50,
                              parallel_chains = chains,
                              iter_warmup = warmup,
                              iter_sampling = sampling,
                              refresh = 10, # print update every 10 iters
                              thin = thin,
                              max_treedepth = 20,
                              adapt_delta = 0.99)
growth_bayes_fit_4= mod4$sample(data = data, seed = 123,
                              chains = chains,
                              #  init = 50,
                              parallel_chains = chains,
                              iter_warmup = warmup,
                              iter_sampling = sampling,
                              refresh = 10, # print update every 10 iters
                              thin = thin,
                              max_treedepth = 20,
                              adapt_delta = 0.99)
growth_bayes_fit_5 = mod5$sample(data = data, seed = 123,
                               chains = chains,
                               #  init = 50,
                               parallel_chains = chains,
                               iter_warmup = warmup,
                               iter_sampling = sampling,
                               refresh = 10, # print update every 10 iters
                               thin = thin,
                               max_treedepth = 20,
                               adapt_delta = 0.99)
growth_bayes_fit_6 = mod6$sample(data = data, seed = 123,
                               chains = chains,
                               #  init = 50,
                               parallel_chains = chains,
                               iter_warmup = warmup,
                               iter_sampling = sampling,
                               refresh = 10, # print update every 10 iters
                               thin = thin,
                               max_treedepth = 20,
                               adapt_delta = 0.99)
growth_bayes_fit_7 = mod7$sample(data = data, seed = 123,
                               chains = chains,
                               #  init = 50,
                               parallel_chains = chains,
                               iter_warmup = warmup,
                               iter_sampling = sampling,
                               refresh = 10, # print update every 10 iters
                               thin = thin,
                               max_treedepth = 20,
                               adapt_delta = 0.99)
growth_bayes_fit_8 = mod8$sample(data = data, seed = 123,
                               chains = chains,
                               #  init = 50,
                               parallel_chains = chains,
                               iter_warmup = warmup,
                               iter_sampling = sampling,
                               refresh = 10, # print update every 10 iters
                               thin = thin,
                               max_treedepth = 20,
                               adapt_delta = 0.99)
  
  df.loo = as.data.frame(loo_compare(
    list(
      m1 = growth_bayes_fit_1$loo(),
      m2 = growth_bayes_fit_2$loo(),
      m3 = growth_bayes_fit_3$loo(),
      m4 = growth_bayes_fit_4$loo(),
      m5 = growth_bayes_fit_5$loo(),
      m6 = growth_bayes_fit_6$loo(),
      m7 = growth_bayes_fit_7$loo(),
      m8 = growth_bayes_fit_8$loo()
    )
  )
  )

  out_bayes_compare.loo_1[[j]] = df.loo

    }
   #save(out_bayes_compare.waic_1,
       # file = paste0("results/fit_results/BSS_exculde_trees_raw/plot_sameages_top50_early_suc/summary/out_bayes_compare",
    #                  f_p,".waic.rdata"))
   save(out_bayes_compare.loo_1,
        file = paste0("results/fit_results/BSS_exculde_trees_raw/plot_sameages_top50_early_suc/summary/out_bayes_compare",
                      f_p, ".loo.rdata"))
   
  }   




