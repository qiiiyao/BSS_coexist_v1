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
load("code/data preparation/transformed data/fit_fp_ages1_35_top50_mergeall.RData")

# run model chains separately then combine then to avoid memory issues
stan.seed = 52
file_bh = file.path("D:/R projects/BSS/code/model_comprison/BH_random_time_loo.stan") 
mod_bh = cmdstan_model(file_bh)
file_logricker = file.path("D:/R projects/BSS/code/model_comprison/Ricker_log(n)_random_time_loo.stan") 
mod_logricker = cmdstan_model(file_logricker)
file_ricker = file.path("D:/R projects/BSS/code/model_comprison/Ricker_random_time_loo.stan") 
mod_ricker = cmdstan_model(file_ricker)
file_null = file.path("D:/R projects/BSS/code/model_comprison/exp_r_random_time.stan") 
mod_null = cmdstan_model(file_null)
file_bh_partb = file.path("D:/R projects/BSS/code/model_comprison/BH_partb_random_time_loo.stan") 
mod_bh_partb = cmdstan_model(file_bh_partb)
file_bh_allb = file.path("D:/R projects/BSS/code/model_comprison/BH_allb_random_time_loo.stan") 
mod_bh_allb = cmdstan_model(file_bh_allb)
file_bh_a_exp = file.path("D:/R projects/BSS/code/model_comprison/BH_n^a_random_time.stan") 
mod_bh_a_exp = cmdstan_model(file_bh_a_exp)

numer_cores = parallel::detectCores(logical = F)
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10

#--------------------------------------------------
# Estimate interactions with a joint NDD*RI model for t1
#--------------------------------------------------
data_t = fit_fp_ages1_35_top50_mergeall

  #i = 1
  datalist = data_t[[3]]
  growD = data_t[[2]]
  sp = colnames(growD)[2:ncol(growD)]
  out_bayes_compare.waic = list() #list for for loop output
  out_bayes_compare.loo = list() #list for for loop output
  
    for(x in c(2:length(datalist))){
      #x = 1
      species = sp[x]

      data = list(N = datalist[[x]]$N,
                  P = datalist[[x]]$P,
                  X = datalist[[x]]$X,
                  ftime = datalist[[x]]$ftime,
                  time = datalist[[x]]$time,
                  time_diff = datalist[[x]]$time_diff,
                  G = datalist[[x]]$G)
      fit_bh = mod_bh$sample(data = data,
                       seed = stan.seed, 
                       chains = chains, 
                       parallel_chains = chains,
                       iter_warmup = warmup,
                       iter_sampling = sampling,
                       refresh = 100, # print update every 100 iters
                       thin = thin,
                       max_treedepth = 20,
                       adapt_delta = 0.99)
      fit_logricker = mod_logricker$sample(data = data,
                             seed = stan.seed, 
                             chains = chains, 
                             parallel_chains = chains,
                             iter_warmup = warmup,
                             iter_sampling = sampling,
                             refresh = 100, # print update every 100 iters
                             thin = thin,
                             max_treedepth = 20,
                             adapt_delta = 0.99)
      fit_ricker = mod_ricker$sample(data = data,
                                     seed = stan.seed, 
                                     chains = chains, 
                                     parallel_chains = chains,
                                     iter_warmup = warmup,
                                     iter_sampling = sampling,
                                     refresh = 100, # print update every 100 iters
                                     thin = thin,
                                     max_treedepth = 20,
                                     adapt_delta = 0.99)
      fit_null = mod_null$sample(data = data,
                             seed = stan.seed, 
                             chains = chains, 
                             parallel_chains = chains,
                             iter_warmup = warmup,
                             iter_sampling = sampling,
                             refresh = 100, # print update every 100 iters
                             thin = thin,
                             max_treedepth = 20,
                             adapt_delta = 0.99)
      fit_bh_partb = mod_bh_partb$sample(data = data,
                                 seed = stan.seed, 
                                 chains = chains, 
                                 parallel_chains = chains,
                                 iter_warmup = warmup,
                                 iter_sampling = sampling,
                                 refresh = 100, # print update every 100 iters
                                 thin = thin,
                                 max_treedepth = 20,
                                 adapt_delta = 0.99)
      fit_bh_allb = mod_bh_allb$sample(data = data,
                                         seed = stan.seed, 
                                         chains = chains, 
                                         parallel_chains = chains,
                                         iter_warmup = warmup,
                                         iter_sampling = sampling,
                                         refresh = 100, # print update every 100 iters
                                         thin = thin,
                                         max_treedepth = 20,
                                         adapt_delta = 0.99)
      fit_bh_a_exp = mod_bh_a_exp$sample(data = data,
                                       seed = stan.seed, 
                                       chains = chains, 
                                       parallel_chains = chains,
                                       iter_warmup = warmup,
                                       iter_sampling = sampling,
                                       refresh = 100, # print update every 100 iters
                                       thin = thin,
                                       max_treedepth = 20,
                                       adapt_delta = 0.99)
      
      
      df.loo = as.data.frame(loo_compare(
        list(
          m1 = fit_bh$loo(),
          m2 = fit_logricker$loo(),
          #m3 = fit_ricker$loo(),
          m4 = fit_null$loo(),
          m5 = fit_bh_partb$loo(),
          m6 = fit_bh_allb$loo(),
          m7 = fit_bh_a_exp$loo()
        )
      )
      )
      df.loo$models = rownames(df.loo)
      
      df.waic = as.data.frame(loo_compare(
        list(
          m1 = waic(fit_bh$draws('log_lik')),
          m2 = waic(fit_logricker$draws('log_lik')),
          #m3 = waic(seed_bayes_fit_3),
          m4 = waic(fit_null$draws('log_lik')),
          m5 = waic(fit_bh_partb$draws('log_lik')),
          m6 = waic(fit_bh_allb$draws('log_lik')),
          m7 = waic(fit_bh_a_exp$draws('log_lik'))
        )
      )
      )
      
      df.waic$models = rownames(df.waic)
      
      out_bayes_compare.waic[[x]] = df.waic
      out_bayes_compare.loo[[x]] = df.loo
      
      # Raw output
      #------------ 
      #fit$save_object(file = paste0('results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35',
                                    #'/posterior/', species, '.RDS'))
    }
    
    save(out_bayes_compare.waic,
    file = 'results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35/summary/out_bayes_compare.waic.RDS')
    save(out_bayes_compare.loo, 
    file = 'results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35/summary/out_bayes_compare.loo.RDS')
    
    
    # calculate delta waic for model comparisons
    for(i in 1:length(out_bayes_compare.waic)){
      out_bayes_compare.waic[[i]]$models = rownames(out_bayes_compare.waic[[i]])
      out_bayes_compare.waic[[i]]$del.waic = out_bayes_compare.waic[[i]]$waic - min(out_bayes_compare.waic[[i]]$waic)
    }
    
    comparisons = do.call(rbind.data.frame, out_bayes_compare.loo)
    meancomps = comparisons %>% group_by(models) %>% summarise_all(mean)
    
    
    # Save results 
    summary_l = foreach(i=1:length(datalist)) %dopar% {
      #i = 3
      species = sp[i]
      fit = readRDS(paste0('results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35', '/posterior/',
                            species, '.RDS'))
      summary = fit$summary()
      }
    save(summary_l,
         file = 'results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35/summary/summary.rdata')
    
    loo_l = lapply(1:length(datalist), function(x){
      #i = x
      species = sp[x]
      fit = readRDS(paste0('results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35', '/posterior/',
                           species, '.RDS'))
      loo = fit$loo()
      if(length(which(loo$diagnostics$pareto_k > 0.5)) < 1) {
        print(paste0(species, "'s some pareto_k > 0.5"))
      }
      })
    
    names(loo_l) = sp
    
    save(loo_l,
         file = 'results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35/summary/loo.rdata')
    ### Test the pareto_k
    lapply(sp, function(x){
      loo = loo_l[[sp[x]]]
      if(length(which(loo$diagnostics$pareto_k > 0.5)) > 0) {
        print(paste0(species, "'s some pareto_k > 0.5"))
      }
    })
    ### All pareto_k is fine
    
    # Extract results
    r_trans_sp = vector()
    a_trans_sp = matrix(NA,length(sp), length(sp))
    r_trans_sp = foreach(i=1:length(datalist), .combine = c, .packages = c("dplyr")) %dopar% {(summary_l[[i]] %>% filter(variable == 'r'))$mean}
    a_trans_sp = foreach(i=1:length(datalist), .combine = rbind, .packages = c("dplyr")) %dopar% {(summary_l[[i]] %>% filter(grepl("a\\[", variable)))$mean}
    colnames(a_trans_sp) = sp
    rownames(a_trans_sp) = sp
    save(a_trans_sp, file=paste0('results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35',
                                       '/parameters/','a', '.rdata'))
    save(r_trans_sp, file=paste0('results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35',
                                 '/parameters/','r', '.rdata'))

  
  
  ### calculate ND FD
  niche.diff_ij = matrix(nrow=length(sp), ncol=length(sp))
  fitness.diff_ij = matrix(nrow=length(sp), ncol=length(sp))
  fitness.diff_ji = matrix(nrow=length(sp), ncol=length(sp))
  
  
  for (i in 1:length(sp)) {
    
    # alphas
    for (j in i+1:(length(sp))){
      if(j > (length(sp))) break()
      
      aii = a_trans_sp[sp[i], sp[i]]
      aij = a_trans_sp[sp[i], sp[j]] ###
      aji = a_trans_sp[sp[j], sp[i]] ###
      ajj = a_trans_sp[sp[j], sp[j]]
      
      niche.over = sqrt( ((aij/r_trans_sp[i])*(aji/r_trans_sp[j]))/((aii/r_trans_sp[i])*(ajj/r_trans_sp[j])) )
      niche.diff = 1 - niche.over
      
      fitness.diff.ij = sqrt( ((ajj/r_trans_sp[j])*(aji/r_trans_sp[j]))/((aii/r_trans_sp[i])*(aij/r_trans_sp[i])) )
      
      fitness.diff.ji = sqrt( ((aii/r_trans_sp[i])*(aij/r_trans_sp[i]))/((ajj/r_trans_sp[j])*(aji/r_trans_sp[j])) )
      
      niche.diff_ij[i, j] = niche.diff
      fitness.diff_ij[i, j] = fitness.diff.ij
      fitness.diff_ji[j, i] = fitness.diff.ji
      
    }
  }
  
  colnames(niche.diff_ij) = sp
  rownames(niche.diff_ij) = sp
  colnames(fitness.diff_ij) = sp
  rownames(fitness.diff_ij) = sp
  colnames(fitness.diff_ji) = sp
  rownames(fitness.diff_ji) = sp
  niche.diff_ijtrans_sp = niche.diff_ij
  fitness.diff_ijtrans_sp = fitness.diff_ij
  fitness.diff_jitrans_sp = fitness.diff_ji
  
  
  ##intra(r aii)
  
  r_trans_sp = as.data.frame(r_trans_sp)
  r_trans_sp$r_trans_sp = as.numeric(r_trans_sp$r_trans_sp)
  
  alpha_intra = data.frame()
  
  for (i in 1:length(sp)){
    if (i > length(sp)) break()
    alpha_intra_1 = c(as.character(sp[i]), a_trans_sp[sp[i], sp[i]]/(r_trans_sp$r_trans_sp[i]))
    alpha_intra_1 = as.data.frame(alpha_intra_1)
    alpha_intra_1 = t(alpha_intra_1)
    alpha_intra = rbind(alpha_intra, alpha_intra_1)
  }
  
  intra_1 = cbind(alpha_intra, r_trans_sp)
  
  colnames(intra_1)[1:2] = c('sp', 'aii')
  
  save(intra_1, file =paste0('results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35',
                                 '/parameters/','intra',
                                  '.rdata'))
  
  ### nf pair arranged
  
  nf_pair_1 = data.frame()
  
  for (i in 1:length(sp)) {
    for (j in i+1:length(sp)){
      if(j > length(sp)) break()
      
      nf_pair_ij = c(as.character(sp[i]), as.character(sp[j]), niche.diff_ij[sp[i],sp[j]], fitness.diff_ij[sp[i],sp[j]])
      nf_pair_ij = as.data.frame(nf_pair_ij)
      nf_pair_ij = t(nf_pair_ij)
      nf_pair_1 = rbind(nf_pair_1, nf_pair_ij)
      
    }
  }
  
  nf_pair_2 = data.frame()
  
  for (i in 1:length(sp)) {
    for (j in i+1:length(sp)){
      if(j > length(sp)) break()
      
      nf_pair_ji = c(as.character(sp[j]), as.character(sp[i]), niche.diff_ij[sp[i],sp[j]], fitness.diff_ji[sp[j],sp[i]])
      nf_pair_ji = as.data.frame(nf_pair_ji)
      nf_pair_ji = t(nf_pair_ji)
      nf_pair_2 = rbind(nf_pair_2, nf_pair_ji)
      
    }
  }
  
  nf_pair = rbind(nf_pair_1, nf_pair_2)
  sp_pair = paste(nf_pair[,1] , nf_pair[,2], sep = '.')
  nf_pair_all = cbind(sp_pair, nf_pair)
  nf_pair_all = cbind(nf_pair_all[,1],
                      nf_pair_all[,2:ncol(nf_pair_all)])
  
  colnames(nf_pair_all) = c("sp_pair","species_i", "species_j", "nd", "fd")
  
  ###inter aij aji
  
  alpha_pair_ij = data.frame()
  
  for (i in 1:length(sp)) {
    for (j in i+1:length(sp)) {
      if (j  > length(sp)) break()
      alpha_pair_ij_1 = c(as.character(sp[i]),
                          as.character(sp[j]),
                          a_trans_sp[sp[i], sp[i]]/(r_trans_sp$r_trans_sp[i]),
                          a_trans_sp[sp[j], sp[j]]/(r_trans_sp$r_trans_sp[j]),
                          mean(c(a_trans_sp[sp[i], sp[i]], a_trans_sp[sp[j], sp[j]])),
                          a_trans_sp[sp[i], sp[j]]/(r_trans_sp$r_trans_sp[i]))
      alpha_pair_ij_1 = as.data.frame(alpha_pair_ij_1)
      alpha_pair_ij_1 = t(alpha_pair_ij_1)
      alpha_pair_ij = rbind(alpha_pair_ij, alpha_pair_ij_1)
      
    }
  }
  
  alpha_pair_ji = data.frame()
  
  for (i in 1:length(sp)) {
    for (j in i+1:length(sp)) {
      if (j > length(sp)) break()
      alpha_pair_ji_1 = c(as.character(sp[j]),
                          as.character(sp[i]),
                          a_trans_sp[sp[j], sp[j]]/(r_trans_sp$r_trans_sp[j]),
                          a_trans_sp[sp[i], sp[i]]/(r_trans_sp$r_trans_sp[i]),
                          mean(c(a_trans_sp[sp[j], sp[j]]/(r_trans_sp$r_trans_sp[j]),
                                 a_trans_sp[sp[i], sp[i]]/(r_trans_sp$r_trans_sp[i]))),
                          a_trans_sp[sp[j], sp[i]]/(r_trans_sp$r_trans_sp[i]))
      alpha_pair_ji_1 = as.data.frame(alpha_pair_ji_1)
      alpha_pair_ji_1 = t(alpha_pair_ji_1)
      alpha_pair_ji = rbind(alpha_pair_ji, alpha_pair_ji_1)
      
    }
  }
  
  alpha_pair = rbind(alpha_pair_ij, alpha_pair_ji)
  sp_pair_alpha = paste(alpha_pair[,1], alpha_pair[,2], sep = '_')
  alpha_pair_all = cbind(sp_pair_alpha, alpha_pair)
  colnames(alpha_pair_all) = c('sp_pair', 'sp_i', 'sp_j', 'aii', 'ajj', 'a_mean', 'aij')
  
  ### aggragate inter 
  
  library('dplyr')
  nf_pair_all = arrange(nf_pair_all, nf_pair_all$sp_pair)
  alpha_pair_all = arrange(alpha_pair_all, alpha_pair_all$sp_pair)
  
  inter_all_trans = cbind(nf_pair_all,
                          alpha_pair_all[, c(4:7)])  
  #colnames(inter_all_trans)
  save(inter_all_trans, file=paste0('results/fit_results/BSS_exclude_trees_raw/plot_mergetop50_into1_ages1_35',
                                    '/parameters/','inter_',
                                    unique(growD$f_p), '.rdata'))
  
  
