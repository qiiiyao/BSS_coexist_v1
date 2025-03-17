############################################################################################################
# This code calculate pair-wise low-density growth rate calculations using the full posterior probabilities
# for each parameter
############################################################################################################
rm(list=ls())

# set up R environment
library(cmdstanr)
options(mc.cores = parallel::detectCores(logical = F)) 
library(dplyr)
library(parallel)
library(doParallel)
library(data.table)

setwd("D:/R projects/BSS")

# Load data
load("code/data preparation/transformed data/fit_fp_top50_ages1_35.RData")

#-----------------------------------------------------------------------
# calculate stabilizing

do.stabilization = function(alpha.ii, alpha.jj, alpha.ij, alpha.ji) {
  rho = sqrt((alpha.ij*alpha.ji)/(alpha.jj*alpha.ii))
  stabilizing = 1- rho
  return(stabilizing)
}

# calculate average fitness differences
do.kj.over.ki = function(lambda.i, lambda.j, alpha.ii, alpha.jj, alpha.ij, alpha.ji) {
  n.i = lambda.i
  n.j = lambda.j
  kj.over.ki = (n.j/n.i)*sqrt((alpha.ij*alpha.ii)/(alpha.jj*alpha.ji))
  #kj.over.ki = sqrt((alpha.ij*alpha.ii)/(alpha.jj*alpha.ji))
  return(kj.over.ki)
}

# calculate coexistence
do.coexistence = function(stabilizing, kj.over.ki) {
  rho = 1 - stabilizing
  ki.over.kj = 1/kj.over.ki
  if (rho < ki.over.kj & (1/rho) > ki.over.kj) {
    outcome = data.frame(coexists = 1, priority = 0, exclusion = 0)
  } else if (rho > ki.over.kj & (1/rho) < ki.over.kj) {
    outcome = data.frame(coexists = 0, priority = 1, exclusion = 0)
  } else {
    outcome = data.frame(coexists = 0, priority = 0, exclusion = 1)
  }
  return(outcome)
}

runs = 1000 # posterior draws
# -----------------------------------------------------------------------------------------
# Make things generic

calculate_coe_prob_all = function(x){
  x=fit_fp_top50_ages1_35[[1]]
  dat = x
  growD = dat[[2]]
  f_p = unique(growD$f_p)
  sp_list = colnames(growD)[10:ncol(growD)]
  
  # Create the lists of intrinsic growth rates and competition coefficients
  list_files = list.files('results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50/posterior',
                          pattern = paste0(f_p, '\\_.'))
  posterior_l = lapply(list_files, function(x){
    #x = list_files[1]
    y = readRDS(paste0('results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50/posterior/',
                       x))
    posterior = as.data.frame(t(y$draws(format = 'data.frame')))
    posterior$variable = rownames(posterior)
    posterior_a = posterior %>% filter(grepl('^a\\[', variable)) %>% 
      as.data.frame()
    posterior_r = posterior %>% filter(grepl('^r$', variable)) %>% 
      as.data.frame()
    posterior_l1 = list(posterior_a = posterior_a,
                        posterior_r = posterior_r)
    
    return(posterior_l1)
  })
  posterior_a_l = lapply(1:runs, function(run){
    #x = 1
    a_run1 = sapply(posterior_l, function(z){
      a = z[['posterior_a']]
      a1 = a[,run]
    })
    a_run1 = t(a_run1)
    colnames(a_run1) = sp_list
    rownames(a_run1) = sp_list
    return(a_run1)
  })
  
  posterior_r_l = lapply(1:runs, function(run){
    #run = 1
    r_run1 = sapply(posterior_l, function(z){
      r = z[['posterior_r']]
      r1 = r[,run]
    })
    names(r_run1) = sp_list
    return(r_run1)
  })
  posterior_r_a_l = lapply(1:runs, function(x){
    r_a = list(r = posterior_r_l[[x]],
               a = posterior_a_l[[x]])
    return(r_a)
  })
  
  total_combos = combn(sp_list, 2)
  
  #p=r # set x to number of runs which equals number of posterior values 
  # -----------------------------------------------------------------------------------------
  # Run through all possible combos
  sp_pair_prob_l = lapply(posterior_r_a_l, function(r_a){
    #r_a = posterior_r_a_l[[1]]
    r_all = r_a[['r']]
    a_all = r_a[['a']]
    co_dat = apply(total_combos, 2, function(x){
      #x = total_combos[,1]
      lambdas.j = r_all[x[1]]
      intras.j = a_all[x[1], x[1]]
      
      lambdas.i = r_all[x[2]]
      intras.i = a_all[x[2], x[2]]
      
      inter.ij = a_all[x[1], x[2]]
      inter.ji = a_all[x[2], x[1]]
      
      stabilizing = do.stabilization(alpha.ii=intras.i,
                                     alpha.jj=intras.j, 
                                     alpha.ij=inter.ij,
                                     alpha.ji=inter.ji)
      
      equalizing = do.kj.over.ki(lambda.i=lambdas.i,
                                 lambda.j=lambdas.j, 
                                 alpha.ii=intras.i,
                                 alpha.jj=intras.j, 
                                 alpha.ij=inter.ij,
                                 alpha.ji=inter.ji)
      
      coexists = do.coexistence(stabilizing=stabilizing,
                                kj.over.ki=equalizing)[1, 'coexists'] # negs in both of these
      prioritys = do.coexistence(stabilizing=stabilizing,
                                 kj.over.ki=equalizing)[1, 'priority'] # negs in both of these
      exclusions = do.coexistence(stabilizing=stabilizing,
                                  kj.over.ki=equalizing)[1, 'exclusion'] # negs in both of these
      dat = data.frame(stabilizing = stabilizing,
                       equalizing = equalizing,
                       coexists = coexists,
                       prioritys = prioritys,
                       exclusions = exclusions)
      return(dat)
    })
    co_dat = rbindlist(co_dat)
    
    return(co_dat)
  })
  
  stabilizings = sapply(sp_pair_prob_l, function(x){
    stab = as.numeric(unlist(x[,'stabilizing']))
  })
  equalizings = sapply(sp_pair_prob_l, function(x){
    equa = as.numeric(unlist(x[,'equalizing']))
  })
  coexists = sapply(sp_pair_prob_l, function(x){
    coex = as.numeric(unlist(x[,'coexists']))
  })
  prioritys = sapply(sp_pair_prob_l, function(x){
    prio = as.numeric(unlist(x[,'prioritys']))
  })
  exclusions = sapply(sp_pair_prob_l, function(x){
    excl = as.numeric(unlist(x[,'exclusions']))
  })
  
  coexist_prob = apply(coexists, 1,
                       function(x){length(x[which(x==1)])/length(x)})
  priority_prob = apply(prioritys, 1,
                        function(x){length(x[which(x==1)])/length(x)})
  exclusion_prob = apply(exclusions, 1,
                         function(x){length(x[which(x==1)])/length(x)})
  
  sp_pair_prob = as.data.frame(cbind(f_p = f_p, t(total_combos),
                               coexist_prob = coexist_prob,
                               priority_prob = priority_prob,
                               exclusion_prob = exclusion_prob))
  
  posterior_nd_fd = list(nd = stabilizings, fd = equalizings)
  file.path_1 = paste0('results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50/summary/',
                     f_p, '_posterior_nd_fd.rdata')
  file.path_2 = paste0('results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50/summary/',
                       f_p, '_posterior_coex_probs.rdata')
  save(posterior_nd_fd, file = file.path_1)
  save(sp_pair_prob, file = file.path_2)
  return(sp_pair_prob)
  
}

# Make things generic
cl = makeCluster(4)      
registerDoParallel(cl) 

sp_pair_prob_l = foreach(x=fit_fp_top50_ages1_35,
                         .packages = c('dplyr', 'cmdstanr', 'data.table')) %dopar% {
                           calculate_coe_prob_all(x)}
stopCluster(cl)



list_files = list.files('D:/R projects/BSS/results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50/summary',
                        pattern = '._nd_fd.rdata')

sp_pair_prob_l = lapply(list_files, function(x){
  x = list_files
  load(paste0('D:/R projects/BSS/results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50/summary/',
              x))
  f_p = regmatches(x, gregexec('[0-9]+\\_[0-9]+', x))[[1]][1]
  f_p_posi = which(names(fit_fp_top50_ages1_35) == f_p)
  dat=fit_fp_top50_ages1_35[[f_p_posi]]
  growD = dat[[2]]
  sp_list = colnames(growD)[10:ncol(growD)]
  total_combos = combn(sp_list, 2)
  
  stabilizing = posterior_nd_fd$nd
  
  equalizing = posterior_nd_fd$fd
  
  needed_probs = mapply(do.coexistence, stabilizing, equalizing,
                        SIMPLIFY = 'matrix') # 16-18
  coexists = matrix(needed_probs['coexists',], ncol = ncol(stabilizing),
                    nrow = nrow(stabilizing), byrow = F) # negs in both of these
  prioritys = matrix(needed_probs['priority',], ncol = ncol(stabilizing),
                     nrow = nrow(stabilizing), byrow = F)
  exclusions = matrix(needed_probs['exclusion',], ncol = ncol(stabilizing),
                      nrow = nrow(stabilizing), byrow = F)
  
  coexist_prob = apply(coexists, 1,
                       function(x){length(x[which(x==1)])/length(x)})
  priority_prob = apply(prioritys, 1,
                        function(x){length(x[which(x==1)])/length(x)})
  exclusion_prob = apply(exclusions, 1,
                         function(x){length(x[which(x==1)])/length(x)})
  
  sp_pair_prob = data.frame(f_p = f_p, t(total_combos),
                            coexist_prob = coexist_prob,
                            priority_prob = priority_prob,
                            exclusion_prob = exclusion_prob)
})

