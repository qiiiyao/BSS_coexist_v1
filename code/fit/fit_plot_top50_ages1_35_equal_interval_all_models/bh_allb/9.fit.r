rm(list=ls())

# set up R environment
library(cmdstanr)
options(mc.cores = parallel::detectCores(logical = F)) 
library(dplyr)
library(parallel)
library(doParallel)
setwd("/data/home/shpli3/R_projects/BSS_exclude_tree_raw")

# Load data
load("/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/fit_fp_top50_ages1_35_equal_interval.RData")

# run model chains separately then combine then to avoid memory issues
stan.seed = 52
file = file.path('/data/home/shpli3/my_pc/BSS/code/fit/stancode/model_comprison/constrained_competition_equal_interval/BH_allb_allposi.stan')
mod = cmdstan_model(file
#, stanc_options = list("allow-undefined")
)
numer_cores = parallel::detectCores(logical = F)
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10


fit_save = function(x) {
  #x = 2
  species = sp[x]
  fit = mod$sample(data = datalist[[x]],
                   seed = stan.seed, 
                   chains = chains, 
                   parallel_chains = chains,
                   iter_warmup = warmup,
                   iter_sampling = sampling,
                   refresh = 100, # print update every 100 iters
                   thin = thin,
                   max_treedepth = 20,
                   adapt_delta = 0.99)
  fit$loo()
  # Raw output
  #------------ 
  fit$save_object(file = paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb', '/posterior/',
                                unique(growD$f_p),'_', species, '.RDS'))
}

#--------------------------------------------------
# Estimate interactions with a joint NDD*RI model for t1
#--------------------------------------------------
data_t = fit_fp_top50_ages1_35_equal_interval
plot_list = sapply(seq(1, 480, 48), function(x){y = c(x:(x+47))})
names(fit_fp_top50_ages1_35_equal_interval)[plot_list]

fitted_plots_1 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb/summary')
m = regexec("[0-9]+\\_[0-9]+", fitted_plots_1)
fitted_plots = unlist(regmatches(fitted_plots_1, m))
unfitted = setdiff(names(fit_fp_top50_ages1_35_equal_interval)[plot_list[,9]], fitted_plots)
fitted_order = which(names(fit_fp_top50_ages1_35_equal_interval) %in% fitted_plots)
unfitted_order = which(names(fit_fp_top50_ages1_35_equal_interval) %in% unfitted)

for (i in unfitted_order) {
  #i = 1
  datalist = data_t[[i]][[3]]
  growD = data_t[[i]][[2]]
  unique(growD$f_p)
  sp = colnames(growD)[10:ncol(growD)]
  
   {cl = makeCluster(numer_cores-2)      
      registerDoParallel(cl)       
      foreach(x=1:length(datalist), .packages = c('dplyr', 'cmdstanr')) %dopar% {fit_save(x)
      }

    # Save results 
    summary_l = foreach(i=1:length(datalist)) %dopar% {
    #i = 2
    species = sp[i]
    fit = readRDS(paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb', '/posterior/',
                     unique(growD$f_p),'_', species, '.RDS'))
    summary = as.data.frame(fit$summary())}

    save(summary_l, file=paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb',
                                  '/summary/','summary_',
                                  unique(growD$f_p), '.rdata'))

    # Extract results
    r_trans_sp = vector()
    a_trans_sp = matrix(NA,length(sp), length(sp))
    r_trans_sp = foreach(i=1:length(datalist), .combine = c, .packages = c("dplyr")) %dopar% {(summary_l[[i]] %>% filter(variable == 'r'))$mean}
    b_trans_sp = foreach(i=1:length(datalist), .combine = c, .packages = c("dplyr")) %dopar% {(summary_l[[i]] %>% filter(variable == 'b'))$mean}
    a_trans_sp = foreach(i=1:length(datalist), .combine = rbind, .packages = c("dplyr")) %dopar% {(summary_l[[i]] %>% filter(grepl("a\\[", variable)))$mean}
    a_trans_sp_inver = a_trans_sp # for calculating ND and FD
    colnames(a_trans_sp_inver) = sp 
    rownames(a_trans_sp_inver) = sp
    colnames(a_trans_sp) = sp # saving the orginal results
    rownames(a_trans_sp) = sp
    names(r_trans_sp) = sp
  save(a_trans_sp, file=paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb', '/parameters/','a_',
                               unique(growD$f_p), '.rdata'))
  save(r_trans_sp, file=paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb', '/parameters/','r_',
                               unique(growD$f_p), '.rdata'))
  save(b_trans_sp, file=paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb', '/parameters/','b_',
                               unique(growD$f_p), '.rdata'))                             
  stopCluster(cl)
}


  ### calculate ND FD
  niche.diff_ij = matrix(nrow=length(sp), ncol=length(sp))
  fitness.diff_ij = matrix(nrow=length(sp), ncol=length(sp))
  fitness.diff_ji = matrix(nrow=length(sp), ncol=length(sp))
  
  for (i in 1:length(sp)) {
    
    # alphas
    for (j in i+1:(length(sp))){
      if(j > (length(sp))) break()
      
      aii = a_trans_sp_inver[sp[i], sp[i]]
      aij = a_trans_sp_inver[sp[i], sp[j]] ###
      aji = a_trans_sp_inver[sp[j], sp[i]] ###
      ajj = a_trans_sp_inver[sp[j], sp[j]]
      ri = r_trans_sp[i]
      rj = r_trans_sp[j]
      bi = b_trans_sp[i]
      bj = b_trans_sp[j]
      
      niche.over = sqrt( (aij*aji)/(aii*ajj) )
      niche.diff = 1 - niche.over
      
      fitness.diff.ij = (((exp(ri)^(1/bi))-1)/((exp(rj)^(1/bj))-1)) * sqrt( (ajj*aji)/(aii*aij) )
      
      fitness.diff.ji = (((exp(rj)^(1/bj))-1)/((exp(ri)^(1/bi))-1)) * sqrt( (aii*aij)/(ajj*aji) )
      
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
    alpha_intra_1 = c(as.character(sp[i]), a_trans_sp[sp[i], sp[i]])
    alpha_intra_1 = as.data.frame(alpha_intra_1)
    alpha_intra_1 = t(alpha_intra_1)
    alpha_intra = rbind(alpha_intra, alpha_intra_1)
  }
  
  intra_1 = cbind(alpha_intra, r_trans_sp, b_trans_sp)
  plot = rep(unique(growD$Plot), length(sp))
  field = rep(unique(growD$Field), length(sp))
  f_p = rep(unique(growD$f_p), length(sp))
  
  intra_trans = cbind(intra_1[, 1], 
                      field = field,
                      plot = plot,
                      f_p = f_p,
                      intra_1[,2: ncol(intra_1)])
  
  colnames(intra_trans)[c(1,5)] = c('sp', 'aii')
  
  save(intra_trans, file =paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb', '/parameters/','intra_',
                                 unique(growD$f_p), '.rdata'))
  
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
  nf_pair_all = cbind(nf_pair_all[,1], plot = rep(unique(growD$Plot), length(nf_pair_all$V1)),
                      field = rep(unique(growD$Field), length(nf_pair_all$V1)),
                      f_p = rep(unique(growD$f_p), length(nf_pair_all$V1)),
                      nf_pair_all[,2:ncol(nf_pair_all)])
  
  colnames(nf_pair_all)[c(1, 5, 6, 7, 8)] = c("sp_pair","species_i", "species_j", "nd", "fd")
  
  ###inter aij aji
  
  alpha_pair_ij = data.frame()
  
  for (i in 1:length(sp)) {
    for (j in i+1:length(sp)) {
      if (j  > length(sp)) break()
      alpha_pair_ij_1 = c(as.character(sp[i]),
                          as.character(sp[j]),
                          a_trans_sp[sp[i], sp[i]],
                          a_trans_sp[sp[j], sp[j]],
                          mean(c(a_trans_sp[sp[j], sp[j]],
                                 a_trans_sp[sp[i], sp[i]])),
                          a_trans_sp[sp[i], sp[j]])
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
                          a_trans_sp[sp[j], sp[j]],
                          a_trans_sp[sp[i], sp[i]],
                          mean(c(a_trans_sp[sp[j], sp[j]],
                                 a_trans_sp[sp[i], sp[i]])),
                          a_trans_sp[sp[j], sp[i]])
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
  save(inter_all_trans, file=paste0('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb', '/parameters/','inter_',
                                    unique(growD$f_p), '.rdata'))
  
}
