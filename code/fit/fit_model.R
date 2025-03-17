# rm(list=ls())

# set up R environment
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = F)) 
library(dplyr)
library(parallel)
library(doParallel)


setwd('/data/home/shpli3/R_projects/BSS')
# Load data
load("/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/fit_stage_field_t_all_filter_1.rdata")

# run model chains separately then combine then to avoid memory issues
stan.seed <- 52
mod = stan_model(file = 'code/fit/fit_filter_1/Ricker_random_time.stan',
                 verbose = T)
numer_cores = parallel::detectCores(logical = F)

fit_save = function(x) {
  #x = 1
  species = sp[x]
  fit = sampling(mod,
                 data = datalist[[x]],seed = stan.seed,  
                 chains = 4,
                 warmup = 1000,          # number of warmup iterations per chain
                 iter = 3000,            # total number of iterations per chain
                 refresh = 100,         # show progress every 'refresh' iterations
                 control = list(max_treedepth = 20,
                                adapt_delta = 0.99)
  )
  
  # Raw output
  #------------
  save(fit, file = paste0('equal_moving_win_field_filter1/','equal_t', 
                          z, '/posterior/',
                          unique(growD$Field),'_', species, '.rdata')) # model fit
  
  #res$summary()
}

#--------------------------------------------------
# Estimate interactions with a joint NDD*RI model for t1
#--------------------------------------------------
for(z in 1:length(fit_stage_field_t_all_filter_1)){
#z = 1
data_t = fit_stage_field_t_all_filter_1[[z]]

field = c(1:length(data_t))

for (i in field) {
  #i = 9
  datalist = data_t[[i]][[5]]
  growD = data_t[[i]][[4]]
  sp = colnames(growD)[7:ncol(growD)]
  { cl = makeCluster(numer_cores-2)      
    registerDoParallel(cl)       
    foreach(x=1:length(datalist), .packages = c('dplyr', 'rstan')) %dopar% {fit_save(x)}
    
    # Save results 
    summary_l = foreach(i=1:length(datalist)) %dopar% {
      #i = 1
      species = sp[i]
      load(paste0('equal_moving_win_field_filter1/','equal_t', 
                  z, '/posterior/',
                  unique(growD$Field),'_', species, '.rdata'))
      summary = as.data.frame(summary(fit)[[1]])}
    
    save(summary_l, file=paste0('equal_moving_win_field_filter1/','equal_t', 
                                z, '/summary/','summary_',
                                unique(growD$Field), '.rdata'))
    
    # Extract results
    r_trans_sp = vector()
    a_trans_sp = matrix(NA,length(sp), length(sp))
    r_trans_sp = foreach(i=1:length(datalist), .combine = c, .packages = c("dplyr")) %dopar% {(summary_l[[i]] %>% filter(rownames(summary_l[[i]]) == 'r'))$mean}
    a_trans_sp = foreach(i=1:length(datalist), .combine = rbind, .packages = c("dplyr")) %dopar% {(summary_l[[i]] %>% filter(grepl("a\\[", rownames(summary_l[[i]]))))$mean}
    a_trans_sp_inver = a_trans_sp*-1
    colnames(a_trans_sp_inver) = sp
    rownames(a_trans_sp_inver) = sp
    save(a_trans_sp_inver, file=paste0('equal_moving_win_field_filter1/','equal_t', 
                                       z, '/parameters/','a_',
                                       unique(growD$Field), '.rdata'))
    save(r_trans_sp, file=paste0('equal_moving_win_field_filter1/','equal_t', 
                                 z, '/parameters/','r_',
                                 unique(growD$Field), '.rdata'))
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
    alpha_intra_1 = c(as.character(sp[i]), a_trans_sp_inver[sp[i], sp[i]]/(r_trans_sp$r_trans_sp[i]))
    alpha_intra_1 = as.data.frame(alpha_intra_1)
    alpha_intra_1 = t(alpha_intra_1)
    alpha_intra = rbind(alpha_intra, alpha_intra_1)
  }
  
  intra_1 = cbind(alpha_intra, r_trans_sp)
  plot = rep(unique(growD$Plot), length(sp))
  field = rep(unique(growD$Field), length(sp))
  f_p = rep(unique(growD$Field), length(sp))
  
  intra_trans = cbind(intra_1[, 1], 
                      field = field,
                      plot = plot,
                      f_p = f_p,
                      t_win = z,
                      intra_1[,2: ncol(intra_1)])
  
  colnames(intra_trans)[c(1,6)] = c('sp', 'aii')
  
  save(intra_trans, file =paste0('equal_moving_win_field_filter1/','equal_t', 
                                 z, '/parameters/','intra_',
                                 unique(growD$Field), '.rdata'))
  
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
                      f_p = rep(unique(growD$Field), length(nf_pair_all$V1)),
                      t_win = z,
                      nf_pair_all[,2:ncol(nf_pair_all)])
  colnames(nf_pair_all)
  colnames(nf_pair_all)[c(1, 6, 7, 8, 9)] = c("sp_pair","species_i", "species_j", "nd", "fd")
  
  ###inter aij aji
  
  alpha_pair_ij = data.frame()
  
  for (i in 1:length(sp)) {
    for (j in i+1:length(sp)) {
      if (j  > length(sp)) break()
      alpha_pair_ij_1 = c(as.character(sp[i]),
                          as.character(sp[j]),
                          a_trans_sp_inver[sp[i], sp[i]]/(r_trans_sp$r_trans_sp[i]),
                          a_trans_sp_inver[sp[j], sp[j]]/(r_trans_sp$r_trans_sp[j]),
                          mean(c(a_trans_sp_inver[sp[i], sp[i]], a_trans_sp_inver[sp[j], sp[j]])),
                          a_trans_sp_inver[sp[i], sp[j]]/(r_trans_sp$r_trans_sp[i]))
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
                          a_trans_sp_inver[sp[j], sp[j]]/(r_trans_sp$r_trans_sp[j]),
                          a_trans_sp_inver[sp[i], sp[i]]/(r_trans_sp$r_trans_sp[i]),
                          mean(c(a_trans_sp_inver[sp[j], sp[j]]/(r_trans_sp$r_trans_sp[j]),
                                 a_trans_sp_inver[sp[i], sp[i]]/(r_trans_sp$r_trans_sp[i]))),
                          a_trans_sp_inver[sp[j], sp[i]]/(r_trans_sp$r_trans_sp[i]))
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
  save(inter_all_trans, file=paste0('equal_moving_win_field_filter1/','equal_t', 
                                    z, '/parameters/','inter_',
                                    unique(growD$Field), '.rdata'))
  unique(inter_all_trans$nd)
  }
}
