
rm(list=ls())

library('cmdstanr')
library("coda")
library("deSolve")
library("bayesplot")
library("parallel")
library("posterior")
library("dplyr")

# Compile outside of loop
#load('D:/BSS study/code/invasion_stage.Rdata')
for (z in 1:length(cover_prop_domin_l)) {
  z = 1
  growD_Estab = re_cover_prop_estab_l_f[[z]]
  growD_Intro = re_cover_prop_intro_l_f[[z]]
  growD_Nati = re_cover_prop_nati_l_f[[z]]
  sp_estab = colnames(growD_Estab)[8:ncol(growD_Estab)]
  sp_intro = colnames(growD_Intro)[8:ncol(growD_Intro)]
  sp_Nati = colnames(growD_Nati)[8:ncol(growD_Nati)]
  
  growD = cbind(growD_Estab, growD_Intro[,8:ncol(growD_Intro)], growD_Nati[,8:ncol(growD_Nati)])
  sp = c(sp_estab, sp_intro, sp_Nati)
  colnames(growD)[8:ncol(growD)] = sp
  growD = arrange(growD, growD$Age)
  growD[growD==0] = 1e-4
  Pred = growD[2:(nrow(growD)),1:ncol(growD)]
  Yrs = nrow(Pred)
  yid = as.numeric(as.factor(Pred$Age))
  sp_matrix = growD[1:(nrow(growD)-1),8:ncol(growD)]
  S = ncol(sp_matrix)
  
  summary_l = list()
  for (j in 1:length(sp)) {
  #j = 1
  Y = unlist(Pred %>% select(sp[j]))
  I = unlist(sp_matrix %>% select(sp[j]))
  G = log(Y/I)
  mod_dat_trans_sp = cbind(G, sp_matrix)
  
  sdscal = sd(residuals(lm(paste('G',
                                 paste(sp, collapse = "+"),
                                 sep = '~'),
                                 data = mod_dat_trans_sp)))
  
  datalist = list(N = Yrs, P = ncol(sp_matrix), ftime=Yrs, time=yid,
                  G = G, X=sp_matrix, sdscal = sdscal)
  
  mod_ricker = cmdstan_model(stan_file = 'D:/BSS study/code/hierarchical bayes lme_no ranr.stan')
  
  RickerFit = mod_ricker$sample(data = datalist, seed = 123, 
                               chains = 3, 
                               parallel_chains = 3,
                               iter_warmup = 1000,
                               iter_sampling = 2000,
                               refresh = 10, # print update every 10 iters
                               max_treedepth = 15,
                               thin = 10,
                               #init = inits,
                               adapt_delta = 0.99)
  
  summary_l[[j]] = RickerFit$summary()
  summary = summary_l[[j]]
  write.table(summary, paste('D:/BSS study/results/summary/',unique(growD$f_p), '_', sp[j], '.txt', sep = ''))
  }
  
  r_trans_sp = vector()
  a_trans_sp = matrix(NA,length(sp), length(sp))
  
  for (l in 1:length(sp)) {
    for (m in 1:length(sp)) {
      
      r = summary_l[[l]] %>% filter(variable == 'r')
      r_trans_sp[l] = r$mean
      a = summary_l[[l]] %>% filter(variable == paste("a[",m,"]",
                                                       sep = ""))
      a_trans_sp[l,m] = a$mean
    }
  }
  r_trans_sp
  a_trans_sp_inver = a_trans_sp*-1
  colnames(a_trans_sp_inver) = sp
  rownames(a_trans_sp_inver) = sp
  write.table(a_trans_sp_inver, file=paste('D:/BSS study/results/parameters/a_',
                                         unique(growD$f_p),
                                         ".txt",sep = ''))
  write.table(r_trans_sp, file=paste('D:/BSS study/results/parameters/r_',
                                   unique(growD$f_p),
                                   ".txt",sep = ''))
  
  
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
  f_p = rep(unique(growD$f_p), length(sp))
  stage = c(rep('establish', ncol(growD_Estab)-7),rep('introduce', ncol(growD_Intro)-7),rep('native', ncol(growD_Nati)-7))
  
  intra_trans = cbind(intra_1[, 1], 
                      field = field,
                      plot = plot,
                      f_p = f_p,
                      stage = stage,
                      intra_1[,2: ncol(intra_1)])
  intra_trans[intra_trans$`intra_1[, 1]` %in% colnames(growD_domin)[8:ncol(growD_domin)],]$stage = 'dominant'
  colnames(intra_trans) = c('sp', 'plot', 'dim', 'pH', 'aii', 'r')
  
  write.csv(intra_trans, file=paste('D:/BSS study/results/parameters/intra_posi_',unique(coexist_trans_sp$plot), "_sp.csv",sep = ''))
  
  ### nf pair arranged
  
  nf_pair_1 <- data.frame()
  
  for (i in 1:length(sp)) {
    for (j in i+1:length(sp)){
      if(j > length(sp)) break()
      
      nf_pair_ij <- c(as.character(sp[i]), as.character(sp[j]), niche.diff_ij[sp[i],sp[j]], fitness.diff_ij[sp[i],sp[j]])
      nf_pair_ij <- as.data.frame(nf_pair_ij)
      nf_pair_ij <- t(nf_pair_ij)
      nf_pair_1 <- rbind(nf_pair_1, nf_pair_ij)
      
    }
  }
  
  nf_pair_2 <- data.frame()
  
  for (i in 1:length(sp)) {
    for (j in i+1:length(sp)){
      if(j > length(sp)) break()
      
      nf_pair_ji <- c(as.character(sp[j]), as.character(sp[i]), niche.diff_ij[sp[i],sp[j]], fitness.diff_ji[sp[j],sp[i]])
      nf_pair_ji <- as.data.frame(nf_pair_ji)
      nf_pair_ji <- t(nf_pair_ji)
      nf_pair_2 <- rbind(nf_pair_2, nf_pair_ji)
      
    }
  }
  
  nf_pair = rbind(nf_pair_1, nf_pair_2)
  sp_pair = paste(nf_pair[,1] , nf_pair[,2], sep = '_')
  nf_pair_all = cbind( sp_pair, nf_pair)
  nf_pair_all = cbind(nf_pair_all[,1], plot = rep(unique(coexist_trans_sp$plot), length(nf_pair_all$V1)), nf_pair_all[,2:ncol(nf_pair_all)])
  colnames(nf_pair_all) = c("sp_pair",'plot', "species_i", "species_j", "nd", "fd")
  nf_pair_trans_sp = nf_pair_all
  
  ###inter aij aji
  
  alpha_pair_ij = data.frame()
  
  for (i in 1:length(sp)) {
    for (j in i+1:length(sp)) {
      if (j  > length(sp)) break()
      alpha_pair_ij_1 = c(as.character(sp[i]),
                          as.character(sp[j]),
                          a_trans_sp_inver[sp[i], sp[i]],
                          a_trans_sp_inver[sp[j], sp[j]],
                          mean(c(a_trans_sp_inver[sp[i], sp[i]], a_trans_sp_inver[sp[j], sp[j]])),
                          a_trans_sp_inver[sp[i], sp[j]])
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
                          a_trans_sp_inver[sp[j], sp[j]],
                          a_trans_sp_inver[sp[i], sp[i]],
                          mean(c(a_trans_sp_inver[sp[j], sp[j]],
                                 a_trans_sp_inver[sp[i], sp[i]])),
                          a_trans_sp_inver[sp[j], sp[i]])
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
  inter_all_trans = cbind(alpha_pair_all,
                          nf_pair_all[, c(2,5:6)],
                          pH = pH,
                          dim = dim)  
  write.csv(inter_all_trans, file=paste('D:/BSS study/results/parameters/inter_posi_',unique(coexist_trans_sp$plot), "_sp.csv",sep = ''))
  
  }
