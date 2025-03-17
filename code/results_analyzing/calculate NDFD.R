
rm(list=ls())

library('cmdstanr')
library("coda")
library("deSolve")
library("bayesplot")
library("parallel")
library("posterior")
library("dplyr")

ls_a = list.files('D:/BSS study/results/parameters', pattern = '^a', full.names = T)
ls_r = list.files('D:/BSS study/results/parameters', pattern = '^r', full.names = T)

# Compile outside of loop
#load('D:/BSS study/code/invasion_stage.Rdata')
for (o in 1:length(ls_a)) {
      #o = 9
      a = read.table(ls_a[o], header = T) 
      f_p = strsplit(strsplit(strsplit(ls_a[o], '/')[[1]][5], '.txt')[[1]][1], 'a_')[[1]][2]
      r = read.table(ls_r[o], header = T)  
      ### calculate ND FD
      for (z in seq(1, 480, by = 1)) {
        #z = 3
      growD_Domin = re_cover_ab_domin_l[[z]]
      growD_estab = re_cover_ab_estab_l_f[[z]]
      growD_intro = re_cover_ab_intro_l_f[[z]]
      growD_nati = re_cover_ab_nati_l_f[[z]]
      
      if (unique(growD_Domin$f_p) == f_p){
        a_trans_sp_inver = a
        r_trans_sp = r
        growD_domin = growD_Domin
        growD_Estab = growD_estab
        growD_Intro = growD_intro
        growD_Nati = growD_nati
        
        sp_estab = colnames(growD_Estab)[8:ncol(growD_Estab)]
        sp_intro = colnames(growD_Intro)[8:ncol(growD_Intro)]
        sp_Nati = colnames(growD_Nati)[8:ncol(growD_Nati)]
        
        growD = cbind(growD_Estab, growD_Intro[,8:ncol(growD_Intro)], growD_Nati[,8:ncol(growD_Nati)])
        sp = c(sp_estab, sp_intro, sp_Nati)
        colnames(growD)[8:ncol(growD)] = sp
      }
        }
      niche.diff_ij = matrix(nrow=length(sp), ncol=length(sp))
      fitness.diff_ij = matrix(nrow=length(sp), ncol=length(sp))
      fitness.diff_ji = matrix(nrow=length(sp), ncol=length(sp))
      r_trans_sp = unlist(r_trans_sp)
      
      sp_domin = colnames(growD_domin)[8:ncol(growD_domin)]
      sp_pure_estab = setdiff(sp_estab, sp_domin)
      
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
      if(is.na(sp_domin[1])){ print('no domiant')
      } else if (sum(intra_trans$sp %in% sp_domin) == 0) {print('no domiant')
          } else {
      intra_trans[intra_trans$sp %in% sp_domin,]$stage = 'dominant'}
      colnames(intra_trans)[c(1,6)] = c('sp', 'aii')
      
      write.table(intra_trans, file=paste('D:/BSS study/results/parameters/intra_',unique(growD_domin$f_p), ".txt",sep = ''))
      
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
      sp_pair = paste(nf_pair[,1] , nf_pair[,2], sep = '.')
      nf_pair_all = cbind(sp_pair, nf_pair)
      nf_pair_all = cbind(nf_pair_all[,1], plot = rep(unique(growD$Plot), length(nf_pair_all$V1)),
                          field = rep(unique(growD$Field), length(nf_pair_all$V1)),
                          f_p = rep(unique(growD$f_p), length(nf_pair_all$V1)),
                          nf_pair_all[,2:ncol(nf_pair_all)])
      colnames(nf_pair_all)[c(1, 5, 6, 7, 8)] = c("sp_pair","species_i", "species_j", "nd", "fd")
      nf_pair_all$stage_i = NA
      nf_pair_all$stage_j = NA
      if(is.na(sp_domin[1])){ print('no domiant')
        } else if (sum(nf_pair_all$species_i %in% sp_domin) == 0) {print('no domiant')
      } else {
      nf_pair_all[nf_pair_all$species_i %in% sp_domin,]$stage_i = 'dominant'}
      nf_pair_all[nf_pair_all$species_i %in% sp_pure_estab,]$stage_i = 'establish'
      nf_pair_all[nf_pair_all$species_i %in% sp_intro,]$stage_i = 'introduce'
      nf_pair_all[nf_pair_all$species_i %in% sp_Nati,]$stage_i = 'native'
      
      if(is.na(sp_domin[1])){ print('no domiant')
        } else if (sum(nf_pair_all$species_j %in% sp_domin) == 0) {print('no domiant')
      } else {
      nf_pair_all[nf_pair_all$species_j %in% sp_domin,]$stage_j = 'dominant'}
      nf_pair_all[nf_pair_all$species_j %in% sp_pure_estab,]$stage_j = 'establish'
      nf_pair_all[nf_pair_all$species_j %in% sp_intro,]$stage_j = 'introduce'
      nf_pair_all[nf_pair_all$species_j %in% sp_Nati,]$stage_j = 'native' 
      unique(nf_pair_all$stage_i)
      nf_pair_all = nf_pair_all %>% relocate(stage_i, .after = species_j)
      nf_pair_all = nf_pair_all %>% relocate(stage_j, .after = stage_i)
      
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
      write.table(inter_all_trans, file=paste('D:/BSS study/results/parameters/inter_',unique(growD_domin$f_p), ".txt",sep = ''))
      #write.csv(inter_all_trans, file=paste('D:/BSS study/results/parameters/inter_',unique(growD$f_p), ".csv",sep = ''))
} 

