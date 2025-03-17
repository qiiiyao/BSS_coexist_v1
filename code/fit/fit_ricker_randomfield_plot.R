# rm(list=ls())
# load packages
library(nlme)
library(dplyr)
library(data.table)
######## nlme frequentist method ##################

#perform the modelling with lme and temporal autocorrelation
load('code/data preparation/fit_ricker_randomspace.RData')
lCtr = lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

mlist = list()
for(i in 1:50){
  #i = 1
  mlist[[i]] = lme(as.formula(paste(yy[i], " ~ ", paste(top50.short, collapse="+"))),
                    data = pchange.all2, 
                    random=~1|Plot/Field,
                    control=lCtr, 
                    correlation=corAR1(form=~Yeart),
                    method='REML', na.action=na.omit)
}

summary(mlist[[1]])
mlist = mlist[!sapply(mlist,is.null)]
coef.list = lapply(mlist, function(x)summary(x)$coef$fixed)

inter.mat = matrix(nrow=length(yy)-1, ncol=length(yy)-1) #matrix of species interactions.
intrinsic.site.lui = vector() #matrix of intrinsic ability to growth and how LUI modifies it. 

for(i in 1:length(coef.list)){
  cc = coef.list[[i]] ## extract coefficients
  
  cc2 = cc[-1]
  cc3 = cc2[c(1:length(coef.list))]
  cc5 = cc[1]
  
  inter.mat[i,] = cc3
  intrinsic.site.lui[i]= cc5
}

yy2 = gsub("_delta", "", yy)
yy2 = yy2[-c(52)] # to remove name "rest"
diag(inter.mat) = diag(inter.mat) -1 #because assuming that the 45 degrees slope means no change.
row.names(inter.mat) = yy2
colnames(inter.mat) = yy2
row.names(lui.mat) = yy2
colnames(lui.mat) = yy2
row.names(intrinsic.site.lui) = yy2
colnames(intrinsic.site.lui) = c("Intrinsic", "LUI")

## Get the real species name
sp = sapply(strsplit(sp, '_'), function(x){
  paste(x[1], x[2], sep = '_') 
})

######## hierachial bayesian method ##################
# load packages
library(cmdstanr)
mod_ricker_1 = cmdstan_model(stan_file = 'code/Ricker_random_space.stan')

fit_trans = list()
summary_l = list()
for (i in 1:50) {
  
  data_list = list(N = length(pchange.all2$Poa_pratensis_delta), # level 1
                   Y = unlist((pchange.all2 %>% select(yy[i]))[,1]),
                   K = 50,
                   X = as.matrix(pchange.all2[,c(57:(57+49))]),
                   N_1 = 49,
                   M_1 = 1,
                   J_1 = dd_data$J_1,
                   Z_1_1 = rep(1, nrow(pchange.all2)),
                   N_2 = 480,
                   M_2 = 1,
                   J_2 = dd_data$J_2,
                   Z_2_1 = rep(1, nrow(pchange.all2)),
                   prior_only = 0)
  
  fit_trans[[i]] = mod_ricker_1$sample(data = data_list, seed = 123, 
                          parallel_chains = 4,
                          iter_warmup = 1000,
                          iter_sampling = 2500,
                          refresh = 500, # print update every 10 iters
                          max_treedepth = 15,
                          thin = 1,
                          cores = parallel::detectCores(),
                          #init = inits,
                          adapt_delta = 0.99)
  RickerFit = fit_trans[[i]]
  summary_l[[i]] = RickerFit$summary()
  posterior = RickerFit$draws(format = "df")
  summary = summary_l[[i]]
  fwrite(summary, paste('results/fit_results/summary_randomspace/',
                        yy[i], '.txt', sep = ''))
  write.csv(posterior, paste('results/fit_results/posterior_randomspace/',
                             yy[i], '.txt', sep = ''))
}

####### extract the parameters
library(data.table)
library(dplyr)
load('code/data preparation/fit_ricker_randomspace.RData')
## Extract r and a
sp = sp
r_trans_sp = vector()
a_trans_sp = matrix(NA,50,50)

for (i in 1:50) {
summary = fread(paste('results/fit_results/summary_randomspace/',
                           yy[i], '.txt', sep = ''), header = T)
  for (m in 1:50) {
      r = summary %>% filter(variable == 'r')
      r_trans_sp[i] = r$median
      a = summary %>% filter(variable == paste("a[",m,"]",
                                                      sep = ""))
      a_trans_sp[i,m] = a$median
  }
}

r_trans_sp_median = vector()
a_trans_sp_median = matrix(NA,50,50)

for (i in 1:50) {
  summary = fread(paste('results/fit_results/summary_randomspace/',
                        yy[i], '.txt', sep = ''), header = T)
  for (m in 1:50) {
    r = summary %>% filter(variable == 'r')
    r_trans_sp_median[i] = r$median
    a = summary %>% filter(variable == paste("a[",m,"]",
                                             sep = ""))
    a_trans_sp_median[i,m] = a$median
  }
}

r_trans_sp
a_trans_sp_inver = a_trans_sp*-1
colnames(a_trans_sp_inver) = sp
rownames(a_trans_sp_inver) = sp
write.table(a_trans_sp_inver, file=paste('results/fit_results/parameters_randomspace/a_randomspace',
                                         ".txt",sep = ''))
write.table(r_trans_sp, file=paste('results/fit_results/parameters_randomspace/r_randomspace',
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
    
    niche.over = sqrt(((aij/r_trans_sp[i])*(aji/r_trans_sp[j]))/((aii/r_trans_sp[i])*(ajj/r_trans_sp[j])))
    niche.diff = 1 - niche.over
    
    fitness.diff.ij = sqrt(((ajj/r_trans_sp[j])*(aji/r_trans_sp[j]))/((aii/r_trans_sp[i])*(aij/r_trans_sp[i])))
    
    fitness.diff.ji = sqrt(((aii/r_trans_sp[i])*(aij/r_trans_sp[i]))/((ajj/r_trans_sp[j])*(aji/r_trans_sp[j])))
    
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
  alpha_intra_1 = c(as.character(sp[i]), a_trans_sp_inver[sp[i], sp[i]])
  alpha_intra_1 = as.data.frame(alpha_intra_1)
  alpha_intra_1 = t(alpha_intra_1)
  alpha_intra = rbind(alpha_intra, alpha_intra_1)
}

intra_1 = cbind(alpha_intra, r_trans_sp)
plot = rep(unique(growD$Plot), length(sp))
field = rep(unique(growD$Field), length(sp))
f_p = rep(unique(growD$f_p), length(sp))
stage = c(rep('establish', length(sp_estab)),rep('introduce', length(sp_intro)),rep('native', length(sp_Nati)))

intra_trans = cbind(intra_1[, 1], 
                    field = field,
                    plot = plot,
                    f_p = f_p,
                    stage = stage,
                    intra_1[,2: ncol(intra_1)])
if (length(sp_domin) > 0){ 
  intra_trans[intra_trans$`intra_1[, 1]` %in% sp_domin,]$stage = 'dominant'
}
colnames(intra_trans)[c(1,6)] = c('sp', 'aii')

write.table(intra_trans, file=paste('D:/BSS study/results/all_prior/parameters/intra_',unique(growD$f_p), ".txt",sep = ''))

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
nf_pair_all$stage_i = NA
nf_pair_all$stage_j = NA

if (length(sp_estab) > 0) {
  nf_pair_all[nf_pair_all$species_i %in% sp_estab,]$stage_i = 'establish'
  nf_pair_all[nf_pair_all$species_j %in% sp_estab,]$stage_j = 'establish'
}
if (length(sp_intro) > 0) {
  nf_pair_all[nf_pair_all$species_i %in% sp_intro,]$stage_i = 'introduce'
  nf_pair_all[nf_pair_all$species_j %in% sp_intro,]$stage_j = 'introduce'
}
if (length(sp_Nati) > 0) {
  nf_pair_all[nf_pair_all$species_i %in% sp_Nati,]$stage_i = 'native'
  nf_pair_all[nf_pair_all$species_j %in% sp_Nati,]$stage_j = 'native'
}
if (length(sp_domin) > 0) {
  nf_pair_all[nf_pair_all$species_i %in% sp_domin,]$stage_i = 'dominant'
  nf_pair_all[nf_pair_all$species_j %in% sp_domin,]$stage_j = 'dominant'
}

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
colnames(alpha_pair_all) = c('sp_pair', 'sp_i', 'sp_j',
                             'aii', 'ajj', 'a_mean', 'aij')

### aggragate inter 

library('dplyr')
nf_pair_all = arrange(nf_pair_all, nf_pair_all$sp_pair)
alpha_pair_all = arrange(alpha_pair_all, alpha_pair_all$sp_pair)

inter_all_trans = cbind(nf_pair_all,
                        alpha_pair_all[, c(4:7)])  
write.table(inter_all_trans, file=paste('D:/BSS study/results/all_prior/parameters/inter_',unique(growD$f_p), ".txt",sep = ''))

#write.csv(inter_all_trans, file=paste('D:/BSS study/results/all_prior/parameters/inter_',unique(growD$f_p), ".csv",sep = ''))
#  } else { print('too low df!')}
#} else { print('need more data!') }
