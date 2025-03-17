# Running the joint model on simulated data 

# set up R environment
rm(list = ls())
library(cmdstanr)
options(mc.cores = parallel::detectCores()) 

library(rethinking)
library(reshape2)
library(bayesplot)
library(posterior)
library(data.table)


# load required functions
source('code/fit/functions/data_prep.R')

setwd("D:/R projects/BSS")
load("code/data preparation/transformed data/fit_fp_top50_ages1_35.RData")

file = file.path("D:/R projects/BSS/code/fit/stancode/RIM_random_time.stan")
mod_cmd = cmdstan_model(file)
file_ode_ricker = file.path("D:/R projects/BSS/code/fit/stancode/ode_Ricker.stan")
mod_ode_ricker = cmdstan_model(file_ode_ricker)
stan.seed = 1234
numer_cores = parallel::detectCores(logical = F)
chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 10

# identify focal and neighbouring species to be matched to parameter estimates
data_t = fit_fp_top50_ages1_35
plot_list = sapply(seq(1, 480, 48), function(x){y = c(x:(x+47))})

for (i in plot_list) {
  #i = 1
df = fit_fp_top50_ages1_35[[i]][[2]]
sps = colnames(df)[10:ncol(df)]
f_p = unique(df$f_p)
df_l = lapply(fit_fp_top50_ages1_35[[i]][[3]], function(x){y = cbind(
  focal = gsub('[0-9]+',
               '',
               names(x$G)),
  time = x$time,
  time_diff = x$time_diff,
  G = x$G, 
  x$X)})

df = as.data.frame(rbindlist(df_l))
focalID = unique(df$focal)  # this should return the names of unique focal groups in the order
# in which they are encountered in the dataframe
neighbourID = colnames(df[,-c(1:3)])

# prepare the data into the format required by STAN and the model code
stan.data = data_prep(perform = 'G', 
                      focal = 'focal', 
                      nonNcols = 3, # number of columns that aren't neighbour abundances
                      df = df)

# Run the model! 
fit_cmd = mod_cmd$sample(data = stan.data,# named list of data 
                         seed = stan.seed, 
                         chains = chains, 
                         parallel_chains = chains,
                         iter_warmup = warmup,
                         iter_sampling = sampling,
                         refresh = 100, # print update every 100 iters
                         thin = thin,
                         max_treedepth = 20,
                         adapt_delta = 0.99
)
fit$save_object(file = paste0('fit_results/plot_ages1_35_top50', '/posterior/',
                              f_p, '.RDS'))

# Save results 
summary = fit_cmd$summary()

save(summary, file=paste0('fit_results/plot_ages1_35_top50',
                            '/summary/','summary_',
                          f_p, '.rdata'))

# Extract results
r_trans_sp = (summary %>% filter(grepl("r\\[", variable)))$mean
a_trans_sp = matrix((summary %>% filter(grepl("a_rim\\[", variable)))$mean,
                    length(sp), length(sp))
colnames(a_trans_sp) = sp
rownames(a_trans_sp) = sp
save(a_trans_sp, file=paste0('fit_results/plot_ages1_35_top50', '/parameters/','a_',
                             f_p, '.rdata'))
save(r_trans_sp, file=paste0('fit_results/plot_ages1_35_top50', '/parameters/','r_',
                             f_p, '.rdata'))



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

intra_trans = cbind(intra_1[, 1], 
                    field = field,
                    plot = plot,
                    f_p = f_p,
                    intra_1[,2: ncol(intra_1)])

colnames(intra_trans)[c(1,5)] = c('sp', 'aii')

save(intra_trans, file =paste0('fit_results/plot_ages1_35_top50', '/parameters/','intra_',
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
save(inter_all_trans, file=paste0('fit_results/plot_ages1_35_top50', '/parameters/','inter_',
                                  unique(growD$f_p), '.rdata'))



}
Sys.time()




