rm(list = ls())

### Load packages
library(dplyr)
#library(betareg)
library(loo)
library(stringr)
library(data.table)
#library(cmdstanr)

### load data
load("/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/fit_fp_top50_ages1_35_equal_interval.RData")

setwd('/data/home/shpli3/R_projects/BSS_exclude_tree_raw')

### LOO for model comparison
plot_list = sapply(seq(1, 480, 48), function(x){y = c(x:(x+47))})
data_t = fit_fp_top50_ages1_35_equal_interval
posterior_posi = '/data/home/shpli3/R_projects/BSS_exclude_tree_raw/fit_results/plot_ages1_35_top50_equal_interval_model_comparison/'

for(i in 1:480){  ## Plot
    #i = 23
    out_bayes_compare.loo_1 = list()
    dat_fit = data_t[[i]][[3]]
    growD = data_t[[i]][[2]]
    f_p = unique(growD$f_p)
    sp_dat = colnames(growD)[10:ncol(growD)]

    for (j in 1:length(dat_fit)) { ## Species
    #j = 1
    sp = sp_dat[j]

    fit_null = readRDS(paste0(posterior_posi, 'null/posterior/',
                       f_p, '_', sp, '.RDS'))
    fit_ricker = readRDS(paste0(posterior_posi, 'ricker/posterior/',
                       f_p, '_', sp, '.RDS'))
    fit_ricker_logn = readRDS(paste0(posterior_posi, 'ricker_logn/posterior/',
                       f_p, '_', sp, '.RDS'))
    fit_ricker_logn1 = readRDS(paste0(posterior_posi, 'ricker_logn1/posterior/',
                       f_p, '_', sp, '.RDS'))
    fit_bh = readRDS(paste0(posterior_posi, 'bh/posterior/',
                       f_p, '_', sp, '.RDS'))
    fit_bh_allb = readRDS(paste0(posterior_posi, 'bh_allb/posterior/',
                       f_p, '_', sp, '.RDS'))
    fit_bh_partialb = readRDS(paste0(posterior_posi, 'bh_partialb/posterior/',
                       f_p, '_', sp, '.RDS'))
    #fit_bh_a = readRDS(paste0(posterior_posi, 'bh_^a/posterior/',
                #       f_p, '_', sp, '.RDS'))

    df.loo = as.data.frame(loo_compare( 
     list(
      null = fit_null$loo(),
      ricker = fit_ricker$loo(),
      ricker_logn = fit_ricker_logn$loo(),
      ricker_logn1 = fit_ricker_logn1$loo(),
      bh = fit_bh$loo()
      ,bh_allb = fit_bh_allb$loo()
      ,bh_partialb = fit_bh_partialb$loo()
      #,bh_a = fit_bh_a$loo()
    )
  )
  )

    df.loo$focal_species = sp
    df.loo$f_p = f_p
    out_bayes_compare.loo_1[[j]] = df.loo
    }
   #save(out_bayes_compare.waic_1,
       # file = paste0("results/fit_results/BSS_exculde_trees_raw/plot_sameages_top50_early_suc/summary/out_bayes_compare",
    #                  f_p,".waic.rdata"))
   save(out_bayes_compare.loo_1,
        file = paste0("/data/home/shpli3/R_projects/BSS_exclude_tree_raw/fit_results/plot_ages1_35_top50_equal_interval_model_comparison/loo_results/",
                      f_p, ".loo.rdata"))
   
  }   


### Calculating summary of these model validation
loo_summary_l = list.files('/data/home/shpli3/R_projects/BSS_exclude_tree_raw/fit_results/plot_ages1_35_top50_equal_interval_model_comparison/loo_results')
setwd('/data/home/shpli3/R_projects/BSS_exclude_tree_raw/fit_results/plot_ages1_35_top50_equal_interval_model_comparison/loo_results')
df_loo_l = lapply(loo_summary_l, function(x){
                         #x = loo_summary_l[1]
                         load(x)
                         out_bayes_compare.loo_2 = lapply(out_bayes_compare.loo_1,
                                                          function(z){
                                                          z$model = NA
                                                          z$model = rownames(z)
                                                          return(z)
                                                          })
                         
                         y = rbindlist(out_bayes_compare.loo_2)
                         
                         return(y)
})
df_loo_all = as_tibble(rbindlist(df_loo_l))
str(df_loo_all)
setwd("/data/home/shpli3/R_projects/BSS_exclude_tree_raw/fit_results/plot_ages1_35_top50_equal_interval_model_comparison")
getwd()
save(df_loo_all, file = 'df_loo_all.Rdata')

str(df_loo_all)
df_loo_all_summary = df_loo_all %>% group_by(model) %>% summarise_at(vars(elpd_diff, se_diff, se_elpd_loo), 
                                                                      c(mean, sd)) %>% arrange(
                                                                      desc(elpd_diff_fn1))
colnames(df_loo_all_summary) = c('model',
'delta.elpd', 'se.delta.elpd', 'se.elpd', 'sd.delta.elpd', 'sd.se.delta.elpd', 'sd.se.elpd')
df_loo_all_summary$delta.elpd = round(df_loo_all_summary$elpd_diff_fn1-df_loo_all_summary$elpd_diff_fn1[1], 2)
df_loo_all_summary$delta_se_elpd_diff = round(df_loo_all_summary$se_diff_fn1-df_loo_all_summary$se_diff_fn1[1], 2)
df_loo_all_summary[,7]






##### Analyse results from BH model ######
# Load packages
rm(list = ls())
library(dplyr)
#library(betareg)
library(stringr)
library(data.table)
#library(cmdstanr)
load("/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/fit_fp_top50_ages1_35_equal_interval.RData")

setwd('/data/home/shpli3/R_projects/BSS_exclude_tree_raw')

# Read data
inter_l_2 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh/parameters',
                       pattern = 'inter',
                       full.names = T)

inter_all = data.frame()
for (i in 1:length(inter_l_2)) {
 load(inter_l_2[i])
 inter_all = rbind(inter_all, inter_all_trans)
}

intra_l_2 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh/parameters',
                           pattern = 'intra',
                           full.names = T)
intra_all = data.frame()
for (i in 1:length(intra_l_2)) {
  #i = 1
 load(intra_l_2[i])
 intra_all = rbind(intra_all, intra_trans)
}

# Filter the species fitting whose rhat < 0.95 / rhat > 1.05
summary_l1 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh/summary', 
                        pattern = '^summary.', full.names = T)

#### Check the rhat
bad_fit_a = data.frame()
bad_fit_r = data.frame()
for (i in 1:length(summary_l1)) {
   #i = 1
      load(summary_l1[[i]])
      load(intra_l_2[[i]])
      sp = intra_trans$sp
      for (k in 1:length(summary_l)) {
        #k = 1
        trans_summary = summary_l[[k]]
        trans_summary_c_a = trans_summary %>% filter(grepl("^a", variable)) %>%
          filter(rhat < 0.95 | rhat > 1.05) #### 使用rhat<1.05的标准
        trans_summary_c_r = trans_summary %>% filter(grepl('^r$', variable)) %>%
          filter(rhat < 0.95 | rhat > 1.05)

        if (nrow(trans_summary_c_a) > 0 & nrow(trans_summary_c_r) > 0) {
          plot = str_split(str_split(str_split(summary_l1[[i]], '/')[[1]][5], '.rdata')[[1]][1],
                           'summary_')[[1]][2]
          species_i = sp[k]
          species_j_posi = as.numeric(gsub("[^0-9]", "", trans_summary_c_a$variable))
          species_j = sp[species_j_posi]
          bad_fit_1 = data.frame(plot = plot, species_i = species_i, species_j = species_j)
          bad_fit_a = rbind(bad_fit_1, bad_fit_a)
          
          bad_fit_2 = data.frame(plot = plot, species = species_i)
          bad_fit_r = rbind(bad_fit_2, bad_fit_r)

        } else if (nrow(trans_summary_c_a) > 0) {
          plot = str_split(str_split(str_split(summary_l1[[i]], '/')[[1]][5], '.rdata')[[1]][1],
                           'summary_')[[1]][2]
          species_i = sp[k]
          species_j_posi = as.numeric(gsub("[^0-9]", "", trans_summary_c_a$variable))
          species_j = sp[species_j_posi]
          bad_fit_1 = data.frame(plot = plot, species_i = species_i, species_j = species_j)
          bad_fit_a = rbind(bad_fit_1, bad_fit_a)
        } else if (nrow(trans_summary_c_r) > 0) {
                  species_i = sp[k]
                  bad_fit_2 = data.frame(plot = plot, species = species_i)
                  bad_fit_r = rbind(bad_fit_2, bad_fit_r)
        }  else {print('no bad fit!')}
      }
}

# Clean the bad fit
inter_all_l = split(inter_all, inter_all$f_p)
inter_all_l_c = list()
for (i in 1:length(inter_all_l)) {
  #i = 3
  trans = inter_all_l[[i]]
  f_plot = unique(trans$f_p)
  sp_c1_a = bad_fit_a %>% filter(plot == f_plot & species_i != species_j)
  if (nrow(sp_c1_a) > 0) {sp_pair_c = c(paste(sp_c1_a$species_i, sp_c1_a$species_j, sep = '.'),
  paste(sp_c1_a$species_j, sp_c1_a$species_i, sep = '.'))} else {sp_pair_c = NULL}

  sp_c1_a_intra = bad_fit_a %>% filter(plot == f_plot & species_i == species_j)
  sp_c1_r = bad_fit_r %>% filter(plot == f_plot)
  if(nrow(sp_c1_a_intra) > 0 | nrow(sp_c1_r) > 0) {sp_intra = unique(c(unique(sp_c1_a_intra$species_i),
   unique(sp_c1_r$species)))} else {sp_intra = NULL}
  
    trans_c1 = trans %>% filter(!(sp_pair %in% sp_pair_c))
    trans_c2 = trans_c1 %>% filter(!(species_i %in% sp_intra) | !(species_j %in% sp_intra) )
    inter_all_l_c[[i]] = trans_c2

}

inter_all_c = rbindlist(inter_all_l_c)
nrow(inter_all_c) / nrow(inter_all)

# good fit propotion 0.9696027

### Data fixed
colnames(inter_all_c)
inter_all_c[,7:ncol(inter_all_c)] = as.data.frame(apply(inter_all_c[,7:ncol(inter_all_c)], 2,
                                               function(x){as.numeric(x)}))
inter_all_c$lgfd = log10(inter_all_c$fd)

### Add the invasion stage
inter_all_c$stage_i = 'native'
inter_all_c$stage_j = 'native'

#sp_domin = 'f_p'
inter_all_c_l = split(inter_all_c, inter_all_c$f_p)
for (i in 1:length(inter_all_c_l)) {
  #i = 344
  trans = inter_all_c_l[[i]]
  f_p_name = unique(trans$f_p)
  i_1 = which(names(fit_fp_top50_ages1_35_equal_interval) == f_p_name)
  spname_stage = as.data.frame(fit_fp_top50_ages1_35_equal_interval[[i_1]]$stage_sp_name)  

  sp_intro = strsplit(strsplit(spname_stage$intro_sp,
                      'Plot, f_p, ')[[1]][1], ', ')[[1]]
  sp_estab = strsplit(spname_stage$estab_sp,
                       ', ')[[1]]
  sp_domin = strsplit(spname_stage$domin_sp,
                      ', ')[[1]]
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_i %in% sp_intro,]$stage_i = 'introduce'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_i %in% sp_estab,]$stage_i = 'establish'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_i %in% sp_domin,]$stage_i = 'dominant'

  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_j %in% sp_intro,]$stage_j = 'introduce'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_j %in% sp_estab,]$stage_j = 'establish'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_j %in% sp_domin,]$stage_j = 'dominant'
}

unique(inter_all_c$stage_i)
inter_all_c$stage_ij = paste(inter_all_c$stage_i, inter_all_c$stage_j, sep = '_')
inter_all_c$sp_pair = paste(inter_all_c$species_i, inter_all_c$species_j, sep = '_')
unique(inter_all_c$species_i)

## Add the different stage probability
inter_all_c$estab_i = 0
inter_all_c$domin_i = 0
unique(inter_all_c$stage_i)
inter_all_c[inter_all_c$stage_i %in% c('establish', 'dominant'), ]$estab_i = 1
inter_all_c[inter_all_c$stage_i %in% c('dominant'), ]$domin_i = 1

inter_all_c$estab_j = 0
inter_all_c$domin_j = 0
unique(inter_all_c$stage_j)
inter_all_c[inter_all_c$stage_j %in% c('establish', 'dominant'), ]$estab_j = 1
inter_all_c[inter_all_c$stage_j %in% c('dominant'), ]$domin_j = 1

## Add the gain or loss, relative abundance, and coexistence probabilities
load('code/Phylo_Func/fd_dat_all.rdata')
load('code/Phylo_Func/pd.rdata')
load('/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/sp_racover_f1_mean_fp_early_suc_1_35.rdata')
colnames(inter_all_c)
colnames(sp_racover_f1_mean_fp_early_suc_1_35)
inter_all_c$f_p_species_i = paste(inter_all_c$f_p,
                                  inter_all_c$species_i,
                                  sep = '_')
inter_all_c$f_p_species_j = paste(inter_all_c$f_p,
                                  inter_all_c$species_j,
                                  sep = '_')
inter_all_c$f_p_sp_pair = paste(inter_all_c$f_p,
                                  inter_all_c$sp_pair,
                                  sep = '_')

sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_i = sp_racover_f1_mean_fp_early_suc_1_35$f_p_species
sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_j = sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_i

colnames(intra_all)
intra_all$f_p_species = paste(intra_all$f_p, intra_all$sp, sep = '_')
intra_all$f_p_species_i = intra_all$f_p_species
intra_all$f_p_species_j = intra_all$f_p_species_i
colnames(intra_all)

inter_all_c_test = inter_all_c %>% left_join(intra_all[,c(6, 8)], 
            by = 'f_p_species_i') %>%
  left_join(intra_all[,c(6, 9)], 
            by = 'f_p_species_j') %>% 
  left_join(sp_racover_f1_mean_fp_early_suc_1_35[,c(4, 5)], 
            by = 'f_p_species_i') %>%
  left_join(sp_racover_f1_mean_fp_early_suc_1_35[,c(4, 6)], 
            by = 'f_p_species_j') %>%    
  left_join(fd_dat_all[,c(3:14)], 
            by = 'sp_pair') %>% 
  left_join(pd_dat[,c(3, 4)], 
            by = 'sp_pair')

colnames(inter_all_c_test)
colname_1 = gsub('\\.x', '_i', colnames(inter_all_c_test))          
colname_2 = gsub('\\.y', '_j', colname_1)
colnames(inter_all_c_test) = colname_2
inter_all_c_alltime = inter_all_c_test
colnames(inter_all_c_alltime)
summary(inter_all_c_alltime)
save(inter_all_c_alltime, file = 'fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh/inter_all_c_alltime_rhat_105.rdata')
getwd()




##### Analyse results from BH_partialb model ######
# Load packages
rm(list = ls())
library(dplyr)
#library(betareg)
library(stringr)
library(data.table)
#library(cmdstanr)
load("/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/fit_fp_top50_ages1_35_equal_interval.RData")

setwd('/data/home/shpli3/R_projects/BSS_exclude_tree_raw')

# Read data
inter_l_2 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_partialb/parameters',
                       pattern = 'inter',
                       full.names = T)

inter_all = data.frame()
for (i in 1:length(inter_l_2)) {
 load(inter_l_2[i])
 inter_all = rbind(inter_all, inter_all_trans)
}

intra_l_2 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_partialb/parameters',
                           pattern = 'intra',
                           full.names = T)
intra_all = data.frame()
for (i in 1:length(intra_l_2)) {
  #i = 1
 load(intra_l_2[i])
 intra_all = rbind(intra_all, intra_trans)
}

# Filter the species fitting whose rhat < 0.95 / rhat > 1.05
summary_l1 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_partialb/summary', 
                        pattern = '^summary.', full.names = T)

#### Check the rhat
bad_fit_a = data.frame()
bad_fit_r = data.frame()
for (i in 1:length(summary_l1)) {
   #i = 1
      load(summary_l1[[i]])
      load(intra_l_2[[i]])
      sp = intra_trans$sp
      for (k in 1:length(summary_l)) {
        #k = 1
        trans_summary = summary_l[[k]]
        trans_summary_c_a = trans_summary %>% filter(grepl("^a", variable)) %>%
          filter(rhat < 0.95 | rhat > 1.05) #### 使用rhat<1.05的标准
        trans_summary_c_r = trans_summary %>% filter(grepl('^r$', variable)) %>%
          filter(rhat < 0.95 | rhat > 1.05)

        if (nrow(trans_summary_c_a) > 0 & nrow(trans_summary_c_r) > 0) {
          plot = str_split(str_split(str_split(summary_l1[[i]], '/')[[1]][5], '.rdata')[[1]][1],
                           'summary_')[[1]][2]
          species_i = sp[k]
          species_j_posi = as.numeric(gsub("[^0-9]", "", trans_summary_c_a$variable))
          species_j = sp[species_j_posi]
          bad_fit_1 = data.frame(plot = plot, species_i = species_i, species_j = species_j)
          bad_fit_a = rbind(bad_fit_1, bad_fit_a)
          
          bad_fit_2 = data.frame(plot = plot, species = species_i)
          bad_fit_r = rbind(bad_fit_2, bad_fit_r)

        } else if (nrow(trans_summary_c_a) > 0) {
          plot = str_split(str_split(str_split(summary_l1[[i]], '/')[[1]][5], '.rdata')[[1]][1],
                           'summary_')[[1]][2]
          species_i = sp[k]
          species_j_posi = as.numeric(gsub("[^0-9]", "", trans_summary_c_a$variable))
          species_j = sp[species_j_posi]
          bad_fit_1 = data.frame(plot = plot, species_i = species_i, species_j = species_j)
          bad_fit_a = rbind(bad_fit_1, bad_fit_a)
        } else if (nrow(trans_summary_c_r) > 0) {
                  species_i = sp[k]
                  bad_fit_2 = data.frame(plot = plot, species = species_i)
                  bad_fit_r = rbind(bad_fit_2, bad_fit_r)
        }  else {print('no bad fit!')}
      }
}

# Clean the bad fit
inter_all_l = split(inter_all, inter_all$f_p)
inter_all_l_c = list()
for (i in 1:length(inter_all_l)) {
  #i = 3
  trans = inter_all_l[[i]]
  f_plot = unique(trans$f_p)
  sp_c1_a = bad_fit_a %>% filter(plot == f_plot & species_i != species_j)
  if (nrow(sp_c1_a) > 0) {sp_pair_c = c(paste(sp_c1_a$species_i, sp_c1_a$species_j, sep = '.'),
  paste(sp_c1_a$species_j, sp_c1_a$species_i, sep = '.'))} else {sp_pair_c = NULL}

  sp_c1_a_intra = bad_fit_a %>% filter(plot == f_plot & species_i == species_j)
  sp_c1_r = bad_fit_r %>% filter(plot == f_plot)
  if(nrow(sp_c1_a_intra) > 0 | nrow(sp_c1_r) > 0) {sp_intra = unique(c(unique(sp_c1_a_intra$species_i),
   unique(sp_c1_r$species)))} else {sp_intra = NULL}
  
    trans_c1 = trans %>% filter(!(sp_pair %in% sp_pair_c))
    trans_c2 = trans_c1 %>% filter(!(species_i %in% sp_intra) | !(species_j %in% sp_intra) )
    inter_all_l_c[[i]] = trans_c2

}

inter_all_c = rbindlist(inter_all_l_c)
nrow(inter_all_c) / nrow(inter_all)

# good fit propotion 0.9599544

### Data fixed
colnames(inter_all_c)
inter_all_c[,7:ncol(inter_all_c)] = as.data.frame(apply(inter_all_c[,7:ncol(inter_all_c)], 2,
                                               function(x){as.numeric(x)}))
inter_all_c$lgfd = log10(inter_all_c$fd)

### Add the invasion stage
inter_all_c$stage_i = 'native'
inter_all_c$stage_j = 'native'

#sp_domin = 'f_p'
inter_all_c_l = split(inter_all_c, inter_all_c$f_p)
for (i in 1:length(inter_all_c_l)) {
  #i = 344
  trans = inter_all_c_l[[i]]
  f_p_name = unique(trans$f_p)
  i_1 = which(names(fit_fp_top50_ages1_35_equal_interval) == f_p_name)
  spname_stage = as.data.frame(fit_fp_top50_ages1_35_equal_interval[[i_1]]$stage_sp_name)  

  sp_intro = strsplit(strsplit(spname_stage$intro_sp,
                      'Plot, f_p, ')[[1]][1], ', ')[[1]]
  sp_estab = strsplit(spname_stage$estab_sp,
                       ', ')[[1]]
  sp_domin = strsplit(spname_stage$domin_sp,
                      ', ')[[1]]
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_i %in% sp_intro,]$stage_i = 'introduce'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_i %in% sp_estab,]$stage_i = 'establish'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_i %in% sp_domin,]$stage_i = 'dominant'

  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_j %in% sp_intro,]$stage_j = 'introduce'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_j %in% sp_estab,]$stage_j = 'establish'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_j %in% sp_domin,]$stage_j = 'dominant'
}

unique(inter_all_c$stage_i)
inter_all_c$stage_ij = paste(inter_all_c$stage_i, inter_all_c$stage_j, sep = '_')
inter_all_c$sp_pair = paste(inter_all_c$species_i, inter_all_c$species_j, sep = '_')
unique(inter_all_c$species_i)

## Add the different stage probability
inter_all_c$estab_i = 0
inter_all_c$domin_i = 0
unique(inter_all_c$stage_i)
inter_all_c[inter_all_c$stage_i %in% c('establish', 'dominant'), ]$estab_i = 1
inter_all_c[inter_all_c$stage_i %in% c('dominant'), ]$domin_i = 1

inter_all_c$estab_j = 0
inter_all_c$domin_j = 0
unique(inter_all_c$stage_j)
inter_all_c[inter_all_c$stage_j %in% c('establish', 'dominant'), ]$estab_j = 1
inter_all_c[inter_all_c$stage_j %in% c('dominant'), ]$domin_j = 1

## Add the gain or loss, relative abundance, and coexistence probabilities
load('code/Phylo_Func/fd_dat_all.rdata')
load('code/Phylo_Func/pd.rdata')
load('/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/sp_racover_f1_mean_fp_early_suc_1_35.rdata')
colnames(inter_all_c)
colnames(sp_racover_f1_mean_fp_early_suc_1_35)
inter_all_c$f_p_species_i = paste(inter_all_c$f_p,
                                  inter_all_c$species_i,
                                  sep = '_')
inter_all_c$f_p_species_j = paste(inter_all_c$f_p,
                                  inter_all_c$species_j,
                                  sep = '_')
inter_all_c$f_p_sp_pair = paste(inter_all_c$f_p,
                                  inter_all_c$sp_pair,
                                  sep = '_')

sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_i = sp_racover_f1_mean_fp_early_suc_1_35$f_p_species
sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_j = sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_i

colnames(intra_all)
intra_all$f_p_species = paste(intra_all$f_p, intra_all$sp, sep = '_')
intra_all$f_p_species_i = intra_all$f_p_species
intra_all$f_p_species_j = intra_all$f_p_species_i
colnames(intra_all)

inter_all_c_test = inter_all_c %>% left_join(intra_all[,c("r_trans_sp",
 "b_trans_sp", "f_p_species_i")], 
            by = 'f_p_species_i') %>%
  left_join(intra_all[,c("r_trans_sp",
 "b_trans_sp", "f_p_species_j")], 
            by = 'f_p_species_j') %>% 
  left_join(sp_racover_f1_mean_fp_early_suc_1_35[,c(4, 5)], 
            by = 'f_p_species_i') %>%
  left_join(sp_racover_f1_mean_fp_early_suc_1_35[,c(4, 6)], 
            by = 'f_p_species_j') %>%    
  left_join(fd_dat_all[,c(3:14)], 
            by = 'sp_pair') %>% 
  left_join(pd_dat[,c(3, 4)], 
            by = 'sp_pair')

colnames(inter_all_c_test)
colname_1 = gsub('\\.x', '_i', colnames(inter_all_c_test))          
colname_2 = gsub('\\.y', '_j', colname_1)
colnames(inter_all_c_test) = colname_2
inter_all_c_alltime = inter_all_c_test
colnames(inter_all_c_alltime)
summary(inter_all_c_alltime)
save(inter_all_c_alltime, file = 'fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_partialb/inter_all_c_alltime_rhat_105.rdata')
getwd()



##### Analyse results from BH_allb model ######
# Load packages
rm(list = ls())
library(dplyr)
#library(betareg)
library(stringr)
library(data.table)
#library(cmdstanr)
load("/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/fit_fp_top50_ages1_35_equal_interval.RData")

setwd('/data/home/shpli3/R_projects/BSS_exclude_tree_raw')

# Read data
inter_l_2 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb/parameters',
                       pattern = 'inter',
                       full.names = T)

inter_all = data.frame()
for (i in 1:length(inter_l_2)) {
 load(inter_l_2[i])
 inter_all = rbind(inter_all, inter_all_trans)
}

intra_l_2 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb/parameters',
                           pattern = 'intra',
                           full.names = T)
intra_all = data.frame()
for (i in 1:length(intra_l_2)) {
  #i = 1
 load(intra_l_2[i])
 intra_all = rbind(intra_all, intra_trans)
}

# Filter the species fitting whose rhat < 0.95 / rhat > 1.05
summary_l1 = list.files('fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb/summary', 
                        pattern = '^summary.', full.names = T)

#### Check the rhat
bad_fit_a = data.frame()
bad_fit_r = data.frame()
for (i in 1:length(summary_l1)) {
   #i = 1
      load(summary_l1[[i]])
      load(intra_l_2[[i]])
      sp = intra_trans$sp
      for (k in 1:length(summary_l)) {
        #k = 1
        trans_summary = summary_l[[k]]
        trans_summary_c_a = trans_summary %>% filter(grepl("^a", variable)) %>%
          filter(rhat < 0.95 | rhat > 1.05) #### 使用rhat<1.05的标准
        trans_summary_c_r = trans_summary %>% filter(grepl('^r$', variable)) %>%
          filter(rhat < 0.95 | rhat > 1.05)

        if (nrow(trans_summary_c_a) > 0 & nrow(trans_summary_c_r) > 0) {
          plot = str_split(str_split(str_split(summary_l1[[i]], '/')[[1]][5], '.rdata')[[1]][1],
                           'summary_')[[1]][2]
          species_i = sp[k]
          species_j_posi = as.numeric(gsub("[^0-9]", "", trans_summary_c_a$variable))
          species_j = sp[species_j_posi]
          bad_fit_1 = data.frame(plot = plot, species_i = species_i, species_j = species_j)
          bad_fit_a = rbind(bad_fit_1, bad_fit_a)
          
          bad_fit_2 = data.frame(plot = plot, species = species_i)
          bad_fit_r = rbind(bad_fit_2, bad_fit_r)

        } else if (nrow(trans_summary_c_a) > 0) {
          plot = str_split(str_split(str_split(summary_l1[[i]], '/')[[1]][5], '.rdata')[[1]][1],
                           'summary_')[[1]][2]
          species_i = sp[k]
          species_j_posi = as.numeric(gsub("[^0-9]", "", trans_summary_c_a$variable))
          species_j = sp[species_j_posi]
          bad_fit_1 = data.frame(plot = plot, species_i = species_i, species_j = species_j)
          bad_fit_a = rbind(bad_fit_1, bad_fit_a)
        } else if (nrow(trans_summary_c_r) > 0) {
                  species_i = sp[k]
                  bad_fit_2 = data.frame(plot = plot, species = species_i)
                  bad_fit_r = rbind(bad_fit_2, bad_fit_r)
        }  else {print('no bad fit!')}
      }
}

# Clean the bad fit
inter_all_l = split(inter_all, inter_all$f_p)
inter_all_l_c = list()
for (i in 1:length(inter_all_l)) {
  #i = 3
  trans = inter_all_l[[i]]
  f_plot = unique(trans$f_p)
  sp_c1_a = bad_fit_a %>% filter(plot == f_plot & species_i != species_j)
  if (nrow(sp_c1_a) > 0) {sp_pair_c = c(paste(sp_c1_a$species_i, sp_c1_a$species_j, sep = '.'),
  paste(sp_c1_a$species_j, sp_c1_a$species_i, sep = '.'))} else {sp_pair_c = NULL}

  sp_c1_a_intra = bad_fit_a %>% filter(plot == f_plot & species_i == species_j)
  sp_c1_r = bad_fit_r %>% filter(plot == f_plot)
  if(nrow(sp_c1_a_intra) > 0 | nrow(sp_c1_r) > 0) {sp_intra = unique(c(unique(sp_c1_a_intra$species_i),
   unique(sp_c1_r$species)))} else {sp_intra = NULL}
  
    trans_c1 = trans %>% filter(!(sp_pair %in% sp_pair_c))
    trans_c2 = trans_c1 %>% filter(!(species_i %in% sp_intra) | !(species_j %in% sp_intra) )
    inter_all_l_c[[i]] = trans_c2

}

inter_all_c = rbindlist(inter_all_l_c)
nrow(inter_all_c) / nrow(inter_all)

# good fit propotion 0.9603334

### Data fixed
colnames(inter_all_c)
inter_all_c[,7:ncol(inter_all_c)] = as.data.frame(apply(inter_all_c[,7:ncol(inter_all_c)], 2,
                                               function(x){as.numeric(x)}))
inter_all_c$lgfd = log10(inter_all_c$fd)

### Add the invasion stage
inter_all_c$stage_i = 'native'
inter_all_c$stage_j = 'native'

#sp_domin = 'f_p'
inter_all_c_l = split(inter_all_c, inter_all_c$f_p)
for (i in 1:length(inter_all_c_l)) {
  #i = 344
  trans = inter_all_c_l[[i]]
  f_p_name = unique(trans$f_p)
  i_1 = which(names(fit_fp_top50_ages1_35_equal_interval) == f_p_name)
  spname_stage = as.data.frame(fit_fp_top50_ages1_35_equal_interval[[i_1]]$stage_sp_name)  

  sp_intro = strsplit(strsplit(spname_stage$intro_sp,
                      'Plot, f_p, ')[[1]][1], ', ')[[1]]
  sp_estab = strsplit(spname_stage$estab_sp,
                       ', ')[[1]]
  sp_domin = strsplit(spname_stage$domin_sp,
                      ', ')[[1]]
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_i %in% sp_intro,]$stage_i = 'introduce'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_i %in% sp_estab,]$stage_i = 'establish'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_i %in% sp_domin,]$stage_i = 'dominant'

  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_j %in% sp_intro,]$stage_j = 'introduce'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_j %in% sp_estab,]$stage_j = 'establish'
  inter_all_c[inter_all_c$f_p == f_p_name & inter_all_c$species_j %in% sp_domin,]$stage_j = 'dominant'
}

unique(inter_all_c$stage_i)
inter_all_c$stage_ij = paste(inter_all_c$stage_i, inter_all_c$stage_j, sep = '_')
inter_all_c$sp_pair = paste(inter_all_c$species_i, inter_all_c$species_j, sep = '_')
unique(inter_all_c$species_i)

## Add the different stage probability
inter_all_c$estab_i = 0
inter_all_c$domin_i = 0
unique(inter_all_c$stage_i)
inter_all_c[inter_all_c$stage_i %in% c('establish', 'dominant'), ]$estab_i = 1
inter_all_c[inter_all_c$stage_i %in% c('dominant'), ]$domin_i = 1

inter_all_c$estab_j = 0
inter_all_c$domin_j = 0
unique(inter_all_c$stage_j)
inter_all_c[inter_all_c$stage_j %in% c('establish', 'dominant'), ]$estab_j = 1
inter_all_c[inter_all_c$stage_j %in% c('dominant'), ]$domin_j = 1

## Add the gain or loss, relative abundance, and coexistence probabilities
load('code/Phylo_Func/fd_dat_all.rdata')
load('code/Phylo_Func/pd.rdata')
load('/data/home/shpli3/my_pc/BSS/code/data preparation/transformed data/sp_racover_f1_mean_fp_early_suc_1_35.rdata')
colnames(inter_all_c)
colnames(sp_racover_f1_mean_fp_early_suc_1_35)
inter_all_c$f_p_species_i = paste(inter_all_c$f_p,
                                  inter_all_c$species_i,
                                  sep = '_')
inter_all_c$f_p_species_j = paste(inter_all_c$f_p,
                                  inter_all_c$species_j,
                                  sep = '_')
inter_all_c$f_p_sp_pair = paste(inter_all_c$f_p,
                                  inter_all_c$sp_pair,
                                  sep = '_')

sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_i = sp_racover_f1_mean_fp_early_suc_1_35$f_p_species
sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_j = sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_i

colnames(intra_all)
intra_all$f_p_species = paste(intra_all$f_p, intra_all$sp, sep = '_')
intra_all$f_p_species_i = intra_all$f_p_species
intra_all$f_p_species_j = intra_all$f_p_species_i
colnames(intra_all)

inter_all_c_test = inter_all_c %>% left_join(intra_all[,c("r_trans_sp",
 "b_trans_sp", "f_p_species_i")], 
            by = 'f_p_species_i') %>%
  left_join(intra_all[,c("r_trans_sp",
 "b_trans_sp", "f_p_species_j")], 
            by = 'f_p_species_j') %>% 
  left_join(sp_racover_f1_mean_fp_early_suc_1_35[,c(4, 5)], 
            by = 'f_p_species_i') %>%
  left_join(sp_racover_f1_mean_fp_early_suc_1_35[,c(4, 6)], 
            by = 'f_p_species_j') %>%    
  left_join(fd_dat_all[,c(3:14)], 
            by = 'sp_pair') %>% 
  left_join(pd_dat[,c(3, 4)], 
            by = 'sp_pair')

colnames(inter_all_c_test)
colname_1 = gsub('\\.x', '_i', colnames(inter_all_c_test))          
colname_2 = gsub('\\.y', '_j', colname_1)
colnames(inter_all_c_test) = colname_2
inter_all_c_alltime = inter_all_c_test
colnames(inter_all_c_alltime)
summary(inter_all_c_alltime)
save(inter_all_c_alltime, file = 'fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_allb/inter_all_c_alltime_rhat_105.rdata')
getwd()
