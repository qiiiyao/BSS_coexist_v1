# Load packages
rm(list = ls())

require(stringr)
require(dplyr)
require(data.table)
setwd("~/BSS_coexist_v1")

load("code/data preparation/transformed data/fit_fp_top50_ages1_35_equal_interval.RData")
load("code/data preparation/transformed data/re_cover_ab_f1_ages1_35_fp.rdata")
dat_all_l = lapply(1:length(fit_fp_top50_ages1_35_equal_interval), function(x){
  #x = 1
  stage_sp = fit_fp_top50_ages1_35_equal_interval[[x]]$stage_sp_name
  re_cover_ab = re_cover_ab_f1_ages1_35_fp[[x]]
  alien_sp = c(str_split_1(stage_sp[[2]], ', '), 
               str_split_1(stage_sp[[3]], ', '))
  native_sp = setdiff(colnames(re_cover_ab)[10:ncol(re_cover_ab)],alien_sp)
  f_p_1 = stage_sp[[1]]
  dat = data.frame(f_p = f_p_1, 
                   species_i = rep(alien_sp, each = length(native_sp)),
                   species_j = rep(native_sp, length(alien_sp)))
  return(dat)
})

dat_all = as.data.frame(rbindlist(dat_all_l))
dat_all$sp_pair = paste(dat_all$species_i, dat_all$species_j, sep = '_')

## Add the gain or loss and relative abundance
load('code/Phylo_Func/fd_dat_all.rdata')
load('code/Phylo_Func/pd.rdata')
load("code/data preparation/transformed data/sp_racover_f1_mean_fp_early_suc_1_35.rdata")

colnames(sp_racover_f1_mean_fp_early_suc_1_35)
dat_all$f_p_species_i = paste(dat_all$f_p,
                              dat_all$species_i,
                              sep = '_')
dat_all$f_p_species_j = paste(dat_all$f_p,
                              dat_all$species_j,
                              sep = '_')

sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_i = sp_racover_f1_mean_fp_early_suc_1_35$f_p_species
sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_j = sp_racover_f1_mean_fp_early_suc_1_35$f_p_species_i
colnames(sp_racover_f1_mean_fp_early_suc_1_35)
dat_all_test = dat_all %>% left_join(sp_racover_f1_mean_fp_early_suc_1_35[,c(4, 5)], 
                                     by = 'f_p_species_i') %>%
  left_join(sp_racover_f1_mean_fp_early_suc_1_35[,c(4, 6)], 
            by = 'f_p_species_j') %>% 
  left_join(fd_dat_all[,c(3:14)], 
            by = 'sp_pair') %>% 
  left_join(pd_dat[,c(3, 4)], 
            by = 'sp_pair')

colname_1 = gsub('\\.x', '_i', colnames(dat_all_test))          
colname_2 = gsub('\\.y', '_j', colname_1)
colnames(dat_all_test) = colname_2
dat_all_alltime = dat_all_test
colnames(dat_all_alltime)
summary(dat_all_alltime)
save(dat_all_alltime, 
     file = 'results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_equal_interval_model_comparison/dat_all_alltime.rdata')
rm(list = ls()) ### free the memory

############### Fast Start ####################
rm(list = ls())
setwd("~/BSS_coexist_v1")
load('results/fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_partialb/inter_all_c_alltime.rdata')
load('results/fit_results/plot_ages1_35_top50_equal_interval_model_comparison/dat_all_alltime.rdata')

### load required packages
library(dplyr)
library(betareg)
library(stringr)
library(data.table)
library(ape)
library(colRoz)
library(lme4)
library(lmerTest)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)
library(ggeffects)
library(stringr)
library(viridis)
library(showtext)
library(extrafont)
windowsFonts(Arial=windowsFont("Arial"))
source('code/function/plot_func.R')

inter_all_c = inter_all_c_alltime
rm(inter_all_c_alltime)
colnames(inter_all_c)

gre = "#55a868ff"
ora = "#dd8452ff"
yel = "#ccb974ff"
blu = "#4c72b0ff"
colors_4d = c(blu,ora,gre,yel)
colors_2d = c(ora, blu)

# get the abs value of lgfd
inter_all_c$lgfd= log10(inter_all_c$fd)

#### Check the correlation between nd and fd
#cor.test(inter_all_c$nd, inter_all_c$ablgfd) 
## Full data but corr coeff just 0.13
## correlation
#plot(inter_all_c$nd, inter_all_c$ablgfd)

## Add stage introduce
inter_all_c$intro_i = 1
inter_all_c$intro_j = 1
inter_all_c[inter_all_c$stage_i == 'native',]$intro_i = 0
inter_all_c[inter_all_c$stage_j == 'native',]$intro_j = 0

inter_all_c$stage_ij_estab = inter_all_c$stage_ij
inter_all_c$stage_ij_domin = inter_all_c$stage_ij

inter_all_c$stage_i_estab = 'native'
inter_all_c$stage_j_estab = 'native'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 0,]$stage_i_estab = 'intro.noestab'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 1,]$stage_i_estab = 'intro.estab'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 0,]$stage_j_estab = 'intro.noestab'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 1,]$stage_j_estab = 'intro.estab'
inter_all_c$stage_ij_estab = paste(inter_all_c$stage_i_estab,
                                   inter_all_c$stage_j_estab,
                                   sep = '_')

inter_all_c$stage_i_domin = 'native'
inter_all_c$stage_j_domin = 'native'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 0&
              inter_all_c$domin_i == 0,]$stage_i_domin = 'intro.nodomin'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 1&
              inter_all_c$domin_i == 0,]$stage_i_domin = 'estab.nodomin'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 1&
              inter_all_c$domin_i == 1,]$stage_i_domin = 'estab.domin'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 0&
              inter_all_c$domin_j == 0,]$stage_j_domin = 'intro.nodomin'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 1&
              inter_all_c$domin_j == 0,]$stage_j_domin = 'estab.nodomin'
inter_all_c[inter_all_c$intro_j == 1&
              inter_all_c$estab_j == 1&
              inter_all_c$domin_j == 1,]$stage_j_domin = 'estab.domin'
inter_all_c$stage_ij_domin = paste(inter_all_c$stage_i_domin,
                                   inter_all_c$stage_j_domin,
                                   sep = '_')

inter_all_c = inter_all_c %>% relocate(intro_i, .after = stage_ij) %>% 
  relocate(intro_j, .after = intro_i) %>% 
  relocate(stage_ij_estab, .after = stage_ij) %>% 
  relocate(stage_i_estab, .after = stage_ij) %>% 
  relocate(stage_j_estab, .after = stage_i_estab) %>% 
  relocate(stage_ij_domin, .after = stage_ij_estab) %>% 
  relocate(stage_i_domin, .after = stage_ij_estab) %>% 
  relocate(stage_j_domin, .after = stage_i_domin)

#### Analyze invasion success ####
inter_all_c$ra_m_real_t_i = as.numeric(inter_all_c$ra_m_real_t_i)
inter_all_c$ra_m_real_t_j = as.numeric(inter_all_c$ra_m_real_t_j)
inter_all_c_l_2 = split(inter_all_c, inter_all_c$f_p)
dat_all_2 = split(dat_all_alltime, dat_all_alltime$f_p)

##### Species level #####
dat_suc_sp = data.frame()

for (i in 1:length(inter_all_c_l_2)) {
  #i = 1
  trans_plot = inter_all_c_l_2[[i]]
  inv_suc = trans_plot %>% filter(stage_i != 'native' & stage_j == 'native')
  inv_sp = unique(inv_suc$species_i)
  inv_suc_all = dat_all_2[[i]]
  
  for (j in 1:length(inv_sp)) {
    #j = 4
    inv_suc_sp = inv_suc %>% filter(species_i == inv_sp[j])
    inv_suc_sp_all = inv_suc_all %>% filter(species_i == inv_sp[j])
    demo_rate = unique(inv_suc_sp$r_trans_sp_i)
    mnd = mean(inv_suc_sp$nd)
    mlgfd = mean(inv_suc_sp$lgfd)
    mablgfd = mean(inv_suc_sp$ablgfd)
    mdemo_ratio = mean(inv_suc_sp$demo_ratio)
    mcomp_ratio = mean(inv_suc_sp$comp_ratio)
    mpd = mean(inv_suc_sp$Phylo_dis)
    mfunc_d = mean(inv_suc_sp$Multi_traits)
    mconti_func_d = mean(inv_suc_sp$Multi_conti_traits)
    mgrowth = mean(inv_suc_sp$growth)
    mspan = mean(inv_suc_sp$span)
    mpollination = mean(inv_suc_sp$pollination)
    mdispersal = mean(inv_suc_sp$dispersal)
    mclonality = mean(inv_suc_sp$clonality)
    mheight = mean(inv_suc_sp$height)
    mldmc = mean(inv_suc_sp$ldmc)
    msla = mean(inv_suc_sp$sla)
    mseedmass = mean(inv_suc_sp$seedmass)
    
    mpd_all = mean(inv_suc_sp_all$Phylo_dis)
    mfunc_d_all = mean(inv_suc_sp_all$Multi_traits)
    mconti_func_d_all = mean(inv_suc_sp_all$Multi_conti_traits)
    mgrowth_all = mean(inv_suc_sp_all$growth)
    mspan_all = mean(inv_suc_sp_all$span)
    mpollination_all = mean(inv_suc_sp_all$pollination)
    mdispersal_all = mean(inv_suc_sp_all$dispersal)
    mclonality_all = mean(inv_suc_sp_all$clonality)
    mheight_all = mean(inv_suc_sp_all$height)
    mldmc_all = mean(inv_suc_sp_all$ldmc)
    msla_all = mean(inv_suc_sp_all$sla)
    mseedmass_all = mean(inv_suc_sp_all$seedmass)
    
    mnd.a = sum(inv_suc_sp$nd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mlgfd.a = sum(inv_suc_sp$lgfd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mablgfd.a = sum(inv_suc_sp$ablgfd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mdemo_ratio.a = sum(inv_suc_sp$demo_ratio*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mcomp_ratio.a = sum(inv_suc_sp$comp_ratio*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mpd.a = sum(inv_suc_sp$Phylo_dis*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mfunc_d.a = sum(inv_suc_sp$Multi_traits*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mconti_func_d.a = sum(inv_suc_sp$Multi_conti_traits*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mgrowth.a = sum(inv_suc_sp$growth*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mspan.a = sum(inv_suc_sp$span*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mpollination.a = sum(inv_suc_sp$pollination*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mdispersal.a = sum(inv_suc_sp$dispersal*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mclonality.a = sum(inv_suc_sp$clonality*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mheight.a = sum(inv_suc_sp$height*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mldmc.a = sum(inv_suc_sp$ldmc*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    msla.a = sum(inv_suc_sp$sla*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mseedmass.a = sum(inv_suc_sp$seedmass*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/
      sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    
    mpd.a_all = sum(inv_suc_sp_all$Phylo_dis*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mfunc_d.a_all = sum(inv_suc_sp_all$Multi_traits*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mconti_func_d.a_all = sum(inv_suc_sp_all$Multi_conti_traits*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mgrowth.a_all = sum(inv_suc_sp_all$growth*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mspan.a_all = sum(inv_suc_sp_all$span*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mpollination.a_all = sum(inv_suc_sp_all$pollination*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mdispersal.a_all = sum(inv_suc_sp_all$dispersal*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mclonality.a_all = sum(inv_suc_sp_all$clonality*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mheight.a_all = sum(inv_suc_sp_all$height*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mldmc.a_all = sum(inv_suc_sp_all$ldmc*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    msla.a_all = sum(inv_suc_sp_all$sla*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mseedmass.a_all = sum(inv_suc_sp_all$seedmass*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/
      sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    
    mnnd = min(inv_suc_sp$nd, na.rm = T)
    mnlgfd = min(inv_suc_sp$lgfd, na.rm = T)
    mnablgfd = min(inv_suc_sp$ablgfd, na.rm = T)
    mndemo_ratio = min(inv_suc_sp$demo_ratio, na.rm = T)
    mncomp_ratio = min(inv_suc_sp$comp_ratio, na.rm = T)
    mntd = min(inv_suc_sp$Phylo_dis, na.rm = T)
    mnfunc_d = min(inv_suc_sp$Multi_traits, na.rm = T)
    mnconti_func_d = min(inv_suc_sp$Multi_conti_traits, na.rm = T)
    mngrowth = min(inv_suc_sp$growth, na.rm = T)
    mnspan = min(inv_suc_sp$span, na.rm = T)
    mnpollination = min(inv_suc_sp$pollination, na.rm = T)
    mndispersal = min(inv_suc_sp$dispersal, na.rm = T)
    mnclonality = min(inv_suc_sp$clonality, na.rm = T)
    mnheight = min(inv_suc_sp$height, na.rm = T)
    mnldmc = min(inv_suc_sp$ldmc, na.rm = T)
    mnsla = min(inv_suc_sp$sla, na.rm = T)
    mnseedmass = min(inv_suc_sp$seedmass, na.rm = T)
    
    mntd_all = min(inv_suc_sp_all$Phylo_dis, na.rm = T)
    mnfunc_d_all = min(inv_suc_sp_all$Multi_traits, na.rm = T)
    mnconti_func_d_all = min(inv_suc_sp_all$Multi_conti_traits, na.rm = T)
    mngrowth_all = min(inv_suc_sp_all$growth, na.rm = T)
    mnspan_all = min(inv_suc_sp_all$span, na.rm = T)
    mnpollination_all = min(inv_suc_sp_all$pollination, na.rm = T)
    mndispersal_all = min(inv_suc_sp_all$dispersal, na.rm = T)
    mnclonality_all = min(inv_suc_sp_all$clonality, na.rm = T)
    mnheight_all = min(inv_suc_sp_all$height, na.rm = T)
    mnldmc_all = min(inv_suc_sp_all$ldmc, na.rm = T)
    mnsla_all = min(inv_suc_sp_all$sla, na.rm = T)
    mnseedmass_all = min(inv_suc_sp_all$seedmass, na.rm = T)
    
    
    dat_suc_sp_1 = data.frame(f_p = unique(trans_plot$f_p), plot = unique(trans_plot$f_p), field = unique(trans_plot$field),
                              species = inv_sp[j], 
                              stage = unique(inv_suc_sp$stage_i),
                              estab = unique(inv_suc_sp$estab_i),
                              domin = unique(inv_suc_sp$domin_i),
                              demo_rate = demo_rate,
                              mnd = mnd, mlgfd = mlgfd, mablgfd = mablgfd,
                              mdemo_ratio = mdemo_ratio,
                              mcomp_ratio = mcomp_ratio,
                              mpd = mpd, mfunc_d = mfunc_d,
                              mconti_func_d = mconti_func_d, mgrowth = mgrowth, mspan = mspan, 
                              mpollination = mpollination, mdispersal = mdispersal, 
                              mclonality = mclonality, mheight = mheight, mldmc = mldmc,
                              msla = msla, mseedmass = mseedmass, 
                              
                              mpd_all = mpd_all, mfunc_d_all = mfunc_d_all,
                              mconti_func_d_all = mconti_func_d_all,
                              mgrowth_all = mgrowth_all, mspan_all = mspan_all, 
                              mpollination_all = mpollination_all,
                              mdispersal_all = mdispersal_all, 
                              mclonality_all = mclonality_all,
                              mheight_all = mheight_all, mldmc_all = mldmc_all,
                              msla_all = msla_all, mseedmass_all = mseedmass_all,
                              
                              mnd.a = mnd.a, mlgfd.a = mlgfd.a, mablgfd.a = mablgfd.a,
                              mdemo_ratio.a = mdemo_ratio.a,
                              mcomp_ratio.a = mcomp_ratio.a,
                              mpd.a = mpd.a,
                              mfunc_d.a = mfunc_d.a, mconti_func_d.a = mconti_func_d.a,
                              mgrowth.a = mgrowth.a, mspan.a = mspan.a,mpollination.a = mpollination.a,
                              mdispersal.a = mdispersal.a, mclonality.a = mclonality.a,
                              mheight.a = mheight.a, mldmc.a = mldmc.a, msla.a = msla.a,
                              mseedmass.a = mseedmass.a, 
                              
                              mpd.a_all = mpd.a_all,
                              mfunc_d.a_all = mfunc_d.a_all,
                              mconti_func_d.a_all = mconti_func_d.a_all,
                              mgrowth.a_all = mgrowth.a_all, mspan.a_all = mspan.a_all,
                              mpollination.a_all = mpollination.a_all,
                              mdispersal.a_all = mdispersal.a_all, mclonality.a_all = mclonality.a_all,
                              mheight.a_all = mheight.a_all, mldmc.a_all = mldmc.a_all,
                              msla.a_all = msla.a_all, mseedmass.a_all = mseedmass.a_all,
                              
                              mnnd = mnnd, mnlgfd = mnlgfd, mnablgfd = mnablgfd,
                              mndemo_ratio = mndemo_ratio,
                              mncomp_ratio = mncomp_ratio,
                              mntd = mntd,
                              mnfunc_d = mnfunc_d, mnconti_func_d = mnconti_func_d, 
                              mngrowth = mngrowth, mnspan = mnspan, 
                              mnpollination = mnpollination, mndispersal = mndispersal, 
                              mnclonality = mnclonality, mnheight = mnheight, mnldmc = mnldmc,
                              mnsla = mnsla, mnseedmass = mnseedmass,
                              
                              mntd_all = mntd_all,
                              mnfunc_d_all = mnfunc_d_all, 
                              mnconti_func_d_all = mnconti_func_d_all, 
                              mngrowth_all = mngrowth_all, mnspan_all = mnspan_all, 
                              mnpollination_all = mnpollination_all,
                              mndispersal_all = mndispersal_all, 
                              mnclonality_all = mnclonality_all, 
                              mnheight_all = mnheight_all, mnldmc_all = mnldmc_all,
                              mnsla_all = mnsla_all, mnseedmass_all = mnseedmass_all)
    dat_suc_sp = rbind(dat_suc_sp, dat_suc_sp_1)
  }
}

summary(dat_suc_sp)
dat_suc_sp$stage_level = NA
dat_suc_sp[dat_suc_sp$stage == 'introduce',]$stage_level = 1
dat_suc_sp[dat_suc_sp$stage == 'establish',]$stage_level = 2
dat_suc_sp[dat_suc_sp$stage == 'dominant',]$stage_level = 3
dat_suc_sp$stage_level = ordered(dat_suc_sp$stage_level)
str(dat_suc_sp)
save(dat_suc_sp,
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/dat_suc_sp.rdata')


#### Single regression results for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for all species ####
##### establishment predictive curves for md.ab #####
setwd("~/BSS_coexist_v1")
library(phyr)
library(tibble)
library(lme4)
require(ape)
library(scales)
library(ggthemes)
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

dat_suc_sp %>% group_by(estab) %>% 
  summarise_at(vars(mnd), c(length))

pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

dat_dom_sp = dat_suc_sp %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_sp %>% group_by(domin) %>% 
  summarise_at(vars(mnd), c(length))

### get establishment predict data for mnd.a
lincombs.data.estab.mnd.a_single = data.frame(mnd.a=seq(min(dat_suc_sp$mnd.a),
                                                        max(dat_suc_sp$mnd.a),
                                                        length=100))

lincombs.matrix.estab.mnd.a_single=model.matrix(~mnd.a,
                                                data=lincombs.data.estab.mnd.a_single)
lincombs.matrix.estab.mnd.a_single=as.data.frame(lincombs.matrix.estab.mnd.a_single)
lincombs.estab.mnd.a_single=inla.make.lincombs(lincombs.matrix.estab.mnd.a_single)

inla.model_lincombs.estab.mnd.a_single = pglmm(estab ~ mnd.a+#(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_suc_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.estab.mnd.a_single,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.estab.mnd.a_single = inla.model_lincombs.estab.mnd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnd.a_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mnd.a_single,
                                                               function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnd.a_single$lower=unlist(lapply(lincombs.posterior.estab.mnd.a_single,
                                                     function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnd.a_single$upper=unlist(lapply(lincombs.posterior.estab.mnd.a_single,
                                                     function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnd.a_single
save(lincombs.data.estab.mnd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnd.a
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnd.a_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(estab.mnd.a_single.logistic=ggplot(data=lincombs.data.estab.mnd.a_single,aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd.a, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",
             x=min(dat_suc_sp$mnd.a)+(max(dat_suc_sp$mnd.a)-min(dat_suc_sp$mnd.a))*0.05,
             y=c(0.90,0.77),
             label=c("italic(β)['ND'[ab]] == '6.15'",
                     "'95% CI' == '[4.60, 7.70]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_1() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)
#mnd.a    6.15    4.60    7.70



### get establishment predict data for mlgfd.a
lincombs.data.estab.mlgfd.a_single = data.frame(mlgfd.a=seq(min(dat_suc_sp$mlgfd.a),
                                                            max(dat_suc_sp$mlgfd.a),
                                                            length=100))

lincombs.matrix.estab.mlgfd.a_single=model.matrix(~mlgfd.a,
                                                  data=lincombs.data.estab.mlgfd.a_single)
lincombs.matrix.estab.mlgfd.a_single=as.data.frame(lincombs.matrix.estab.mlgfd.a_single)
lincombs.estab.mlgfd.a_single=inla.make.lincombs(lincombs.matrix.estab.mlgfd.a_single)

inla.model_lincombs.estab.mlgfd.a_single = pglmm(estab ~ mlgfd.a+#(1|species) + 
                                                   (1|f_p) + (1|field), data = dat_suc_sp,
                                                 family = "binomial", cov_ranef = list(species = tree),
                                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                             config = TRUE),
                                                                      quantiles=c(0.025,0.5,0.975),
                                                                      lincomb=lincombs.estab.mlgfd.a_single,
                                                                      control.predictor=list(compute=T)),
                                                 bayes = T)

lincombs.posterior.estab.mlgfd.a_single = inla.model_lincombs.estab.mlgfd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mlgfd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mlgfd.a_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mlgfd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mlgfd.a_single$lower=unlist(lapply(lincombs.posterior.estab.mlgfd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mlgfd.a_single$upper=unlist(lapply(lincombs.posterior.estab.mlgfd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mlgfd.a_single
save(lincombs.data.estab.mlgfd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mlgfd.a
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd.a_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(estab.mlgfd.a_single.logistic=ggplot(data=lincombs.data.estab.mlgfd.a_single,aes(x=mlgfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mlgfd.a, y=estab),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y=' ')+
    annotate(geom="text",
             x=min(dat_suc_sp$mlgfd.a)+(max(dat_suc_sp$mlgfd.a)-min(dat_suc_sp$mlgfd.a))*0.05,
             y=c(0.90,0.77),
             label=c("italic(β)['RFD'[ab]] == '2.93'",
                     "'95% CI' == '[2.49, 3.37]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1()
)
#mlgfd.a      2.93       2.49       3.37


### get establishment predict data for mpd.a_all
lincombs.data.estab.mpd.a_single = data.frame(mpd.a_all=seq(min(dat_suc_sp$mpd.a_all),
                                                            max(dat_suc_sp$mpd.a_all),
                                                            length=100))

lincombs.matrix.estab.mpd.a_single=model.matrix(~mpd.a_all,
                                                data=lincombs.data.estab.mpd.a_single)
lincombs.matrix.estab.mpd.a_single=as.data.frame(lincombs.matrix.estab.mpd.a_single)
lincombs.estab.mpd.a_single=inla.make.lincombs(lincombs.matrix.estab.mpd.a_single)

inla.model_lincombs.estab.mpd.a_single = pglmm(estab ~ mpd.a_all+#(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_suc_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.estab.mpd.a_single,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.estab.mpd.a_single = inla.model_lincombs.estab.mpd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mpd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mpd.a_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mpd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mpd.a_single$lower=unlist(lapply(lincombs.posterior.estab.mpd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mpd.a_single$upper=unlist(lapply(lincombs.posterior.estab.mpd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mpd.a_single
save(lincombs.data.estab.mpd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mpd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mpd.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mpd.a_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(estab.mpd.a_single.logistic=ggplot(data=lincombs.data.estab.mpd.a_single,aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd.a_all, y=estab),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",
             x=min(dat_suc_sp$mpd.a_all)+(max(dat_suc_sp$mpd.a_all)-min(dat_suc_sp$mpd.a_all))*0.1,
             y=c(0.90,0.77),
             label=c("italic(β)['MPD'[ab]] == '-0.0180'",
                     "'95% CI' == '[-0.0196, -0.0165]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)
# mpd.a_all   -0.0180    -0.0196    -0.0165


### get establishment predict data for mfunc_d.a_all
lincombs.data.estab.mfunc_d.a_single = data.frame(mfunc_d.a_all = seq(min(dat_suc_sp$mfunc_d.a_all),
                                                                      max(dat_suc_sp$mfunc_d.a_all),
                                                                      length=100))

lincombs.matrix.estab.mfunc_d.a_single=model.matrix(~mfunc_d.a_all,
                                                    data=lincombs.data.estab.mfunc_d.a_single)
lincombs.matrix.estab.mfunc_d.a_single=as.data.frame(lincombs.matrix.estab.mfunc_d.a_single)
lincombs.estab.mfunc_d.a_single=inla.make.lincombs(lincombs.matrix.estab.mfunc_d.a_single)

inla.model_lincombs.estab.mfunc_d.a_single = pglmm(estab ~ mfunc_d.a_all+#(1|species) + 
                                                     (1|f_p) + (1|field), data = dat_suc_sp,
                                                   family = "binomial", cov_ranef = list(species = tree),
                                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                               config = TRUE),
                                                                        quantiles=c(0.025,0.5,0.975),
                                                                        lincomb=lincombs.estab.mfunc_d.a_single,
                                                                        control.predictor=list(compute=T)),
                                                   bayes = T)

lincombs.posterior.estab.mfunc_d.a_single = inla.model_lincombs.estab.mfunc_d.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mfunc_d.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mfunc_d.a_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mfunc_d.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mfunc_d.a_single$lower=unlist(lapply(lincombs.posterior.estab.mfunc_d.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mfunc_d.a_single$upper=unlist(lapply(lincombs.posterior.estab.mfunc_d.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mfunc_d.a_single
save(lincombs.data.estab.mfunc_d.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mfunc_d.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfunc_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mfunc_d.a_single.rdata")
(estab.mfunc_d.a_single.logistic=ggplot(data=lincombs.data.estab.mfunc_d.a_single,
                                        aes(x=mfunc_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1,
              linetype= 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d.a_all, y=estab),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=c(0.25,0.25),y=c(0.9,0.80),
             label=c("italic(β)['MFD'[ab]] == '-3.3781'",
                     "'95%CI' == '[-4.2299, -2.5261]'"),
             parse=T,size=3.5)+
    theme_regular()
)


### get establishment predict data for mconti_func_d.a_all
lincombs.data.estab.mconti_func_d.a_single = data.frame(mconti_func_d.a_all = seq(min(dat_suc_sp$mconti_func_d.a_all),
                                                                                  max(dat_suc_sp$mconti_func_d.a_all),
                                                                                  length=100))

lincombs.matrix.estab.mconti_func_d.a_single=model.matrix(~mconti_func_d.a_all,
                                                          data=lincombs.data.estab.mconti_func_d.a_single)
lincombs.matrix.estab.mconti_func_d.a_single=as.data.frame(lincombs.matrix.estab.mconti_func_d.a_single)
lincombs.estab.mconti_func_d.a_single=inla.make.lincombs(lincombs.matrix.estab.mconti_func_d.a_single)

inla.model_lincombs.estab.mconti_func_d.a_single = pglmm(estab ~ mconti_func_d.a_all+#(1|species) + 
                                                           (1|f_p) + (1|field), data = dat_suc_sp,
                                                         family = "binomial", cov_ranef = list(species = tree),
                                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                     config = TRUE),
                                                                              quantiles=c(0.025,0.5,0.975),
                                                                              lincomb=lincombs.estab.mconti_func_d.a_single,
                                                                              control.predictor=list(compute=T)),
                                                         bayes = T)

lincombs.posterior.estab.mconti_func_d.a_single = inla.model_lincombs.estab.mconti_func_d.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mconti_func_d.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mconti_func_d.a_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mconti_func_d.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mconti_func_d.a_single$lower=unlist(lapply(lincombs.posterior.estab.mconti_func_d.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mconti_func_d.a_single$upper=unlist(lapply(lincombs.posterior.estab.mconti_func_d.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mconti_func_d.a_single
save(lincombs.data.estab.mconti_func_d.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mconti_func_d.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mconti_func_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mconti_func_d.a_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(estab.mconti_func_d.a_single.logistic=ggplot(data=lincombs.data.estab.mconti_func_d.a_single,
                                              aes(x=mconti_func_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1,
              linetype= 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d.a_all, y=estab),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",
             x=min(dat_suc_sp$mconti_func_d.a_all)+(max(dat_suc_sp$mconti_func_d.a_all)-min(dat_suc_sp$mconti_func_d.a_all))*0.1,
             y=c(0.9,0.77),
             label=c("italic(β)['MFD'[ab]] == '-0.6965'",
                     "'95% CI' == '[-0.8422, -0.5507]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1()
)

#mconti_func_d.a_all -0.6965    -0.8422    -0.5507


##### dominance predictive curves for md.ab #####
### get dominance predict data for mnd.a
lincombs.data.domin.mnd.a_single = data.frame(mnd.a=seq(min(dat_dom_sp$mnd.a),
                                                        max(dat_dom_sp$mnd.a),
                                                        length=100))

lincombs.matrix.domin.mnd.a_single=model.matrix(~mnd.a,
                                                data=lincombs.data.domin.mnd.a_single)
lincombs.matrix.domin.mnd.a_single=as.data.frame(lincombs.matrix.domin.mnd.a_single)
lincombs.domin.mnd.a_single=inla.make.lincombs(lincombs.matrix.domin.mnd.a_single)

inla.model_lincombs.domin.mnd.a_single = pglmm(domin ~ mnd.a+#(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_dom_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.domin.mnd.a_single,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.domin.mnd.a_single = inla.model_lincombs.domin.mnd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnd.a_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mnd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnd.a_single$lower=unlist(lapply(lincombs.posterior.domin.mnd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnd.a_single$upper=unlist(lapply(lincombs.posterior.domin.mnd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnd.a_single
save(lincombs.data.domin.mnd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnd.a
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnd.a_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin.mnd.a_single.logistic=ggplot(data=lincombs.data.domin.mnd.a_single,aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd.a, y=domin),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=expression("Niche difference (" * ND["ab"] * ")"), y='Dominance probability')+
    annotate(geom="text",
             x=min(dat_dom_sp$mnd.a)+(max(dat_dom_sp$mnd.a)-min(dat_dom_sp$mnd.a))*0.32,
             y=c(0.20,0.07),
             label=c("italic(β)['ND'[ab]] == '5.08'",
                     "'95% CI' == '[2.35, 7.81]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)
#mnd.a        5.08       2.35       7.81


### get dominance predict data for mlgfd.a
lincombs.data.domin.mlgfd.a_single = data.frame(mlgfd.a=seq(min(dat_dom_sp$mlgfd.a),
                                                            max(dat_dom_sp$mlgfd.a),
                                                            length=100))

lincombs.matrix.domin.mlgfd.a_single=model.matrix(~mlgfd.a,
                                                  data=lincombs.data.domin.mlgfd.a_single)
lincombs.matrix.domin.mlgfd.a_single=as.data.frame(lincombs.matrix.domin.mlgfd.a_single)
lincombs.domin.mlgfd.a_single=inla.make.lincombs(lincombs.matrix.domin.mlgfd.a_single)

inla.model_lincombs.domin.mlgfd.a_single = pglmm(domin ~ mlgfd.a+#(1|species) + 
                                                   (1|f_p) + (1|field), data = dat_dom_sp,
                                                 family = "binomial", cov_ranef = list(species = tree),
                                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                             config = TRUE),
                                                                      quantiles=c(0.025,0.5,0.975),
                                                                      lincomb=lincombs.domin.mlgfd.a_single,
                                                                      control.predictor=list(compute=T)),
                                                 bayes = T)

lincombs.posterior.domin.mlgfd.a_single = inla.model_lincombs.domin.mlgfd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mlgfd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mlgfd.a_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mlgfd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mlgfd.a_single$lower=unlist(lapply(lincombs.posterior.domin.mlgfd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mlgfd.a_single$upper=unlist(lapply(lincombs.posterior.domin.mlgfd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mlgfd.a_single
save(lincombs.data.domin.mlgfd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mlgfd.a
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd.a_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin.mlgfd.a_single.logistic=ggplot(data=lincombs.data.domin.mlgfd.a_single,aes(x=mlgfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mlgfd.a, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=expression("Relative fitness difference (" * RFD["ab"] * ")"), y='')+
    annotate(geom="text",
             x=min(dat_dom_sp$mlgfd.a)+(max(dat_dom_sp$mlgfd.a)-min(dat_dom_sp$mlgfd.a))*0.32,
             y=c(0.2,0.07),
             label=c("italic(β)['RFD'[ab]] == '1.32'",
                     "'95% CI' == '[0.70, 1.94]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)
#mlgfd.a      1.32       0.70       1.94



### get dominance predict data for mpd.a_all
lincombs.data.domin.mpd.a_single = data.frame(mpd.a_all=seq(min(dat_dom_sp$mpd.a_all),
                                                            max(dat_dom_sp$mpd.a_all),
                                                            length=100))

lincombs.matrix.domin.mpd.a_single=model.matrix(~mpd.a_all,
                                                data=lincombs.data.domin.mpd.a_single)
lincombs.matrix.domin.mpd.a_single=as.data.frame(lincombs.matrix.domin.mpd.a_single)
lincombs.domin.mpd.a_single=inla.make.lincombs(lincombs.matrix.domin.mpd.a_single)

inla.model_lincombs.domin.mpd.a_single = pglmm(domin ~ mpd.a_all
                                               #+(1|species) 
                                               +(1|f_p) + (1|field), data = dat_dom_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.domin.mpd.a_single,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.domin.mpd.a_single = inla.model_lincombs.domin.mpd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mpd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mpd.a_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mpd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mpd.a_single$lower=unlist(lapply(lincombs.posterior.domin.mpd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mpd.a_single$upper=unlist(lapply(lincombs.posterior.domin.mpd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mpd.a_single
save(lincombs.data.domin.mpd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mpd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mpd.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mpd.a_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin.mpd.a_single.logistic=ggplot(data=lincombs.data.domin.mpd.a_single,aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 3)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mpd.a_all, y=domin),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=expression("Phylogenetic difference (" * MPD["ab"] * ")"), y='Dominance probability')+
    annotate(geom="text",
             x=min(dat_dom_sp$mpd.a_all)+(max(dat_dom_sp$mpd.a_all)-min(dat_dom_sp$mpd.a_all))*0.1,
             y=c(0.90,0.77),
             label=c("italic(β)['MPD'[ab]] == '-0.0028'",
                     "'95% CI' == '[-0.0058, 0.0002]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)
# -0.0028    -0.0058     0.0002

### get dominance predict data for mfunc_d.a_all
lincombs.data.domin.mfunc_d.a_single = data.frame(mfunc_d.a_all = seq(min(dat_dom_sp$mfunc_d.a_all),
                                                                      max(dat_dom_sp$mfunc_d.a_all),length=100))

lincombs.matrix.domin.mfunc_d.a_single=model.matrix(~mfunc_d.a_all,
                                                    data=lincombs.data.domin.mfunc_d.a_single)
lincombs.matrix.domin.mfunc_d.a_single=as.data.frame(lincombs.matrix.domin.mfunc_d.a_single)
lincombs.domin.mfunc_d.a_single=inla.make.lincombs(lincombs.matrix.domin.mfunc_d.a_single)

inla.model_lincombs.domin.mfunc_d.a_single = pglmm(domin ~ mfunc_d.a_all+#(1|species) + 
                                                     (1|f_p) + (1|field), data = dat_dom_sp,
                                                   family = "binomial", cov_ranef = list(species = tree),
                                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                               config = TRUE),
                                                                        quantiles=c(0.025,0.5,0.975),
                                                                        lincomb=lincombs.domin.mfunc_d.a_single,
                                                                        control.predictor=list(compute=T)),
                                                   bayes = T)

lincombs.posterior.domin.mfunc_d.a_single = inla.model_lincombs.domin.mfunc_d.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mfunc_d.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mfunc_d.a_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mfunc_d.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mfunc_d.a_single$lower=unlist(lapply(lincombs.posterior.domin.mfunc_d.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mfunc_d.a_single$upper=unlist(lapply(lincombs.posterior.domin.mfunc_d.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mfunc_d.a_single
save(lincombs.data.domin.mfunc_d.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mfunc_d.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfunc_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mfunc_d.a_single.rdata")
(domin.mfunc_d.a_single.logistic=ggplot(data=lincombs.data.domin.mfunc_d.a_single,
                                        aes(x=mfunc_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_x_continuous(breaks = c(0.1, 0.3, 0.5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfunc_d.a_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=expression("Functional difference (" * MFD["ab"] * ")"), y='  ')+
    annotate(geom="text",x=c(0.25,0.25),y=c(0.90,0.80),
             label=c("italic(β)['MFD'[ab]] == '5.9743'",
                     "'95%CI' == '[4.2340, 7.7029]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mfunc_d.a_all  5.9743     4.2340     7.7029



### get dominance predict data for mconti_func_d.a_all
lincombs.data.domin.mconti_func_d.a_single = data.frame(mconti_func_d.a_all = seq(min(dat_dom_sp$mconti_func_d.a_all),
                                                                                  max(dat_dom_sp$mconti_func_d.a_all),
                                                                                  length=100))

lincombs.matrix.domin.mconti_func_d.a_single=model.matrix(~mconti_func_d.a_all,
                                                          data=lincombs.data.domin.mconti_func_d.a_single)
lincombs.matrix.domin.mconti_func_d.a_single=as.data.frame(lincombs.matrix.domin.mconti_func_d.a_single)
lincombs.domin.mconti_func_d.a_single=inla.make.lincombs(lincombs.matrix.domin.mconti_func_d.a_single)

inla.model_lincombs.domin.mconti_func_d.a_single = pglmm(domin ~ mconti_func_d.a_all+#(1|species) + 
                                                           (1|f_p) + (1|field), data = dat_dom_sp,
                                                         family = "binomial",
                                                         cov_ranef = list(species = tree),
                                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                     config = TRUE),
                                                                              quantiles=c(0.025,0.5,0.975),
                                                                              lincomb=lincombs.domin.mconti_func_d.a_single,
                                                                              control.predictor=list(compute=T)),
                                                         bayes = T)

lincombs.posterior.domin.mconti_func_d.a_single = inla.model_lincombs.domin.mconti_func_d.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mconti_func_d.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mconti_func_d.a_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mconti_func_d.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mconti_func_d.a_single$lower=unlist(lapply(lincombs.posterior.domin.mconti_func_d.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mconti_func_d.a_single$upper=unlist(lapply(lincombs.posterior.domin.mconti_func_d.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mconti_func_d.a_single
save(lincombs.data.domin.mconti_func_d.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mconti_func_d.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mconti_func_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mconti_func_d.a_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin.mconti_func_d.a_single.logistic=ggplot(data=lincombs.data.domin.mconti_func_d.a_single,
                                              aes(x=mconti_func_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size = 1,
              linetype= 3)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mconti_func_d.a_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=expression("Functional difference (" * MFD["ab"] * ")"), y='  ')+
    annotate(geom="text",
             x=min(dat_dom_sp$mconti_func_d.a_all)+(max(dat_dom_sp$mconti_func_d.a_all)-min(dat_dom_sp$mconti_func_d.a_all))*0.1,
             y=c(0.9,0.77),
             label=c("italic(β)['MFD'[ab]] == '-0.0776'",
                     "'95% CI' == '[-0.4201, 0.2635]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)
#mconti_func_d.a_all -0.0776    -0.4201     0.2635


##### establishment predictive curves for md #####
setwd("~/BSS_coexist_v1")
library(phyr)
library(tibble)
library(lme4)
require(ape)
library(scales)
library(ggthemes)
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

dat_dom_sp = dat_suc_sp %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c('establish', 'dominant'))


### get establishment predict data for mnd
lincombs.data.estab.mnd_single = data.frame(mnd=seq(min(dat_suc_sp$mnd),
                                                    max(dat_suc_sp$mnd),
                                                    length=100))

lincombs.matrix.estab.mnd_single=model.matrix(~mnd,
                                              data=lincombs.data.estab.mnd_single)
lincombs.matrix.estab.mnd_single=as.data.frame(lincombs.matrix.estab.mnd_single)
lincombs.estab.mnd_single=inla.make.lincombs(lincombs.matrix.estab.mnd_single)

inla.model_lincombs.estab.mnd_single = pglmm(estab ~ mnd+#(1|species) + 
                                               (1|f_p) + (1|field), data = dat_suc_sp,
                                             family = "binomial", cov_ranef = list(species = tree),
                                             bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                         config = TRUE),
                                                                  quantiles=c(0.025,0.5,0.975),
                                                                  lincomb=lincombs.estab.mnd_single,
                                                                  control.predictor=list(compute=T)),
                                             bayes = T)

lincombs.posterior.estab.mnd_single = inla.model_lincombs.estab.mnd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnd_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mnd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnd_single$lower=unlist(lapply(lincombs.posterior.estab.mnd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnd_single$upper=unlist(lapply(lincombs.posterior.estab.mnd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnd_single
save(lincombs.data.estab.mnd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnd_single.rdata")

ggthemr::ggthemr(palette = "fresh", layout = "clean")

(estab.mnd_single.logistic=ggplot(data=lincombs.data.estab.mnd_single,aes(x=mnd,
                                                                            y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_x_continuous(limits = c(ifelse(min(dat_suc_sp$mnd) > 0,
                                         min(dat_suc_sp$mnd) * 0.1,
                                         min(dat_suc_sp$mnd) * 1.9  
    ),
    max(dat_suc_sp$mnd)))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",
             x=min(dat_suc_sp$mnd) * 1.9   +
               (max(dat_suc_sp$mnd)-min(dat_suc_sp$mnd) * 1.9 )*0.01,
             y=c(0.90,0.77),
             label=c("italic(β)['ND'] == '12.87'",
                     "'95% CI' == '[9.98, 15.76]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_1() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)
#mnd         12.87       9.98      15.76



### get establishment predict data for mlgfd
lincombs.data.estab.mlgfd_single = data.frame(mlgfd=seq(min(dat_suc_sp$mlgfd),
                                                        max(dat_suc_sp$mlgfd),
                                                        length=100))

lincombs.matrix.estab.mlgfd_single=model.matrix(~mlgfd,
                                                data=lincombs.data.estab.mlgfd_single)
lincombs.matrix.estab.mlgfd_single=as.data.frame(lincombs.matrix.estab.mlgfd_single)
lincombs.estab.mlgfd_single=inla.make.lincombs(lincombs.matrix.estab.mlgfd_single)

inla.model_lincombs.estab.mlgfd_single = pglmm(estab ~ mlgfd+#(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_suc_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.estab.mlgfd_single,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.estab.mlgfd_single = inla.model_lincombs.estab.mlgfd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mlgfd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mlgfd_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mlgfd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mlgfd_single$lower=unlist(lapply(lincombs.posterior.estab.mlgfd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mlgfd_single$upper=unlist(lapply(lincombs.posterior.estab.mlgfd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mlgfd_single
save(lincombs.data.estab.mlgfd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mlgfd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(estab.mlgfd_single.logistic=ggplot(data=lincombs.data.estab.mlgfd_single,aes(x=mlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mlgfd, y=estab),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y=' ')+
    scale_x_continuous(limits = c(ifelse(min(dat_suc_sp$mlgfd) > 0,
                                         min(dat_suc_sp$mlgfd)*0.2,
                                         min(dat_suc_sp$mlgfd)*1.8   
    ),
                                  max(dat_suc_sp$mlgfd)))+
    annotate(geom="text",
             x=min(dat_suc_sp$mlgfd)*1.8+
               (max(dat_suc_sp$mlgfd)-min(dat_suc_sp$mlgfd)*1.8)*0.01,
             y=c(0.90,0.77),
             label=c("italic(β)['RFD'] == '3.79'",
                     "'95% CI' == '[3.27, 4.32]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1()
)
#mlgfd        3.79       3.27       4.32


### get establishment predict data for mpd_all
lincombs.data.estab.mpd_single = data.frame(mpd_all=seq(min(dat_suc_sp$mpd_all),
                                                        max(dat_suc_sp$mpd_all),
                                                        length=100))

lincombs.matrix.estab.mpd_single=model.matrix(~mpd_all,
                                              data=lincombs.data.estab.mpd_single)
lincombs.matrix.estab.mpd_single=as.data.frame(lincombs.matrix.estab.mpd_single)
lincombs.estab.mpd_single=inla.make.lincombs(lincombs.matrix.estab.mpd_single)

inla.model_lincombs.estab.mpd_single = pglmm(estab ~ mpd_all+#(1|species) + 
                                               (1|f_p) + (1|field), data = dat_suc_sp,
                                             family = "binomial", cov_ranef = list(species = tree),
                                             bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                         config = TRUE),
                                                                  quantiles=c(0.025,0.5,0.975),
                                                                  lincomb=lincombs.estab.mpd_single,
                                                                  control.predictor=list(compute=T)),
                                             bayes = T)

lincombs.posterior.estab.mpd_single = inla.model_lincombs.estab.mpd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mpd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mpd_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mpd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mpd_single$lower=unlist(lapply(lincombs.posterior.estab.mpd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mpd_single$upper=unlist(lapply(lincombs.posterior.estab.mpd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mpd_single
save(lincombs.data.estab.mpd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mpd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mpd_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mpd_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(estab.mpd_single.logistic=ggplot(data=lincombs.data.estab.mpd_single,aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=estab),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",
             x=min(dat_suc_sp$mpd_all)+(max(dat_suc_sp$mpd_all)-min(dat_suc_sp$mpd_all))*0.1,
             y=c(0.90,0.77),
             label=c("italic(β)['MPD'] == '-0.0261'",
                     "'95% CI' == '[-0.0284, -0.0238]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)

#mpd_all     -0.0261    -0.0284    -0.0238



### get establishment predict data for mconti_func_d_all
lincombs.data.estab.mconti_func_d_single = data.frame(mconti_func_d_all = seq(min(dat_suc_sp$mconti_func_d_all),
                                                                              max(dat_suc_sp$mconti_func_d_all),
                                                                              length=100))

lincombs.matrix.estab.mconti_func_d_single=model.matrix(~mconti_func_d_all,
                                                        data=lincombs.data.estab.mconti_func_d_single)
lincombs.matrix.estab.mconti_func_d_single=as.data.frame(lincombs.matrix.estab.mconti_func_d_single)
lincombs.estab.mconti_func_d_single=inla.make.lincombs(lincombs.matrix.estab.mconti_func_d_single)

inla.model_lincombs.estab.mconti_func_d_single = pglmm(estab ~ mconti_func_d_all+#(1|species) + 
                                                         (1|f_p) + (1|field), data = dat_suc_sp,
                                                       family = "binomial", cov_ranef = list(species = tree),
                                                       bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                   config = TRUE),
                                                                            quantiles=c(0.025,0.5,0.975),
                                                                            lincomb=lincombs.estab.mconti_func_d_single,
                                                                            control.predictor=list(compute=T)),
                                                       bayes = T)

lincombs.posterior.estab.mconti_func_d_single = inla.model_lincombs.estab.mconti_func_d_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mconti_func_d_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mconti_func_d_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mconti_func_d_single,
                                                                       function(x) inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mconti_func_d_single$lower=unlist(lapply(lincombs.posterior.estab.mconti_func_d_single,
                                                             function(x) inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mconti_func_d_single$upper=unlist(lapply(lincombs.posterior.estab.mconti_func_d_single,
                                                             function(x) inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mconti_func_d_single
save(lincombs.data.estab.mconti_func_d_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mconti_func_d_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mconti_func_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mconti_func_d_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(estab.mconti_func_d_single.logistic=ggplot(data=lincombs.data.estab.mconti_func_d_single,
                                              aes(x=mconti_func_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1,
              linetype= 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d_all, y=estab),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",
             x=min(dat_suc_sp$mconti_func_d_all)+(max(dat_suc_sp$mconti_func_d_all)-min(dat_suc_sp$mconti_func_d_all))*0.1,
             y=c(0.9,0.77),
             label=c("italic(β)['MFD'] == '-0.1117'",
                     "'95% CI' == '[-0.3285, -0.1051]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1()
)
#mconti_func_d_all -0.1117    -0.3285     0.1051


### get establishment predict data for mfunc_d_all
lincombs.data.estab.mfunc_d_single = data.frame(mfunc_d_all = seq(min(dat_suc_sp$mfunc_d_all),
                                                                  max(dat_suc_sp$mfunc_d_all),
                                                                  length=100))

lincombs.matrix.estab.mfunc_d_single=model.matrix(~mfunc_d_all,
                                                  data=lincombs.data.estab.mfunc_d_single)
lincombs.matrix.estab.mfunc_d_single=as.data.frame(lincombs.matrix.estab.mfunc_d_single)
lincombs.estab.mfunc_d_single=inla.make.lincombs(lincombs.matrix.estab.mfunc_d_single)

inla.model_lincombs.estab.mfunc_d_single = pglmm(estab ~ mfunc_d_all+#(1|species) + 
                                                   (1|f_p) + (1|field), data = dat_suc_sp,
                                                 family = "binomial", cov_ranef = list(species = tree),
                                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                             config = TRUE),
                                                                      quantiles=c(0.025,0.5,0.975),
                                                                      lincomb=lincombs.estab.mfunc_d_single,
                                                                      control.predictor=list(compute=T)),
                                                 bayes = T)

lincombs.posterior.estab.mfunc_d_single = inla.model_lincombs.estab.mfunc_d_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mfunc_d_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mfunc_d_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mfunc_d_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mfunc_d_single$lower=unlist(lapply(lincombs.posterior.estab.mfunc_d_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mfunc_d_single$upper=unlist(lapply(lincombs.posterior.estab.mfunc_d_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mfunc_d_single
save(lincombs.data.estab.mfunc_d_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mfunc_d_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mfunc_d_single.rdata")
(estab.mfunc_d_single.logistic=ggplot(data=lincombs.data.estab.mfunc_d_single,
                                      aes(x=mfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1,
              linetype= 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d_all, y=estab),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=c(0.3,0.3),y=c(0.9,0.80),
             label=c("italic(β)['MFD'] == '7.0488'",
                     "'95%CI' == '[6.0936, 8.0041]'"),
             parse=T,size=3.5)+
    theme_regular()
)


##### dominance predictive curves for md #####
### get dominance predict data for mnd
lincombs.data.domin.mnd_single = data.frame(mnd=seq(min(dat_dom_sp$mnd),
                                                    max(dat_dom_sp$mnd),
                                                    length=100))

lincombs.matrix.domin.mnd_single=model.matrix(~mnd,
                                              data=lincombs.data.domin.mnd_single)
lincombs.matrix.domin.mnd_single=as.data.frame(lincombs.matrix.domin.mnd_single)
lincombs.domin.mnd_single=inla.make.lincombs(lincombs.matrix.domin.mnd_single)

inla.model_lincombs.domin.mnd_single = pglmm(domin ~ mnd+#(1|species) + 
                                               (1|f_p) + (1|field), data = dat_dom_sp,
                                             family = "binomial", cov_ranef = list(species = tree),
                                             bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                         config = TRUE),
                                                                  quantiles=c(0.025,0.5,0.975),
                                                                  lincomb=lincombs.domin.mnd_single,
                                                                  control.predictor=list(compute=T)),
                                             bayes = T)

lincombs.posterior.domin.mnd_single = inla.model_lincombs.domin.mnd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnd_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mnd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnd_single$lower=unlist(lapply(lincombs.posterior.domin.mnd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnd_single$upper=unlist(lapply(lincombs.posterior.domin.mnd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnd_single
save(lincombs.data.domin.mnd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnd_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin.mnd_single.logistic=ggplot(data=lincombs.data.domin.mnd_single,
                                  aes(x=mnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_x_continuous(limits = c(ifelse(min(dat_suc_sp$mnd) > 0,
                                         min(dat_suc_sp$mnd) * 0.1,
                                         min(dat_suc_sp$mnd) * 1.9  
    ),
    max(dat_suc_sp$mnd)))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd, y=domin),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Niche difference (ND)', y='Dominance probability')+
    annotate(geom="text",
             x=min(dat_suc_sp$mnd) * 1.9   +
               (max(dat_dom_sp$mnd)-min(dat_suc_sp$mnd) * 1.9)*0.001,
             y=c(0.90,0.77),
             label=c("italic(β)['ND'] == '11.65'",
                     "'95% CI' == '[6.57, 16.72]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)
#mnd         11.65       6.57      16.72



### get dominance predict data for mlgfd
lincombs.data.domin.mlgfd_single = data.frame(mlgfd=seq(min(dat_dom_sp$mlgfd),
                                                        max(dat_dom_sp$mlgfd),
                                                        length=100))

lincombs.matrix.domin.mlgfd_single=model.matrix(~mlgfd,
                                                data=lincombs.data.domin.mlgfd_single)
lincombs.matrix.domin.mlgfd_single=as.data.frame(lincombs.matrix.domin.mlgfd_single)
lincombs.domin.mlgfd_single=inla.make.lincombs(lincombs.matrix.domin.mlgfd_single)

inla.model_lincombs.domin.mlgfd_single = pglmm(domin ~ mlgfd+#(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_dom_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.domin.mlgfd_single,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.domin.mlgfd_single = inla.model_lincombs.domin.mlgfd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mlgfd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mlgfd_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mlgfd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mlgfd_single$lower=unlist(lapply(lincombs.posterior.domin.mlgfd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mlgfd_single$upper=unlist(lapply(lincombs.posterior.domin.mlgfd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mlgfd_single
save(lincombs.data.domin.mlgfd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mlgfd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin.mlgfd_single.logistic=ggplot(data=lincombs.data.domin.mlgfd_single,aes(x=mlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mlgfd, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    scale_x_continuous(limits = c(ifelse(min(dat_dom_sp$mlgfd) > 0,
                                         min(dat_dom_sp$mlgfd)*0.2,
                                         min(dat_dom_sp$mlgfd)*2.5   
    ),
    max(dat_dom_sp$mlgfd)))+
    labs(x='Relative fitness difference (RFD)', y='')+
    annotate(geom="text",
             x=min(dat_dom_sp$mlgfd)*2.5+
               (max(dat_dom_sp$mlgfd)-min(dat_dom_sp$mlgfd)*2.5)*0.32,
             y=c(0.2,0.07),
             label=c("italic(β)['RFD'] == '1.88'",
                     "'95% CI' == '[1.10, 2.67]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)
#mlgfd        1.88       1.10       2.67



### get dominance predict data for mpd_all
lincombs.data.domin.mpd_single = data.frame(mpd_all=seq(min(dat_dom_sp$mpd_all),
                                                        max(dat_dom_sp$mpd_all),
                                                        length=100))

lincombs.matrix.domin.mpd_single=model.matrix(~mpd_all,
                                              data=lincombs.data.domin.mpd_single)
lincombs.matrix.domin.mpd_single=as.data.frame(lincombs.matrix.domin.mpd_single)
lincombs.domin.mpd_single=inla.make.lincombs(lincombs.matrix.domin.mpd_single)

inla.model_lincombs.domin.mpd_single = pglmm(domin ~ mpd_all+#(1|species) + 
                                               (1|f_p) + (1|field), data = dat_dom_sp,
                                             family = "binomial", cov_ranef = list(species = tree),
                                             bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                         config = TRUE),
                                                                  quantiles=c(0.025,0.5,0.975),
                                                                  lincomb=lincombs.domin.mpd_single,
                                                                  control.predictor=list(compute=T)),
                                             bayes = T)

lincombs.posterior.domin.mpd_single = inla.model_lincombs.domin.mpd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mpd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mpd_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mpd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mpd_single$lower=unlist(lapply(lincombs.posterior.domin.mpd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mpd_single$upper=unlist(lapply(lincombs.posterior.domin.mpd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mpd_single
save(lincombs.data.domin.mpd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mpd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mpd_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mpd_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin.mpd_single.logistic=ggplot(data=lincombs.data.domin.mpd_single,aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 3)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mpd_all, y=domin),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Phylogenetic difference (MPD)', y='Dominance probability')+
    annotate(geom="text",
             x=min(dat_dom_sp$mpd_all)+(max(dat_dom_sp$mpd_all)-min(dat_dom_sp$mpd_all))*0.1,
             y=c(0.90,0.77),
             label=c("italic(β)['MPD'] == '-0.0015'",
                     "'95% CI' == '[-0.0062, 0.0032]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)
#mpd_all     -0.0015    -0.0062     0.0032


### get dominance predict data for mconti_func_d_all
lincombs.data.domin.mconti_func_d_single = data.frame(mconti_func_d_all = seq(min(dat_dom_sp$mconti_func_d_all),
                                                                              max(dat_dom_sp$mconti_func_d_all),length=100))

lincombs.matrix.domin.mconti_func_d_single=model.matrix(~mconti_func_d_all,
                                                        data=lincombs.data.domin.mconti_func_d_single)
lincombs.matrix.domin.mconti_func_d_single=as.data.frame(lincombs.matrix.domin.mconti_func_d_single)
lincombs.domin.mconti_func_d_single=inla.make.lincombs(lincombs.matrix.domin.mconti_func_d_single)

inla.model_lincombs.domin.mconti_func_d_single = pglmm(domin ~ mconti_func_d_all+#(1|species) + 
                                                         (1|f_p) + (1|field), data = dat_dom_sp,
                                                       family = "binomial", cov_ranef = list(species = tree),
                                                       bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                   config = TRUE),
                                                                            quantiles=c(0.025,0.5,0.975),
                                                                            lincomb=lincombs.domin.mconti_func_d_single,
                                                                            control.predictor=list(compute=T)),
                                                       bayes = T)

lincombs.posterior.domin.mconti_func_d_single = inla.model_lincombs.domin.mconti_func_d_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mconti_func_d_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mconti_func_d_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mconti_func_d_single,
                                                                       function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mconti_func_d_single$lower=unlist(lapply(lincombs.posterior.domin.mconti_func_d_single,
                                                             function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mconti_func_d_single$upper=unlist(lapply(lincombs.posterior.domin.mconti_func_d_single,
                                                             function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mconti_func_d_single
save(lincombs.data.domin.mconti_func_d_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mconti_func_d_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mconti_func_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mconti_func_d_single.rdata")
ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin.mconti_func_d_single.logistic=ggplot(data=lincombs.data.domin.mconti_func_d_single,
                                              aes(x=mconti_func_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size = 1,
              linetype= 3)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mconti_func_d_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Functional difference (MFD)', y='  ')+
    annotate(geom="text",
             x=min(dat_dom_sp$mconti_func_d_all)+(max(dat_dom_sp$mconti_func_d_all)-min(dat_dom_sp$mconti_func_d_all))*0.1,
             y=c(0.9,0.77),
             label=c("italic(β)['MFD'] == '-0.2371'",
                     "'95% CI' == '[-0.7147, 0.2411]'"),
             parse=T,size=5, hjust = 0)+
    theme_regular_1() 
)
#mconti_func_d_all -0.2371    -0.7147     0.2411



### get dominance predict data for mfunc_d_all
lincombs.data.domin.mfunc_d_single = data.frame(mfunc_d_all = seq(min(dat_dom_sp$mfunc_d_all),
                                                                  max(dat_dom_sp$mfunc_d_all),length=100))

lincombs.matrix.domin.mfunc_d_single=model.matrix(~mfunc_d_all,
                                                  data=lincombs.data.domin.mfunc_d_single)
lincombs.matrix.domin.mfunc_d_single=as.data.frame(lincombs.matrix.domin.mfunc_d_single)
lincombs.domin.mfunc_d_single=inla.make.lincombs(lincombs.matrix.domin.mfunc_d_single)

inla.model_lincombs.domin.mfunc_d_single = pglmm(domin ~ mfunc_d_all+#(1|species) + 
                                                   (1|f_p) + (1|field), data = dat_dom_sp,
                                                 family = "binomial", cov_ranef = list(species = tree),
                                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                             config = TRUE),
                                                                      quantiles=c(0.025,0.5,0.975),
                                                                      lincomb=lincombs.domin.mfunc_d_single,
                                                                      control.predictor=list(compute=T)),
                                                 bayes = T)

lincombs.posterior.domin.mfunc_d_single = inla.model_lincombs.domin.mfunc_d_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mfunc_d_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mfunc_d_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mfunc_d_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mfunc_d_single$lower=unlist(lapply(lincombs.posterior.domin.mfunc_d_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mfunc_d_single$upper=unlist(lapply(lincombs.posterior.domin.mfunc_d_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mfunc_d_single
save(lincombs.data.domin.mfunc_d_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mfunc_d_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mfunc_d_single.rdata")
(domin.mfunc_d_single.logistic=ggplot(data=lincombs.data.domin.mfunc_d_single,
                                      aes(x=mfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_x_continuous(breaks = c(0.2, 0.3, 0.4))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfunc_d_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Functional difference (FD)', y='  ')+
    annotate(geom="text",x=c(0.28,0.28),y=c(0.90,0.80),
             label=c("italic(β)['MFD'] == '9.0397'",
                     "'95%CI' == '[7.1790, 10.9020]'"),
             parse=T,size=3.5)+
    theme_regular()
)

estab.mnd_single.logistic
estab.mlgfd_single.logistic
estab.mpd_single.logistic
estab.mconti_func_d_single.logistic
domin.mnd_single.logistic
domin.mlgfd_single.logistic
domin.mpd_single.logistic
domin.mconti_func_d_single.logistic


#### Fast start for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for all species ####
setwd("~/BSS_coexist_v1")
library(phyr)
library(tibble)
library(lme4)
require(ape)
library(scales)
library(ggthemes)
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

##### Estab #####
pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

## Check the co-linearity
### only continuous functional trait distance
car::vif(
  glmer(estab ~ mnd + mlgfd + mpd_all + mconti_func_d_all + (1|f_p) + (1|field),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnd.a + mlgfd.a + mpd.a_all + mconti_func_d.a_all + (1|f_p)+ (1|field),
        family = binomial, data = dat_suc_sps)
)


###no co-linearity problem, all VIF < 3

estab_model_md_conti_func_d_all = pglmm(estab~mnd+mlgfd+mpd_all+mconti_func_d_all
                                        #+(1|species) 
                                        + (1|f_p) + (1|field), data = dat_suc_sps,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975)),
                                        bayes = T)



estab_model_md.a_conti_func_d_all = pglmm(estab~mnd.a+mlgfd.a+mpd.a_all+mconti_func_d.a_all
                                          #+(1|species) 
                                          + (1|f_p) + (1|field), data = dat_suc_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)


summary(estab_model_md_conti_func_d_all)
summary(estab_model_md.a_conti_func_d_all)

# plot 
estab_data.inla.md_conti_func_d.all_intercept1 = estab_model_md_conti_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.md_conti_func_d.all_intercept = estab_data.inla.md_conti_func_d.all_intercept1%>%
  mutate(rowname=c("MND","MRFD","MPD_all","MFD_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

estab_data.inla.md.a_conti_func_d.all_intercept1 = estab_model_md.a_conti_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.md.a_conti_func_d.all_intercept = estab_data.inla.md.a_conti_func_d.all_intercept1%>%
  mutate(rowname=c("MND.ab","MRFD.ab","MPD.ab_all","MFD.ab_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))


### Effect size plot for mean differences
# point + effect size
require(ggplot2)
require(ggpubr)

estab_data.inla.md_conti_func_d.all_intercept_ndrfd = 
  estab_data.inla.md_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MRFD', 'MND'))
x_min = min(estab_data.inla.md_conti_func_d.all_intercept_ndrfd$lower)
x_max = max(estab_data.inla.md_conti_func_d.all_intercept_ndrfd$upper)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
(
  estab_nd_rfd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=estab_data.inla.md_conti_func_d.all_intercept_ndrfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=estab_data.inla.md_conti_func_d.all_intercept_ndrfd,
               mapping = aes(x=mean,y=rowname),
               color = 'black',
               size = 3.5)+
    geom_errorbar(data=estab_data.inla.md_conti_func_d.all_intercept_ndrfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper),
                  color = 'black',
                  width = 0.12, linewidth = 1)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD', 'MND'),
                     label = c(expression(RFD), expression(ND)))+
    #scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
        #                       x_max+0.1
   # ))+
    scale_x_continuous(limits=c(-0.15,
                                0.61
    ))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND" = "MND",
                                  "MRFD" = "MRFD"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[1],
                                 colors_4d[2]),
                      labels = c("MND" = "MND",
                                 "MRFD" = "MRFD"),
                      name = ' ')+
    labs(x = '  ', y = '  ')+
    guides(color="none")+
    theme_regular_1() 
)

estab_data.inla.md_conti_func_d.all_intercept_pdfd = 
  estab_data.inla.md_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MFD_all', 'MPD_all'))
x_min = min(estab_data.inla.md_conti_func_d.all_intercept_pdfd$lower)
x_max = max(estab_data.inla.md_conti_func_d.all_intercept_pdfd$upper)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
(estab_pd_fd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=estab_data.inla.md_conti_func_d.all_intercept_pdfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=estab_data.inla.md_conti_func_d.all_intercept_pdfd,
               mapping = aes(x=mean,y=rowname),
               color = 'black',
               size = 3.5)+
    geom_errorbar(data=estab_data.inla.md_conti_func_d.all_intercept_pdfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12,
                  color = 'black',linewidth = 1)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD_all','MPD_all'),
                     label=c(expression(MFD), expression(MPD)))+
    #scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
       #                        ifelse(x_max<0,0.15,x_max+0.15)
    #))+
    scale_x_continuous(limits=c(-0.83,
                                0.2
    ), breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD_all" = "MPD_all",
                                  "MFD_all" = "MFD_all"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[4],
                                 colors_4d[3]),
                      labels = c("MPD_all" = "MPD_all",
                                 "MFD_all" = "MFD_all"),
                      name = ' ')+
    labs(x = '  ', y = ' ')+
    guides(color="none")+
    theme_regular_1() 
)


### Effect size plot for abundance weighted mean differences
# point + effect size
estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd = 
  estab_data.inla.md.a_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MRFD.ab', 'MND.ab'))
x_min = min(estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd$lower)
x_max = max(estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd$upper)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
(
  estab_nd_rfd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
               mapping = aes(x=mean,y=rowname),
               color = 'black',
               size = 3.5)+
    geom_errorbar(data=estab_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper),
                  color = 'black',
                  width = 0.12, linewidth = 1)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD.ab', 'MND.ab'),
                     label = c(expression(RFD['ab']), expression(ND['ab'])))+
    #scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
      #                         x_max+0.1
    #))+
    scale_x_continuous(limits=c(-0.2,
                                0.6
    ))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND.ab" = "MND.ab",
                                  "MRFD.ab" = "MRFD.ab"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[1],
                                 colors_4d[2]),
                      labels = c("MND.ab" = "MND.ab",
                                 "MRFD.ab" = "MRFD.ab"),
                      name = ' ')+
    labs(x = '  ', y = '  ')+
    guides(color="none")+
    theme_regular_1() 
)

estab_data.inla.md.a_conti_func_d.all_intercept_pdfd = 
  estab_data.inla.md.a_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all'))
x_min = min(estab_data.inla.md.a_conti_func_d.all_intercept_pdfd$lower)
x_max = max(estab_data.inla.md.a_conti_func_d.all_intercept_pdfd$upper)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
(estab_pd_fd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=estab_data.inla.md.a_conti_func_d.all_intercept_pdfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=estab_data.inla.md.a_conti_func_d.all_intercept_pdfd,
               mapping = aes(x=mean,y=rowname),
               color = 'black',
               size = 3.5)+
    geom_errorbar(data=estab_data.inla.md.a_conti_func_d.all_intercept_pdfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12,
                  color = 'black',linewidth = 1)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all'),
                     label=c(expression(MFD['ab']), expression(MPD['ab'])))+
    #scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
      #                         ifelse(x_max<0,0.15,x_max+0.15)
   # ))+
    scale_x_continuous(limits=c(-0.83,
                                0.2
    ), breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD.ab_all" = "MPD.ab_all",
                                  "MFD.ab_all" = "MFD.ab_all"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[4],
                                 colors_4d[3]),
                      labels = c("MPD.ab_all" = "MPD.ab_all",
                                 "MFD.ab_all" = "MFD.ab_all"),
                      name = ' ')+
    labs(x = '  ', y = ' ')+
    guides(color="none")+
    theme_regular_1() 
)


##### Domin #####
dat_dom_sp = dat_suc_sp %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c('establish', 'dominant'))

## Check the co-linearity
car::vif(
  glmer(domin ~ mnd + mlgfd + mpd_all + mconti_func_d_all + 
          (1|f_p) + (1|field),
        family=binomial,data=dat_dom_sps)
)

car::vif(
  glmer(domin ~ mnd.a + mlgfd.a + mpd.a_all+mconti_func_d.a_all + 
          (1|f_p)+ (1|field),
        family=binomial,data=dat_dom_sps)
)


### no colinearity problem

domin_model_md_conti_func_d_all = pglmm(domin~mnd+mlgfd+mpd_all+mconti_func_d_all
                                        #+(1|species) 
                                        +  (1|f_p) + (1|field), data = dat_dom_sps,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975)),
                                        bayes = T)

domin_model_md.a_conti_func_d_all = pglmm(domin~mnd.a+mlgfd.a+mpd.a_all+mconti_func_d.a_all
                                          #+(1|species)
                                          + (1|f_p)
                                          + (1|field), data = dat_dom_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)

domin_model_mnd_conti_func_d_all  = pglmm(domin~mnnd+mnlgfd+mntd_all+mnconti_func_d_all
                                          #+(1|species)
                                          + (1|f_p) + (1|field), data = dat_dom_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,
                                                                                      waic=T,
                                                                                      cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)
summary(domin_model_md_conti_func_d_all)
summary(domin_model_md.a_conti_func_d_all)


# plot 
domin_data.inla.md_conti_func_d.all_intercept1 = domin_model_md_conti_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.md_conti_func_d.all_intercept = domin_data.inla.md_conti_func_d.all_intercept1%>%
  mutate(rowname=c("MND","MRFD","MPD_all","MFD_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean))) %>% 
  mutate(significant = ifelse(lower / abs(lower) == upper / abs(upper), 'yes', 'no'))

domin_data.inla.md.a_conti_func_d.all_intercept1 = domin_model_md.a_conti_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.md.a_conti_func_d.all_intercept = domin_data.inla.md.a_conti_func_d.all_intercept1%>%
  mutate(rowname=c("MND.ab","MRFD.ab","MPD.ab_all","MFD.ab_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))%>% 
  mutate(significant = ifelse(lower / abs(lower) == upper / abs(upper), 'yes', 'no'))

### Effect size plot for mean differences
# point + effect size
domin_data.inla.md_conti_func_d.all_intercept_ndrfd = 
  domin_data.inla.md_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MRFD', 'MND'))
x_min = min(domin_data.inla.md_conti_func_d.all_intercept_ndrfd$lower)
x_max = max(domin_data.inla.md_conti_func_d.all_intercept_ndrfd$upper)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
(
  domin_nd_rfd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=domin_data.inla.md_conti_func_d.all_intercept_ndrfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=domin_data.inla.md_conti_func_d.all_intercept_ndrfd,
               mapping = aes(x=mean,y=rowname),
               color = 'black',
               size = 3.5)+
    geom_errorbar(data=domin_data.inla.md_conti_func_d.all_intercept_ndrfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper),
                  color = 'black',
                  width = 0.12, linewidth = 1)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD', 'MND'),
                     label = c(expression(RFD), expression(ND)))+
    #scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
       #                         x_max+0.1
    #))+
    scale_x_continuous(limits=c(-0.15,
                             0.61
    ))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND" = "MND",
                                  "MRFD" = "MRFD"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[1],
                                 colors_4d[2]),
                      labels = c("MND" = "MND",
                                 "MRFD" = "MRFD"),
                      name = ' ')+
    labs(x = 'Standardized effects', y = '  ')+
    guides(color="none")+
    theme_regular_1() 
)

domin_data.inla.md_conti_func_d.all_intercept_pdfd = 
  domin_data.inla.md_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MFD_all', 'MPD_all'))
x_min = min(domin_data.inla.md_conti_func_d.all_intercept_pdfd$lower)
x_max = max(domin_data.inla.md_conti_func_d.all_intercept_pdfd$upper)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin_pd_fd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=domin_data.inla.md_conti_func_d.all_intercept_pdfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=domin_data.inla.md_conti_func_d.all_intercept_pdfd,
               mapping = aes(x=mean,y=rowname, color = significant),
               size = 3.5)+
    geom_errorbar(data=domin_data.inla.md_conti_func_d.all_intercept_pdfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper, color = significant), width = 0.12,
                  linewidth = 1)+
    scale_color_manual(values = c('grey',
                                  'black'),
                       name = ' ')+
    ggnewscale::new_scale_color()+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD_all','MPD_all'),
                     label=c(expression(MFD), expression(MPD)))+
    #scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
       #                         ifelse(x_max<0,0.15,x_max+0.15)
    #))+
    scale_x_continuous(limits=c(-0.83,
                                0.2
    ), breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD_all" = "MPD_all",
                                  "MFD_all" = "MFD_all"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[4],
                                 colors_4d[3]),
                      labels = c("MPD_all" = "MPD_all",
                                 "MFD_all" = "MFD_all"),
                      name = ' ')+
    labs(x = 'Standardized effects', y = ' ')+
    guides(color="none")+
    theme_regular_1() 
)



### Effect size plot for abundance weighted mean differences
# point + effect size
domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd = 
  domin_data.inla.md.a_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MRFD.ab', 'MND.ab'))
x_min = min(domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd$lower)
x_max = max(domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd$upper)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin_nd_rfd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND.ab" = "MND.ab",
                                  "MRFD.ab" = "MRFD.ab"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[1],
                                 colors_4d[2]),
                      labels = c("MND.ab" = "MND.ab",
                                 "MRFD.ab" = "MRFD.ab"),
                      name = ' ')+
    geom_point(data=domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
               mapping = aes(x=mean,y=rowname, color = significant),
               size = 3.5)+
    geom_errorbar(data=domin_data.inla.md.a_conti_func_d.all_intercept_ndrfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper, color = significant),
                  width = 0.12, linewidth = 1)+
    scale_color_manual(values = c('black',
                                  'grey'),
                       name = ' ')+
    ggnewscale::new_scale_color()+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD.ab', 'MND.ab'),
                     label =c(c(expression(RFD['ab']),
                                expression(ND['ab']))))+
    #scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
        #                        ifelse(x_max<0,0.15,x_max+0.15)))+
    scale_x_continuous(limits=c(-0.2,
                                0.6))+
    labs(x = 'Standardized effects', y = ' ')+
    guides(color="none")+
    theme_regular_1()
)

domin_data.inla.md.a_conti_func_d.all_intercept_pdfd = 
  domin_data.inla.md.a_conti_func_d.all_intercept %>% 
  filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all'))
x_min = min(domin_data.inla.md.a_conti_func_d.all_intercept_pdfd$lower)
x_max = max(domin_data.inla.md.a_conti_func_d.all_intercept_pdfd$upper)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
(domin_pd_fd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=domin_data.inla.md.a_conti_func_d.all_intercept_pdfd,
             mapping = aes(x=mean,y=rowname,fill=rowname), 
             width = 0.4)+
    geom_point(data=domin_data.inla.md.a_conti_func_d.all_intercept_pdfd,
               mapping = aes(x=mean,y=rowname), color = 'grey',
               size = 3.5)+
    geom_errorbar(data=domin_data.inla.md.a_conti_func_d.all_intercept_pdfd,
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), color = 'grey',
                  width = 0.12, linewidth = 1)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all'),
                     label=c(expression(MFD['ab']), expression(MPD['ab'])))+
    #scale_x_continuous(limits=c(ifelse(x_min>0,-0.15,x_min-0.15),
        #                        ifelse(x_max<0,0.15,x_max+0.15)
   # ))+
    scale_x_continuous(limits=c(-0.83,
                                0.2
    ), breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD.ab_all" = "MPD.ab_all",
                                  "MFD.ab_all" = "MFD.ab_all"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[4],
                                 colors_4d[3]),
                      labels = c("MPD.ab_all" = "MPD.ab_all",
                                 "MFD.ab_all" = "MFD.ab_all"),
                      name = ' ')+
    labs(x = 'Standardized effects', y = ' ')+
    guides(color="none")+
    theme_regular_1() 
)




#### Fig.1：Demonstration of how Weighted Mean PD/FD affects successful colonization and dominance of species, including presentation of raw data and single & multi-regression results ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height =  0.47  # slightly smaller for more gap

# Slight vertical gap
top_row_y = 0.49
bottom_row_y = 0.01

# Create plot
Fig.1.all = ggdraw() +
  # Top row
  draw_plot(estab.mpd.a_single.logistic, x = 0, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mconti_func_d.a_single.logistic, x = plot_width, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(estab_pd_fd.a.all.varied.intercept.plot, x = 2 * plot_width,
            y = top_row_y, width = plot_width, height = plot_height) +
  
  # Bottom row
  draw_plot(domin_pd_fd.a.all.varied.intercept.plot, x = 2 * plot_width,
            y = bottom_row_y, width = plot_width, height = plot_height) +
  draw_plot(domin.mconti_func_d.a_single.logistic, x = plot_width, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mpd.a_single.logistic, x = 0, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  
  # Slightly higher labels
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, plot_width + 0.01, 2 * plot_width + 0.01,
          0.01, plot_width + 0.01, 2 * plot_width + 0.01),
    y = c(0.995, 0.995, 0.995, 0.515, 0.515, 0.515),  # slight bump upward
    hjust = 0, vjust = 1.1, size = 18, color = 'black'
  )


Fig.1.all
emf('results/figures/Fig.1.all.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.1.all
dev.off() #turn off device and finalize file




#### Fig.2：Demonstration of how Weighted Mean ND/RFD affects successful colonization and dominance of species, including presentation of raw data and single & multi-regression results ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height =  0.47  # slightly smaller for more gap

# Slight vertical gap
top_row_y = 0.49
bottom_row_y = 0.01

# Create plot
Fig.2.all = ggdraw() +
  # Top row
  draw_plot(estab.mlgfd.a_single.logistic, x = plot_width, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mnd.a_single.logistic, x = 0, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(estab_nd_rfd.a.all.varied.intercept.plot, x = 2 * plot_width,
            y = top_row_y, width = plot_width, height = plot_height) +
  
  # Bottom row
  draw_plot(domin_nd_rfd.a.all.varied.intercept.plot, x = 2 * plot_width,
            y = bottom_row_y, width = plot_width, height = plot_height) +
  draw_plot(domin.mlgfd.a_single.logistic, x = plot_width, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mnd.a_single.logistic, x = 0, y = bottom_row_y,
            width = plot_width, height = plot_height) +

  # Slightly higher labels
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, plot_width + 0.01, 2 * plot_width + 0.01,
          0.01, plot_width + 0.01, 2 * plot_width + 0.01),
    y = c(0.995, 0.995, 0.995, 0.515, 0.515, 0.515),  # slight bump upward
    hjust = 0, vjust = 1.1, size = 18, color = 'black'
  )


Fig.2.all
emf('results/figures/Fig.2.all.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.2.all
dev.off() #turn off device and finalize file



#### Fig.S9：Demonstration of how Mean ND/RFD affects successful colonization and dominance of species ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height =  0.47  # slightly smaller for more gap

# Slight vertical gap
top_row_y = 0.49
bottom_row_y = 0.01

# Create plot
Fig.S9.all = ggdraw() +
  # Top row
  draw_plot(estab.mnd_single.logistic, x = 0, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mlgfd_single.logistic, x = plot_width, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(estab_nd_rfd.all.varied.intercept.plot, x = 2 * plot_width,
            y = top_row_y, width = plot_width, height = plot_height) +
  # Bottom row
  draw_plot(domin.mnd_single.logistic, x = 0, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(domin_nd_rfd.all.varied.intercept.plot, x = 2 * plot_width,
            y = bottom_row_y, width = plot_width, height = plot_height) +
  draw_plot(domin.mlgfd_single.logistic, x = plot_width, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  # Slightly higher labels
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, plot_width + 0.01, 2 * plot_width + 0.01,
          0.01, plot_width + 0.01, 2 * plot_width + 0.01),
    y = c(0.995, 0.995, 0.995, 0.515, 0.515, 0.515),  # slight bump upward
    hjust = 0, vjust = 1.1, size = 18, color = 'black'
  )


Fig.S9.all
emf('results/figures/Fig.S9.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S9.all
dev.off() #turn off device and finalize file




#### Fig.S5：Demonstration of how Mean PD/FD affects successful colonization and dominance of species ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height =  0.47  # slightly smaller for more gap

# Slight vertical gap
top_row_y = 0.49
bottom_row_y = 0.01

# Create plot
Fig.S5.all = ggdraw() +
  # Top row
  draw_plot(estab.mpd_single.logistic, x = 0, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mconti_func_d_single.logistic, x = plot_width, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(estab_pd_fd.all.varied.intercept.plot, x = 2 * plot_width,
            y = top_row_y, width = plot_width, height = plot_height) +
  # Bottom row
  draw_plot(domin.mpd_single.logistic, x = 0, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(domin_pd_fd.all.varied.intercept.plot, x = 2 * plot_width,
            y = bottom_row_y, width = plot_width, height = plot_height) +
  draw_plot(domin.mconti_func_d_single.logistic, x = plot_width, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  # Slightly higher labels
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, plot_width + 0.01, 2 * plot_width + 0.01,
          0.01, plot_width + 0.01, 2 * plot_width + 0.01),
    y = c(0.995, 0.995, 0.995, 0.515, 0.515, 0.515),  # slight bump upward
    hjust = 0, vjust = 1.1, size = 18, color = 'black'
  )



#Fig.S5.all
emf('results/figures/Fig.S5.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S5.all
dev.off() #turn off device and finalize file


