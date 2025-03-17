############### Fast Start ####################
rm(list = ls())
load('results/fit_results/BSS_exculde_trees_raw/plot_sameages_top50/dat_all_alltime.rdata')
load("D:/R projects/BSS/code/data preparation/transformed data/fit_fp_same_ages_0.20_raw.RData")

library(dplyr)
library(betareg)
library(stringr)
library(data.table)
library(ape)
library(reticulate)

dat_all_l = lapply(fit_fp_same_ages_0.20_raw, function(x){
  alien_sp = c(str_split_1(x$stage_sp_name[[2]], ', '), 
               str_split_1(x$stage_sp_name[[3]], ', '))
  intro_sp = str_split_1(x$stage_sp_name[[2]], ', ')
  estab_sp = setdiff(str_split_1(x$stage_sp_name[[3]], ', '),
                     str_split_1(x$stage_sp_name[[4]], ', '))
  domin_sp = str_split_1(x$stage_sp_name[[4]], ', ')
  native_sp = setdiff(colnames(x$re_cover_ab)[10:ncol(x$re_cover_ab)],
                      alien_sp)
  f_p_1 = x$stage_sp_name[[1]]
  stage_infor = data.frame(sp = c(intro_sp, estab_sp, domin_sp, native_sp),
                           stage = c(rep('intro', length(intro_sp)),
                             rep('estab', length(estab_sp)),
                             rep('domin', length(domin_sp)),
                             rep('native', length(native_sp))))
  dat = data.frame(f_p = f_p_1, 
                   species_i = rep(alien_sp, each = length(native_sp)),
                   species_j = rep(native_sp, length(alien_sp)))
  dat_all = dat %>% left_join(stage_infor, by = join_by(species_i == sp)) %>% 
                      left_join(stage_infor, by = join_by(species_j == sp))
  colname_1 = gsub('\\.x', '_i', colnames(dat_all))          
  colname_2 = gsub('\\.y', '_j', colname_1)
  colnames(dat_all) = colname_2
  return(dat_all)
})

dat_all = as.data.frame(rbindlist(dat_all_l))
dat_all$sp_pair = paste(dat_all$species_i, dat_all$species_j, sep = '_')

## Add the gain or loss and relative abundance
load('code/Phylo_Func/fd_dat_all.rdata')
load('code/Phylo_Func/pd.rdata')
load('code/data preparation/transformed data/sp_racover_f1_2_mean_fp_alltime.rdata')

colnames(sp_racover_f1_2_mean_fp_alltime)
dat_all$f_p_species_i = paste(dat_all$f_p,
                              dat_all$species_i,
                              sep = '_')
dat_all$f_p_species_j = paste(dat_all$f_p,
                              dat_all$species_j,
                              sep = '_')

sp_racover_f1_2_mean_fp_alltime$f_p_species_i = sp_racover_f1_2_mean_fp_alltime$f_p_species
sp_racover_f1_2_mean_fp_alltime$f_p_species_j = sp_racover_f1_2_mean_fp_alltime$f_p_species_i

dat_all_test = dat_all %>% left_join(sp_racover_f1_2_mean_fp_alltime[,c(4, 6)], 
                                     by = 'f_p_species_i') %>%
  left_join(sp_racover_f1_2_mean_fp_alltime[,c(4, 7)], 
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
     file = 'results/fit_results/BSS_exculde_trees_raw/plot_sameages_top50/dat_all_alltime.rdata')
rm(list = ls()) ### free the memory


######## Analyse invasion success ########
load('results/fit_results/BSS_exculde_trees_raw/plot_sameages_top50/dat_all_alltime.rdata')
dat_all_alltime$estab_i = 0
dat_all_alltime$estab_j = 0
dat_all_alltime$domin_i = 0
dat_all_alltime$domin_j = 0
dat_all_alltime[dat_all_alltime$stage_i %in% c('estab', 'domin'),]$estab_i = 1
dat_all_alltime[dat_all_alltime$stage_i %in% c('domin'),]$domin_i = 1

dat_all_alltime$ra_m_real_t_i = as.numeric(dat_all_alltime$ra_m_real_t_i)
dat_all_alltime$ra_m_real_t_j = as.numeric(dat_all_alltime$ra_m_real_t_j)
dat_all_alltime_l_2 = split(dat_all_alltime, dat_all_alltime$f_p)

##### Species level #####
dat_suc_sp_nofit = data.frame()

for (i in 1:length(dat_all_alltime_l_2)) {
  #i = 1
  trans_plot = dat_all_alltime_l_2[[i]]
  inv_suc = trans_plot %>% filter(stage_i != 'native' & stage_j == 'native')
  inv_sp = unique(inv_suc$species_i)
  
  for (j in 1:length(inv_sp)) {
    #j = 1
    inv_suc_sp_all = inv_suc %>% filter(species_i == inv_sp[j])
    
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
 
    mpd.a_all = sum(inv_suc_sp_all$Phylo_dis*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mfunc_d.a_all = sum(inv_suc_sp_all$Multi_traits*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mconti_func_d.a_all = sum(inv_suc_sp_all$Multi_conti_traits*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mgrowth.a_all = sum(inv_suc_sp_all$growth*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mspan.a_all = sum(inv_suc_sp_all$span*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mpollination.a_all = sum(inv_suc_sp_all$pollination*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mdispersal.a_all = sum(inv_suc_sp_all$dispersal*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mclonality.a_all = sum(inv_suc_sp_all$clonality*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mheight.a_all = sum(inv_suc_sp_all$height*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mldmc.a_all = sum(inv_suc_sp_all$ldmc*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    msla.a_all = sum(inv_suc_sp_all$sla*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
    mseedmass.a_all = sum(inv_suc_sp_all$seedmass*inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)/sum(inv_suc_sp_all$ra_m_real_t_i*inv_suc_sp_all$ra_m_real_t_j)
  
    mntd_all = min(inv_suc_sp_all$Phylo_dis)
    mnfunc_d_all = min(inv_suc_sp_all$Multi_traits)
    mnconti_func_d_all = min(inv_suc_sp_all$Multi_conti_traits)
    mngrowth_all = min(inv_suc_sp_all$growth)
    mnspan_all = min(inv_suc_sp_all$span)
    mnpollination_all = min(inv_suc_sp_all$pollination)
    mndispersal_all = min(inv_suc_sp_all$dispersal)
    mnclonality_all = min(inv_suc_sp_all$clonality)
    mnheight_all = min(inv_suc_sp_all$height)
    mnldmc_all = min(inv_suc_sp_all$ldmc)
    mnsla_all = min(inv_suc_sp_all$sla)
    mnseedmass_all = min(inv_suc_sp_all$seedmass)
    
    
    dat_suc_sp_nofit_1 = data.frame(f_p = unique(trans_plot$f_p), plot = unique(trans_plot$f_p), field = str_split_1(unique(trans_plot$f_p),
                                                                                                               '_')[1],
                              species = inv_sp[j], 
                              stage = unique(inv_suc_sp_all$stage_i),
                              estab = unique(inv_suc_sp_all$estab_i),
                              domin = unique(inv_suc_sp_all$domin_i),
                              
                              mpd_all = mpd_all, mfunc_d_all = mfunc_d_all,
                              mconti_func_d_all = mconti_func_d_all,
                              mgrowth_all = mgrowth_all, mspan_all = mspan_all, 
                              mpollination_all = mpollination_all,
                              mdispersal_all = mdispersal_all, 
                              mclonality_all = mclonality_all,
                              mheight_all = mheight_all, mldmc_all = mldmc_all,
                              msla_all = msla_all, mseedmass_all = mseedmass_all,
                              
                              mpd.a_all = mpd.a_all,
                              mfunc_d.a_all = mfunc_d.a_all,
                              mconti_func_d.a_all = mconti_func_d.a_all,
                              mgrowth.a_all = mgrowth.a_all, mspan.a_all = mspan.a_all,
                              mpollination.a_all = mpollination.a_all,
                              mdispersal.a_all = mdispersal.a_all, mclonality.a_all = mclonality.a_all,
                              mheight.a_all = mheight.a_all, mldmc.a_all = mldmc.a_all,
                              msla.a_all = msla.a_all, mseedmass.a_all = mseedmass.a_all,
                             
                              mntd_all = mntd_all,
                              mnfunc_d_all = mnfunc_d_all, 
                              mnconti_func_d_all = mnconti_func_d_all, 
                              mngrowth_all = mngrowth_all, mnspan_all = mnspan_all, 
                              mnpollination_all = mnpollination_all,
                              mndispersal_all = mndispersal_all, 
                              mnclonality_all = mnclonality_all, 
                              mnheight_all = mnheight_all, mnldmc_all = mnldmc_all,
                              mnsla_all = mnsla_all, mnseedmass_all = mnseedmass_all)
    dat_suc_sp_nofit = rbind(dat_suc_sp_nofit, dat_suc_sp_nofit_1)
  }
}

summary(dat_suc_sp_nofit)
dat_suc_sp_nofit$stage_level = NA
dat_suc_sp_nofit[dat_suc_sp_nofit$stage == 'intro',]$stage_level = 1
dat_suc_sp_nofit[dat_suc_sp_nofit$stage == 'estab',]$stage_level = 2
dat_suc_sp_nofit[dat_suc_sp_nofit$stage == 'domin',]$stage_level = 3
dat_suc_sp_nofit$stage_level = ordered(dat_suc_sp_nofit$stage_level)
str(dat_suc_sp_nofit)
save(dat_suc_sp_nofit,
     file = 'code/results_analyzing/analysing_sameages_top50_data/dat_suc_sp_nofit.rdata')



#### domin ####
numcols = grep("^m",names(dat_suc_sp_nofit))
dat_suc_sp_nofits = dat_suc_sp_nofit
dat_suc_sp_nofits[,numcols] = scale(dat_suc_sp_nofits[,numcols])

dat_dom_sp_nofit = dat_suc_sp_nofit %>% filter(stage %in% c('estab', 'domin'))
dat_dom_sp_nofits = dat_suc_sp_nofits %>% filter(stage %in% c('estab', 'domin'))

## Check the co-linearity
car::vif(
  glmer(domin ~ mnd + mlgfd + mpd_all + mfunc_d_all + (1|species) + (1|f_p) + (1|field),
        family=binomial,data=dat_dom_sp_nofits)
)

car::vif(
  glmer(domin ~ mnd.a + mlgfd.a + mpd.a_all + mfunc_d.a_all + (1|species) + (1|f_p)+ (1|field),
        family=binomial,data=dat_dom_sp_nofits)
)

car::vif(
  glmer(domin ~ mnnd + mnlgfd + mntd_all + mnfunc_d_all + (1|species) + (1|f_p)+ (1|field),
        family=binomial,data=dat_dom_sp_nofits)
)

summary(glmer(estab ~ mpd_all + mfunc_d_all + 
               (1|f_p) + (1|field),
              family=binomial,data=dat_suc_sp_nofits))

summary(glmer(domin ~ mpd_all + mfunc_d_all + 
                (1|f_p) + (1|field),
              family=binomial,data=dat_dom_sp_nofits))


