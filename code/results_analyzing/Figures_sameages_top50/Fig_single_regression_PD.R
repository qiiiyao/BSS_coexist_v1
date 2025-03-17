############ Fast start for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for all species #########
rm(list = ls())
load("results/fit_results/BSS_exculde_trees_raw/plot_sameages_top50/inter_all_c_alltime.rdata")
load('results/fit_results/BSS_exculde_trees_raw/plot_sameages_top50/dat_all_alltime.rdata')

library(dplyr)
library(betareg)
library(stringr)
library(data.table)
library(ape)
library(reticulate)
inter_all_c = inter_all_c_alltime

## Functions for plot 
theme_for_coe_plot = function(x){
  ggplot2::theme_test() + 
    theme(legend.position = c(0.12, 0.15),
          legend.background = element_blank(),
          legend.key = element_blank(),
          text = element_text(family = 'Arial'),
          axis.title = element_text(face="bold"),
          axis.title.x = element_text(margin = margin(t = 0.5, r = 0,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.3,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0.5, l = 0.2, 
                               unit = "cm"))
}

theme_regular = function(x){
  ggplot2::theme_test() + 
    theme(text = element_text(family = 'Arial'),
          legend.background = element_blank(),
          legend.key = element_blank(),
          axis.title = element_text(face="bold"),
          axis.title.x = element_text(margin = margin(t = 0.5, r = 0,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0, r = 0.5,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0.5, l = 0.2, 
                               unit = "cm"))
}

# get the abs value of lgfd
inter_all_c$ablgfd = inter_all_c$lgfd
inter_all_c$ablgfd[inter_all_c$ablgfd < 0] = inter_all_c$ablgfd[inter_all_c$ablgfd < 0]*-1

#### Check the correlation between nd and fd
cor.test(inter_all_c$nd, inter_all_c$ablgfd) 
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
              inter_all_c$estab_i == 1&
              inter_all_c$domin_i == 0,]$stage_i_domin = 'estab.nodomin'
inter_all_c[inter_all_c$intro_i == 1&
              inter_all_c$estab_i == 1&
              inter_all_c$domin_i == 1,]$stage_i_domin = 'estab.domin'
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


######## Analyze invasion success ########
inter_all_c$ra_m_real_t_i = as.numeric(inter_all_c$ra_m_real_t_i)
inter_all_c$ra_m_real_t_j = as.numeric(inter_all_c$ra_m_real_t_j)
inter_all_c_l_2 = split(inter_all_c, inter_all_c$f_p)
dat_all_2 = split(dat_all_alltime, dat_all_alltime$f_p)

##### Plot level #####
dat_suc_plot = data.frame()

for (i in 1:length(inter_all_c_l_2)) {
  #i = 1
  trans_plot = inter_all_c_l_2[[i]]
  inv_suc = trans_plot %>% filter(stage_i != 'native' & stage_j == 'native')
  inv_stage = unique(inv_suc$stage_i)
  inv_suc_all = dat_all_2[[i]]
  
  for (j in 1:length(inv_stage)) {
    #j = 1
    inv_suc_stage = inv_suc %>% filter(stage_i == inv_stage[j])
    stage_sps = unique(inv_suc_stage$species_i)
    inv_suc_stage_all = inv_suc_all %>% filter(species_i %in% stage_sps)
    mnd = mean(inv_suc_stage$nd)
    mlgfd = mean(inv_suc_stage$lgfd)
    mablgfd = mean(inv_suc_stage$ablgfd)
    mpd = mean(inv_suc_stage$Phylo_dis)
    mfunc_d = mean(inv_suc_stage$Multi_traits)
    mconti_func_d = mean(inv_suc_stage$Multi_conti_traits)
    mgrowth = mean(inv_suc_stage$growth)
    mspan = mean(inv_suc_stage$span)
    mpollination = mean(inv_suc_stage$pollination)
    mdispersal = mean(inv_suc_stage$dispersal)
    mclonality = mean(inv_suc_stage$clonality)
    mheight = mean(inv_suc_stage$height)
    mldmc = mean(inv_suc_stage$ldmc)
    msla = mean(inv_suc_stage$sla)
    mseedmass = mean(inv_suc_stage$seedmass)
    
    mpd_all = mean(inv_suc_stage_all$Phylo_dis)
    mfunc_d_all = mean(inv_suc_stage_all$Multi_traits)
    mconti_func_d_all = mean(inv_suc_stage_all$Multi_conti_traits)
    mgrowth_all = mean(inv_suc_stage_all$growth)
    mspan_all = mean(inv_suc_stage_all$span)
    mpollination_all = mean(inv_suc_stage_all$pollination)
    mdispersal_all = mean(inv_suc_stage_all$dispersal)
    mclonality_all = mean(inv_suc_stage_all$clonality)
    mheight_all = mean(inv_suc_stage_all$height)
    mldmc_all = mean(inv_suc_stage_all$ldmc)
    msla_all = mean(inv_suc_stage_all$sla)
    mseedmass_all = mean(inv_suc_stage_all$seedmass)
    
    mnd.a = sum(inv_suc_stage$nd*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mlgfd.a = sum(inv_suc_stage$lgfd*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mablgfd.a = sum(inv_suc_stage$ablgfd*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mpd.a = sum(inv_suc_stage$Phylo_dis*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mfunc_d.a = sum(inv_suc_stage$Multi_traits*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mconti_func_d.a = sum(inv_suc_stage$Multi_conti_traits*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mgrowth.a = sum(inv_suc_stage$growth*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mspan.a = sum(inv_suc_stage$span*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mpollination.a = sum(inv_suc_stage$pollination*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mdispersal.a = sum(inv_suc_stage$dispersal*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mclonality.a = sum(inv_suc_stage$clonality*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mheight.a = sum(inv_suc_stage$height*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mldmc.a = sum(inv_suc_stage$ldmc*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    msla.a = sum(inv_suc_stage$sla*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    mseedmass.a = sum(inv_suc_stage$seedmass*inv_suc_stage$ra_m_real_t_j)/sum(inv_suc_stage$ra_m_real_t_j)
    
    mpd.a_all = sum(inv_suc_stage_all$Phylo_dis*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mfunc_d.a_all = sum(inv_suc_stage_all$Multi_traits*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mconti_func_d.a_all = sum(inv_suc_stage_all$Multi_conti_traits*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mgrowth.a_all = sum(inv_suc_stage_all$growth*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mspan.a_all = sum(inv_suc_stage_all$span*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mpollination.a_all = sum(inv_suc_stage_all$pollination*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mdispersal.a_all = sum(inv_suc_stage_all$dispersal*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mclonality.a_all = sum(inv_suc_stage_all$clonality*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mheight.a_all = sum(inv_suc_stage_all$height*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mldmc.a_all = sum(inv_suc_stage_all$ldmc*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    msla.a_all = sum(inv_suc_stage_all$sla*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    mseedmass.a_all = sum(inv_suc_stage_all$seedmass*inv_suc_stage_all$ra_m_real_t_j)/sum(inv_suc_stage_all$ra_m_real_t_j)
    
    inv_suc_stage_esp_l = split(inv_suc_stage, inv_suc_stage$species_i)
    mnnd = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$nd)}))
    mnlgfd = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$lgfd)}))
    mnablgfd = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$ablgfd)}))
    mntd = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$Phylo_dis)}))
    mnfunc_d = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$Multi_traits)}))
    mnconti_func_d = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$Multi_conti_traits)}))
    mngrowth = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$growth)}))
    mnspan = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$span)}))
    mnpollination = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$pollination)}))
    mndispersal = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$dispersal)}))
    mnclonality = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$clonality)}))
    mnheight = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$height)}))
    mnldmc = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$ldmc)}))
    mnsla = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$sla)}))
    mnseedmass = mean(sapply(inv_suc_stage_esp_l, function(x){y = min(x$seedmass)}))
    
    inv_suc_stage_all_esp_l = split(inv_suc_stage_all, inv_suc_stage_all$species_i)
    mntd_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$Phylo_dis)}))
    mnfunc_d_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$Multi_traits)}))
    mnconti_func_d_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$Multi_conti_traits)}))
    mngrowth_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$growth)}))
    mnspan_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$span)}))
    mnpollination_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$pollination)}))
    mndispersal_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$dispersal)}))
    mnclonality_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$clonality)}))
    mnheight_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$height)}))
    mnldmc_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$ldmc)}))
    mnsla_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$sla)}))
    mnseedmass_all = mean(sapply(inv_suc_stage_all_esp_l, function(x){y = min(x$seedmass)}))
    
    
    dat_suc_plot_1 = data.frame(f_p = unique(trans_plot$f_p), plot = unique(trans_plot$f_p), field = unique(trans_plot$field),
                              stage = inv_stage[j], 
                              estab = unique(inv_suc_stage$estab_i),
                              domin = unique(inv_suc_stage$domin_i),
                              mnd = mnd, mlgfd = mlgfd, mablgfd = mablgfd, mpd = mpd, mfunc_d = mfunc_d,
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
                              
                              mnd.a = mnd.a, mlgfd.a = mlgfd.a, mablgfd.a = mablgfd.a, mpd.a = mpd.a,
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
                              
                              mnnd = mnnd, mnlgfd = mnlgfd, mnablgfd = mnablgfd, mntd = mntd,
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
    dat_suc_plot = rbind(dat_suc_plot, dat_suc_plot_1)
  }
}

summary(dat_suc_plot)
dat_suc_plot$stage_level = NA
dat_suc_plot[dat_suc_plot$stage == 'introduce',]$stage_level = 1
dat_suc_plot[dat_suc_plot$stage == 'establish',]$stage_level = 2
dat_suc_plot[dat_suc_plot$stage == 'dominant',]$stage_level = 3
dat_suc_plot$stage_level = ordered(dat_suc_plot$stage_level)
str(dat_suc_plot)
save(dat_suc_plot,
     file = 'code/results_analyzing/analysing_sameages_top50_data/dat_suc_plot.rdata')


load('code/results_analyzing/analysing_sameages_top50_data/dat_suc_plot.rdata')
numcols = grep("^m",names(dat_suc_plot))
dat_suc_plots = dat_suc_plot
dat_suc_plots[,numcols] = scale(dat_suc_plots[,numcols])

#### estab ####
pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
#dat_suc_plots$species_1 = as.factor(dat_suc_plots$species)
estab_sp_names = unique(dat_suc_plots$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

### all functional trait distance
estab_model_mpd_all = pglmm(estab~mpd_all+#(1|species) + 
                                    (1|f_p) + (1|field), data = dat_suc_plots,
                                  family = "binomial", #cov_ranef = list(species = tree),
                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                              config = TRUE),
                                                       quantiles=c(0.025,0.5,0.975)),
                                  bayes = T)

estab_model_mpd.a_all = pglmm(estab~mpd.a_all+#(1|species) + 
                                      (1|f_p) + (1|field), data = dat_suc_plots,
                                    family = "binomial", #cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975)),
                                    bayes = T)

estab_model_mntd_all = pglmm(estab~mntd_all+#(1|species) + 
                                     (1|f_p) + (1|field), data = dat_suc_plots,
                                   family = "binomial", #cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T)
summary(estab_model_mpd_all)
summary(estab_model_mpd.a_all)
summary(estab_model_mntd_all)


###### establishment predictive curves for md ####
### get establishment predict data for mpd_all
lincombs.data.estab.mpd_all = data.frame(mpd_all=seq(min(dat_suc_sp$mpd_all),max(dat_suc_sp$mpd_all),length=100),
                                         mnd=mean(dat_suc_sp$mnd),
                                         mlgfd = mean(dat_suc_sp$mlgfd),
                                         mfunc_d_all = mean(dat_suc_sp$mfunc_d_all))

lincombs.matrix.estab.mpd_all=model.matrix(~mnd+mlgfd+mpd_all+mfunc_d_all,
                                           data=lincombs.data.estab.mpd_all)
lincombs.matrix.estab.mpd_all=as.data.frame(lincombs.matrix.estab.mpd_all)
lincombs.estab.mpd_all=inla.make.lincombs(lincombs.matrix.estab.mpd_all)

inla.model_lincombs.estab.mpd_all = pglmm(estab ~ mnd+mlgfd+mpd_all+mfunc_d_all+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_suc_sp,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.estab.mpd_all,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)

lincombs.posterior.estab.mpd_all = inla.model_lincombs.estab.mpd_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mpd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mpd_all$predicted.value=unlist(lapply(lincombs.posterior.estab.mpd_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mpd_all$lower=unlist(lapply(lincombs.posterior.estab.mpd_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mpd_all$upper=unlist(lapply(lincombs.posterior.estab.mpd_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mpd_all
save(lincombs.data.estab.mpd_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mpd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mpd_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mpd_all.rdata")
(estab.mpd_all.partial.logistic=ggplot(data=lincombs.data.estab.mpd_all,aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[4],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[4],size=1,
              linetype = 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=estab),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native phylogenetic difference', y='Establishment probability')+
    annotate(geom="text",x=c(200,200),y=c(0.80,0.70),
             label=c("italic(β)['I-N MPD'] == '-0.02'",
                     "'95%CI' == '[-0.03, 0.00]'"),
             parse=T,size=3.5)
)



### get establishment predict data for mfunc_d_all
lincombs.data.estab.mfunc_d_all = data.frame(mpd_all=mean(dat_suc_sp$mpd_all),
                                             mnd=mean(dat_suc_sp$mnd),
                                             mlgfd = mean(dat_suc_sp$mlgfd),
                                             mfunc_d_all = seq(min(dat_suc_sp$mfunc_d_all),
                                                               max(dat_suc_sp$mfunc_d_all),length=100))

lincombs.matrix.estab.mfunc_d_all=model.matrix(~mnd+mlgfd+mpd_all+mfunc_d_all,
                                               data=lincombs.data.estab.mfunc_d_all)
lincombs.matrix.estab.mfunc_d_all=as.data.frame(lincombs.matrix.estab.mfunc_d_all)
lincombs.estab.mfunc_d_all=inla.make.lincombs(lincombs.matrix.estab.mfunc_d_all)

inla.model_lincombs.estab.mfunc_d_all = pglmm(estab ~ mnd+mlgfd+mpd_all+mfunc_d_all+(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "binomial", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.estab.mfunc_d_all,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)

lincombs.posterior.estab.mfunc_d_all = inla.model_lincombs.estab.mfunc_d_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mfunc_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mfunc_d_all$predicted.value=unlist(lapply(lincombs.posterior.estab.mfunc_d_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mfunc_d_all$lower=unlist(lapply(lincombs.posterior.estab.mfunc_d_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mfunc_d_all$upper=unlist(lapply(lincombs.posterior.estab.mfunc_d_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mfunc_d_all
save(lincombs.data.estab.mfunc_d_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfunc_d_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mfunc_d_all.rdata")
(estab.mfunc_d_all.partial.logistic=ggplot(data=lincombs.data.estab.mfunc_d_all,
                                           aes(x=mfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[1],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[1],size=1,
              linetype= 1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d_all, y=estab),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native functional difference', y='  ')+
    annotate(geom="text",x=c(0.3,0.3),y=c(0.80,0.70),
             label=c("italic(β)['I-N MFD'] == 8.45",
                     "'95%CI' == '[3.93, 12.94]'"),
             parse=T,size=3.5)
)



###### establishment predictive curves for mean nearest difference ####
### get establishment predict data for mntd_all
lincombs.data.estab.mntd_all = data.frame(mntd_all=seq(min(dat_suc_sp$mntd_all),max(dat_suc_sp$mntd_all),length=100),
                                          mnnd=mean(dat_suc_sp$mnnd),
                                          mnlgfd = mean(dat_suc_sp$mnlgfd),
                                          mnfunc_d_all = mean(dat_suc_sp$mnfunc_d_all))

lincombs.matrix.estab.mntd_all=model.matrix(~mnnd+mnlgfd+mntd_all+mnfunc_d_all,
                                            data=lincombs.data.estab.mntd_all)
lincombs.matrix.estab.mntd_all=as.data.frame(lincombs.matrix.estab.mntd_all)
lincombs.estab.mntd_all=inla.make.lincombs(lincombs.matrix.estab.mntd_all)

inla.model_lincombs.estab.mntd_all = pglmm(estab ~ mnnd+mnlgfd+mntd_all+mnfunc_d_all+(1|species) + 
                                             (1|f_p) + (1|field), data = dat_suc_sp,
                                           family = "binomial", cov_ranef = list(species = tree),
                                           bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                       config = TRUE),
                                                                quantiles=c(0.025,0.5,0.975),
                                                                lincomb=lincombs.estab.mntd_all,
                                                                control.predictor=list(compute=T)),
                                           bayes = T)

lincombs.posterior.estab.mntd_all = inla.model_lincombs.estab.mntd_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mntd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(6)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mntd_all$predicted.value=unlist(lapply(lincombs.posterior.estab.mntd_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mntd_all$lower=unlist(lapply(lincombs.posterior.estab.mntd_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mntd_all$upper=unlist(lapply(lincombs.posterior.estab.mntd_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mntd_all
save(lincombs.data.estab.mntd_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mntd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mntd_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mntd_all.rdata")
(estab.mntd_all.partial.logistic=ggplot(data=lincombs.data.estab.mntd_all,aes(x=mntd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[4],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[4],size=1,
              linetype = 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mntd_all, y=estab),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native phylogenetic difference', y='Establishment probability')+
    annotate(geom="text",x=c(200,200),y=c(0.80,0.70),
             label=c("italic(β)['I-N MNTD'] == '-0.00'",
                     "'95%CI' == '[-0.00, 0.00]'"),
             parse=T,size=3.5)
)



### get establishment predict data for mnfunc_d_all
lincombs.data.estab.mnfunc_d_all = data.frame(mntd_all=mean(dat_suc_sp$mntd_all),
                                              mnnd=mean(dat_suc_sp$mnnd),
                                              mnlgfd = mean(dat_suc_sp$mnlgfd),
                                              mnfunc_d_all = seq(min(dat_suc_sp$mnfunc_d_all),
                                                                 max(dat_suc_sp$mnfunc_d_all),length=100))

lincombs.matrix.estab.mnfunc_d_all=model.matrix(~mnnd+mnlgfd+mntd_all+mnfunc_d_all,
                                                data=lincombs.data.estab.mnfunc_d_all)
lincombs.matrix.estab.mnfunc_d_all=as.data.frame(lincombs.matrix.estab.mnfunc_d_all)
lincombs.estab.mnfunc_d_all=inla.make.lincombs(lincombs.matrix.estab.mnfunc_d_all)

inla.model_lincombs.estab.mnfunc_d_all = pglmm(estab ~ mnnd+mnlgfd+mntd_all+mnfunc_d_all+(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_suc_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.estab.mnfunc_d_all,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.estab.mnfunc_d_all = inla.model_lincombs.estab.mnfunc_d_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnfunc_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnfunc_d_all$predicted.value=unlist(lapply(lincombs.posterior.estab.mnfunc_d_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnfunc_d_all$lower=unlist(lapply(lincombs.posterior.estab.mnfunc_d_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnfunc_d_all$upper=unlist(lapply(lincombs.posterior.estab.mnfunc_d_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnfunc_d_all
save(lincombs.data.estab.mnfunc_d_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mnfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnfunc_d_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mnfunc_d_all.rdata")
(estab.mnfunc_d_all.partial.logistic=ggplot(data=lincombs.data.estab.mnfunc_d_all,
                                            aes(x=mnfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[1],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[1],size=1,
              linetype= 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnfunc_d_all, y=estab),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native functional difference', y='  ')+
    annotate(geom="text",x=c(0.15,0.15),y=c(0.80,0.70),
             label=c("italic(β)['I-N MNFD'] == 0.97",
                     "'95%CI' == '[-1.24, 3.19]'"),
             parse=T,size=3.5)
)


###### establishment predictive curves for md.ab ####


### get establishment predict data for mpd.a_all
lincombs.data.estab.mpd.a_all = data.frame(mpd.a_all=seq(min(dat_suc_sp$mpd.a_all),max(dat_suc_sp$mpd.a_all),length=100),
                                           mnd.a=mean(dat_suc_sp$mnd.a),
                                           mlgfd.a = mean(dat_suc_sp$mlgfd.a),
                                           mfunc_d.a_all = mean(dat_suc_sp$mfunc_d.a_all))

lincombs.matrix.estab.mpd.a_all=model.matrix(~mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all,
                                             data=lincombs.data.estab.mpd.a_all)
lincombs.matrix.estab.mpd.a_all=as.data.frame(lincombs.matrix.estab.mpd.a_all)
lincombs.estab.mpd.a_all=inla.make.lincombs(lincombs.matrix.estab.mpd.a_all)

inla.model_lincombs.estab.mpd.a_all = pglmm(estab ~ mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "binomial", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.estab.mpd.a_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

lincombs.posterior.estab.mpd.a_all = inla.model_lincombs.estab.mpd.a_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mpd.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mpd.a_all$predicted.value=unlist(lapply(lincombs.posterior.estab.mpd.a_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mpd.a_all$lower=unlist(lapply(lincombs.posterior.estab.mpd.a_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mpd.a_all$upper=unlist(lapply(lincombs.posterior.estab.mpd.a_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mpd.a_all
save(lincombs.data.estab.mpd.a_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mpd.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mpd.a_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mpd.a_all.rdata")
(estab.mpd.a_all.partial.logistic=ggplot(data=lincombs.data.estab.mpd.a_all,aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[4],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[4],size=1,
              linetype = 3)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd.a_all, y=estab),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",x=c(200,200),y=c(0.80,0.70),
             label=c("italic(β)['MPD'[ab]] == '0.00'",
                     "'95%CI' == '[-0.01, 0.00]'"),
             parse=T,size=3.5)+
    theme_regular()
)



### get establishment predict data for mfunc_d.a_all
lincombs.data.estab.mfunc_d.a_all = data.frame(mpd.a_all=mean(dat_suc_sp$mpd.a_all),
                                               mnd.a=mean(dat_suc_sp$mnd.a),
                                               mlgfd.a = mean(dat_suc_sp$mlgfd.a),
                                               mfunc_d.a_all = seq(min(dat_suc_sp$mfunc_d.a_all),
                                                                   max(dat_suc_sp$mfunc_d.a_all),length=100))

lincombs.matrix.estab.mfunc_d.a_all=model.matrix(~mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all,
                                                 data=lincombs.data.estab.mfunc_d.a_all)
lincombs.matrix.estab.mfunc_d.a_all=as.data.frame(lincombs.matrix.estab.mfunc_d.a_all)
lincombs.estab.mfunc_d.a_all=inla.make.lincombs(lincombs.matrix.estab.mfunc_d.a_all)

inla.model_lincombs.estab.mfunc_d.a_all = pglmm(estab ~ mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_suc_sp,
                                                family = "binomial", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.estab.mfunc_d.a_all,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

lincombs.posterior.estab.mfunc_d.a_all = inla.model_lincombs.estab.mfunc_d.a_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mfunc_d.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mfunc_d.a_all$predicted.value=unlist(lapply(lincombs.posterior.estab.mfunc_d.a_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mfunc_d.a_all$lower=unlist(lapply(lincombs.posterior.estab.mfunc_d.a_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mfunc_d.a_all$upper=unlist(lapply(lincombs.posterior.estab.mfunc_d.a_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mfunc_d.a_all
save(lincombs.data.estab.mfunc_d.a_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mfunc_d.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfunc_d.a_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mfunc_d.a_all.rdata")
(estab.mfunc_d.a_all.partial.logistic=ggplot(data=lincombs.data.estab.mfunc_d.a_all,
                                             aes(x=mfunc_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[1],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[1],size=1,
              linetype= 3)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d.a_all, y=estab),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=c(0.3,0.3),y=c(0.80,0.70),
             label=c("italic(β)['MFD'[ab]] == -1.93", "'95%CI' == '[-3.50, 0.35]'"),
             parse=T,size=3.5)+
    theme_regular()
)


#### domin ####
dat_dom_plot = dat_suc_plot %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_plots = dat_suc_plots %>% filter(stage %in% c('establish', 'dominant'))


### all functional trait distance
domin_model_mpd = pglmm(domin~mpd_all+#(1|species) + 
                                    (1|f_p) + (1|field), data = dat_dom_plots,
                                  family = "binomial", #cov_ranef = list(species = tree),
                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                              config = TRUE),
                                                       quantiles=c(0.025,0.5,0.975)),
                                  bayes = T)

domin_model_mpd.a = pglmm(domin~mpd.a_all+#(1|species) + 
                                      (1|f_p) + (1|field), data = dat_dom_plots,
                                    family = "binomial", #cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975)),
                                    bayes = T)

domin_model_mntd = pglmm(domin~mntd_all+#(1|species) + 
                                     (1|f_p) + (1|field), data = dat_dom_plots,
                                   family = "binomial", #cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T)
summary(domin_model_mpd)
summary(domin_model_mpd.a)
summary(domin_model_mntd)

###### dominance predictive curves for md ####
### get dominance predict data for mpd_all
lincombs.data.domin.mpd_all = data.frame(mpd_all=seq(min(dat_dom_sp$mpd_all),max(dat_dom_sp$mpd_all),length=100),
                                         mnd=mean(dat_dom_sp$mnd),
                                         mlgfd = mean(dat_dom_sp$mlgfd),
                                         mfunc_d_all = mean(dat_dom_sp$mfunc_d_all))

lincombs.matrix.domin.mpd_all=model.matrix(~mnd+mlgfd+mpd_all+mfunc_d_all,
                                           data=lincombs.data.domin.mpd_all)
lincombs.matrix.domin.mpd_all=as.data.frame(lincombs.matrix.domin.mpd_all)
lincombs.domin.mpd_all=inla.make.lincombs(lincombs.matrix.domin.mpd_all)

inla.model_lincombs.domin.mpd_all = pglmm(domin ~ mnd+mlgfd+mpd_all+mfunc_d_all+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_dom_sp,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.domin.mpd_all,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)

lincombs.posterior.domin.mpd_all = inla.model_lincombs.domin.mpd_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mpd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mpd_all$predicted.value=unlist(lapply(lincombs.posterior.domin.mpd_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mpd_all$lower=unlist(lapply(lincombs.posterior.domin.mpd_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mpd_all$upper=unlist(lapply(lincombs.posterior.domin.mpd_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mpd_all
save(lincombs.data.domin.mpd_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mpd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mpd_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mpd_all.rdata")
(domin.mpd_all.partial.logistic=ggplot(data=lincombs.data.domin.mpd_all,aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[4],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[4],size=1,
              linetype = 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mpd_all, y=domin),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native phylogenetic difference', y='Dominance probability')+
    annotate(geom="text",x=c(260,260),y=c(0.80,0.70),
             label=c("italic(β)['E-N MPD'] == '-0.02'",
                     "'95%CI' == '[-0.03, 0.00]'"),
             parse=T,size=3.5)
)



### get dominance predict data for mfunc_d_all
lincombs.data.domin.mfunc_d_all = data.frame(mpd_all=mean(dat_dom_sp$mpd_all),
                                             mnd=mean(dat_dom_sp$mnd),
                                             mlgfd = mean(dat_dom_sp$mlgfd),
                                             mfunc_d_all = seq(min(dat_dom_sp$mfunc_d_all),
                                                               max(dat_dom_sp$mfunc_d_all),length=100))

lincombs.matrix.domin.mfunc_d_all=model.matrix(~mnd+mlgfd+mpd_all+mfunc_d_all,
                                               data=lincombs.data.domin.mfunc_d_all)
lincombs.matrix.domin.mfunc_d_all=as.data.frame(lincombs.matrix.domin.mfunc_d_all)
lincombs.domin.mfunc_d_all=inla.make.lincombs(lincombs.matrix.domin.mfunc_d_all)

inla.model_lincombs.domin.mfunc_d_all = pglmm(domin ~ mnd+mlgfd+mpd_all+mfunc_d_all+(1|species) + 
                                                (1|f_p) + (1|field), data = dat_dom_sp,
                                              family = "binomial", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.domin.mfunc_d_all,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)

lincombs.posterior.domin.mfunc_d_all = inla.model_lincombs.domin.mfunc_d_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mfunc_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mfunc_d_all$predicted.value=unlist(lapply(lincombs.posterior.domin.mfunc_d_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mfunc_d_all$lower=unlist(lapply(lincombs.posterior.domin.mfunc_d_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mfunc_d_all$upper=unlist(lapply(lincombs.posterior.domin.mfunc_d_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mfunc_d_all
save(lincombs.data.domin.mfunc_d_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfunc_d_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mfunc_d_all.rdata")
(domin.mfunc_d_all.partial.logistic=ggplot(data=lincombs.data.domin.mfunc_d_all,
                                           aes(x=mfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[1],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[1],size=1,
              linetype= 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfunc_d_all, y=domin),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native functional difference', y='  ')+
    annotate(geom="text",x=c(0.3,0.3),y=c(0.80,0.70),
             label=c("italic(β)['E-N MFD'] == 7.40",
                     "'95%CI' == '[-0.45, 15.07]'"),
             parse=T,size=3.5)
)


###### dominance predictive curves for mean nearest difference ####

### get dominance predict data for mntd_all
lincombs.data.domin.mntd_all = data.frame(mntd_all=seq(min(dat_dom_sp$mntd_all),max(dat_dom_sp$mntd_all),length=100),
                                          mnnd=mean(dat_dom_sp$mnnd),
                                          mnlgfd = mean(dat_dom_sp$mnlgfd),
                                          mnfunc_d_all = mean(dat_dom_sp$mnfunc_d_all))

lincombs.matrix.domin.mntd_all=model.matrix(~mnnd+mnlgfd+mntd_all+mnfunc_d_all,
                                            data=lincombs.data.domin.mntd_all)
lincombs.matrix.domin.mntd_all=as.data.frame(lincombs.matrix.domin.mntd_all)
lincombs.domin.mntd_all=inla.make.lincombs(lincombs.matrix.domin.mntd_all)

inla.model_lincombs.domin.mntd_all = pglmm(domin ~ mnnd+mnlgfd+mntd_all+mnfunc_d_all+(1|species) + 
                                             (1|f_p) + (1|field), data = dat_dom_sp,
                                           family = "binomial", cov_ranef = list(species = tree),
                                           bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                       config = TRUE),
                                                                quantiles=c(0.025,0.5,0.975),
                                                                lincomb=lincombs.domin.mntd_all,
                                                                control.predictor=list(compute=T)),
                                           bayes = T)

lincombs.posterior.domin.mntd_all = inla.model_lincombs.domin.mntd_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mntd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(6)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mntd_all$predicted.value=unlist(lapply(lincombs.posterior.domin.mntd_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mntd_all$lower=unlist(lapply(lincombs.posterior.domin.mntd_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mntd_all$upper=unlist(lapply(lincombs.posterior.domin.mntd_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mntd_all
save(lincombs.data.domin.mntd_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mntd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mntd_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mntd_all.rdata")
(domin.mntd_all.partial.logistic=ggplot(data=lincombs.data.domin.mntd_all,aes(x=mntd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[4],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[4],size=1,
              linetype = 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mntd_all, y=domin),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native phylogenetic difference', y='Dominance probability')+
    annotate(geom="text",x=c(150,150),y=c(0.80,0.70),
             label=c("italic(β)['E-N MNTD'] == '-0.00'",
                     "'95%CI' == '[-0.00, 0.00]'"),
             parse=T,size=3.5)
)



### get dominance predict data for mnfunc_d_all
lincombs.data.domin.mnfunc_d_all = data.frame(mntd_all=mean(dat_dom_sp$mntd_all),
                                              mnnd=mean(dat_dom_sp$mnnd),
                                              mnlgfd = mean(dat_dom_sp$mnlgfd),
                                              mnfunc_d_all = seq(min(dat_dom_sp$mnfunc_d_all),
                                                                 max(dat_dom_sp$mnfunc_d_all),length=100))

lincombs.matrix.domin.mnfunc_d_all=model.matrix(~mnnd+mnlgfd+mntd_all+mnfunc_d_all,
                                                data=lincombs.data.domin.mnfunc_d_all)
lincombs.matrix.domin.mnfunc_d_all=as.data.frame(lincombs.matrix.domin.mnfunc_d_all)
lincombs.domin.mnfunc_d_all=inla.make.lincombs(lincombs.matrix.domin.mnfunc_d_all)

inla.model_lincombs.domin.mnfunc_d_all = pglmm(domin ~ mnnd+mnlgfd+mntd_all+mnfunc_d_all+(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_dom_sp,
                                               family = "binomial", cov_ranef = list(species = tree),
                                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                           config = TRUE),
                                                                    quantiles=c(0.025,0.5,0.975),
                                                                    lincomb=lincombs.domin.mnfunc_d_all,
                                                                    control.predictor=list(compute=T)),
                                               bayes = T)

lincombs.posterior.domin.mnfunc_d_all = inla.model_lincombs.domin.mnfunc_d_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnfunc_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnfunc_d_all$predicted.value=unlist(lapply(lincombs.posterior.domin.mnfunc_d_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnfunc_d_all$lower=unlist(lapply(lincombs.posterior.domin.mnfunc_d_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnfunc_d_all$upper=unlist(lapply(lincombs.posterior.domin.mnfunc_d_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnfunc_d_all
save(lincombs.data.domin.mnfunc_d_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mnfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnfunc_d_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mnfunc_d_all.rdata")
(domin.mnfunc_d_all.partial.logistic=ggplot(data=lincombs.data.domin.mnfunc_d_all,
                                            aes(x=mnfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[1],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[1],size=1,
              linetype= 1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnfunc_d_all, y=domin),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native functional difference', y='  ')+
    annotate(geom="text",x=c(0.15,0.15),y=c(0.80,0.70),
             label=c("italic(β)['E-N MNFD'] == 4.17",
                     "'95%CI' == '[0.22, 8.11]'"),
             parse=T,size=3.5)
)



###### dominance predictive curves for md.ab ####

### get dominance predict data for mpd.a_all
lincombs.data.domin.mpd.a_all = data.frame(mpd.a_all=seq(min(dat_dom_sp$mpd.a_all),max(dat_dom_sp$mpd.a_all),length=100),
                                           mnd.a=mean(dat_dom_sp$mnd.a),
                                           mlgfd.a = mean(dat_dom_sp$mlgfd.a),
                                           mfunc_d.a_all = mean(dat_dom_sp$mfunc_d.a_all))

lincombs.matrix.domin.mpd.a_all=model.matrix(~mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all,
                                             data=lincombs.data.domin.mpd.a_all)
lincombs.matrix.domin.mpd.a_all=as.data.frame(lincombs.matrix.domin.mpd.a_all)
lincombs.domin.mpd.a_all=inla.make.lincombs(lincombs.matrix.domin.mpd.a_all)

inla.model_lincombs.domin.mpd.a_all = pglmm(domin ~ mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                              (1|f_p) + (1|field), data = dat_dom_sp,
                                            family = "binomial", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.domin.mpd.a_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

lincombs.posterior.domin.mpd.a_all = inla.model_lincombs.domin.mpd.a_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mpd.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mpd.a_all$predicted.value=unlist(lapply(lincombs.posterior.domin.mpd.a_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mpd.a_all$lower=unlist(lapply(lincombs.posterior.domin.mpd.a_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mpd.a_all$upper=unlist(lapply(lincombs.posterior.domin.mpd.a_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mpd.a_all
save(lincombs.data.domin.mpd.a_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mpd.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mpd.a_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mpd.a_all.rdata")
(domin.mpd.a_all.partial.logistic=ggplot(data=lincombs.data.domin.mpd.a_all,aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[4],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[4],size=1,
              linetype = 3)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mpd.a_all, y=domin),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Phylogenetic difference (PD)', y='Dominance probability')+
    annotate(geom="text",x=c(200,200),y=c(0.80,0.70),
             label=c("italic(β)['MPD'[ab]] == '-0.01'", "'95%CI' == '[-0.01, 0.00]'"),
             parse=T,size=3.5)+
    theme_regular()
)



### get dominance predict data for mfunc_d.a_all
lincombs.data.domin.mfunc_d.a_all = data.frame(mpd.a_all=mean(dat_dom_sp$mpd.a_all),
                                               mnd.a=mean(dat_dom_sp$mnd.a),
                                               mlgfd.a = mean(dat_dom_sp$mlgfd.a),
                                               mfunc_d.a_all = seq(min(dat_dom_sp$mfunc_d.a_all),
                                                                   max(dat_dom_sp$mfunc_d.a_all),length=100))

lincombs.matrix.domin.mfunc_d.a_all=model.matrix(~mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all,
                                                 data=lincombs.data.domin.mfunc_d.a_all)
lincombs.matrix.domin.mfunc_d.a_all=as.data.frame(lincombs.matrix.domin.mfunc_d.a_all)
lincombs.domin.mfunc_d.a_all=inla.make.lincombs(lincombs.matrix.domin.mfunc_d.a_all)

inla.model_lincombs.domin.mfunc_d.a_all = pglmm(domin ~ mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_dom_sp,
                                                family = "binomial", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.domin.mfunc_d.a_all,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

lincombs.posterior.domin.mfunc_d.a_all = inla.model_lincombs.domin.mfunc_d.a_all$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mfunc_d.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mfunc_d.a_all$predicted.value=unlist(lapply(lincombs.posterior.domin.mfunc_d.a_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mfunc_d.a_all$lower=unlist(lapply(lincombs.posterior.domin.mfunc_d.a_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mfunc_d.a_all$upper=unlist(lapply(lincombs.posterior.domin.mfunc_d.a_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mfunc_d.a_all
save(lincombs.data.domin.mfunc_d.a_all, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mfunc_d.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfunc_d.a_all
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mfunc_d.a_all.rdata")
(domin.mfunc_d.a_all.partial.logistic=ggplot(data=lincombs.data.domin.mfunc_d.a_all,
                                             aes(x=mfunc_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[1],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[1],size=1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfunc_d.a_all, y=domin),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Functional difference (FD)', y='  ')+
    annotate(geom="text",x=c(0.3,0.3),y=c(0.80,0.70),
             label=c("italic(β)['MFD'[ab]] == 5.89", "'95%CI' == '[3.28, 8.49]'"),
             parse=T,size=3.5)+
    theme_regular()
)


