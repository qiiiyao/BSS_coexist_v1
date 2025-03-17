# Load packages
rm(list = ls())

require(stringr)
require(dplyr)
require(data.table)

load("D:/R projects/BSS/code/data preparation/transformed data/fit_fp_top50_ages1_35_equal_interval.RData")
load("D:/R projects/BSS/code/data preparation/transformed data/re_cover_ab_f1_ages1_35_fp.rdata")
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
setwd("D:/R projects/BSS")
load("D:/R projects/BSS/results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_equal_interval_model_comparison/bh_partialb/inter_all_c_alltime_rhat_105_newfd.rdata")
load('results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_equal_interval_model_comparison/dat_all_alltime.rdata')

library(dplyr)
library(betareg)
library(stringr)
library(data.table)
library(ape)
library(colRoz)

## function for calculating coexistence probabilities
do.coexistence = function(stabilizing, ki.over.kj) {
  rho = 1 - stabilizing
  if (rho < ki.over.kj & (1/rho) > ki.over.kj) {
    outcome = data.frame(coexists = 1, priority = 0, i_exclusion_j = 0, j_exclusion_i = 0)
  } else if (rho > ki.over.kj & (1/rho) < ki.over.kj) {
    outcome = data.frame(coexists = 0, priority = 1, i_exclusion_j = 0, j_exclusion_i = 0)
  } else if (ki.over.kj > 1) {
    outcome = data.frame(coexists = 0, priority = 0, i_exclusion_j = 1, j_exclusion_i = 0)
  } else {
    outcome = data.frame(coexists = 0, priority = 0, i_exclusion_j = 0, j_exclusion_i = 1)
  }
  return(outcome)
}

#library(ochRe)
#library(reticulate)
inter_all_c = inter_all_c_alltime
rm(inter_all_c_alltime)
colnames(inter_all_c)

# Add the intrinsic growth rate and interaction coefficients 
# in the standard Lotka-Volterra model: dNi/Ni*dt = r-aii*Ni-aij*Nj
inter_all_c = inter_all_c %>% mutate(r_lv_i = (exp(r_trans_sp_i)-1)^(1/b_trans_sp_i),
                                     a_lv_ii = aii/((exp(r_trans_sp_i)-1)^(1/b_trans_sp_i)),
                                     r_lv_j = (exp(r_trans_sp_j)-1)^(1/b_trans_sp_j),
                                     a_lv_jj = ajj/((exp(r_trans_sp_j)-1)^(1/b_trans_sp_j)),
                                     a_lv_ij = aij/(((exp(r_trans_sp_i)-1)^(1/b_trans_sp_i)))) 

#outcomes = rbindlist((apply(inter_all_c, 1, function(x){
#outcomes = do.coexistence(stabilizing = as.numeric(x['nd']), ki.over.kj = as.numeric(x['fd']))     
#})))

#inter_all_c = cbind(inter_all_c, outcomes)

#save(inter_all_c,
#  file = 'results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_equal_interval_model_comparison/bh/inter_all_c_outcomes.rdata')

gre = "#55a868ff"
ora = "#dd8452ff"
yel = "#ccb974ff"
blu = "#4c72b0ff"
colors_4d = c(blu,ora,gre,yel)
colors_2d = c(ora, blu)

# get the abs value of lgfd
inter_all_c$lgfd= log10(inter_all_c$fd)
inter_all_c$ablgfd= abs(inter_all_c$lgfd)

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

### Draw the coexistence plot
library(lme4)
library(lmerTest)
library(dplyr)
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
library(scatterpie)
source('code/function/plot_func.R')
#loadfonts(device = "win")

#### Draw coexsit_probs_pie plots ####
inter_all_forexotics = inter_all_c %>%
  filter(stage_i != 'native', stage_j == 'native')
length(unique(inter_all_forexotics$species_i))
length(unique(inter_all_forexotics$species_j))

inter_all_forestab = inter_all_c %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))

inter_all_forestab_msd = inter_all_forestab %>% 
  group_by(stage_ij_estab) %>% summarise_at(vars(nd, lgfd),
                                            c(mean), na.rm = T)
inter_all_forestab_outcome = inter_all_forestab %>% 
  group_by(stage_ij_estab) %>% summarise_at(vars(coexists,
                                                 priority,
                                                 i_exclusion_j,
                                                 j_exclusion_i),
                                            c(probs = 
                                                function(x){
                                                  y = length(x[which(x == 1)])/length(x)}))
inter_all_forestab_outcomes = cbind(inter_all_forestab_msd,
                                    inter_all_forestab_outcome[,2:ncol(inter_all_forestab_outcome)])
colnames(inter_all_forestab_msd)[1:3] = c('stage', 'nd_mean', 'lgfd_mean')

BoundaryValues = data.frame(x = seq(0,0.99,l=100), 
                            y1 =  log10(1/ (1-seq(0,0.99,l=100))),
                            y2 = log10((1-seq(0,0.99,l=100))))      

BoundaryValues_P = data.frame(x = seq(0,-1,l=100), 
                              y1 =  log10( 1/ (1-seq(0,-1,l=100))),
                              y2 = log10((1-seq(0,-1,l=100))))      
OneLine = data.frame(x= c(rev(BoundaryValues$x),BoundaryValues$x),
                     y= c(rev(BoundaryValues$y2),BoundaryValues$y1 ))


#"#4c72b0ff" "#dd8452ff" "#55a868ff" "#ccb974ff"
Fig.1_compare_estab_original_pie = ggplot()+
  geom_rect(aes( xmin = -1, xmax = 1, ymax = 1, ymin =0), fill = '#ccb974ff',
            alpha = 0.4)+
  geom_rect(aes( xmin = -1, xmax = 1, ymax = 0, ymin =-1), fill = '#55a868ff',
            alpha = 0.4)+
  geom_path(data = OneLine,
            aes(x, y), col = 'red', size = 1.5)+
  geom_line(data = BoundaryValues,
            aes(x, y2), col = 'red', size = 1.5)+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_ribbon(data = BoundaryValues,
              aes(x = x, ymin = y2, ymax = y1), inherit.aes = FALSE, fill = 'grey75')+
  geom_ribbon(data = BoundaryValues_P,
              aes(x = x, ymin = y1, ymax = y2), inherit.aes = FALSE, fill = 'grey90')+
  geom_scatterpie(aes(x = nd, y = lgfd, group = stage_ij_estab), 
                  data = inter_all_forestab_outcomes,
                  pie_scale = 80,
                  cols= c('j_exclusion_i_probs', 'i_exclusion_j_probs',
                          'coexists_probs' ,'priority_probs'))+
  coord_fixed(xlim =c(-0.05, 0.07), ylim= c(-0.02, 0.1))+
  geom_point(aes(x = nd, y = lgfd, col = stage_ij_estab), size = 6,
             inter_all_forestab_outcomes)+ 
  geom_path(aes(nd, lgfd), arrow = arrow(length = unit(0.3, "cm"), ends = 'first'), size =1.5,
            data = inter_all_forestab_outcomes)+
  labs(x = 'Niche difference (ND)', 
       y = 'Relative fitness difference\n(RFD)') +
  annotate(geom="text", x=-0.03, y=0.097, family = 'Arial',
           label='Establishment', fontface = 'bold')+
  scale_colour_manual(values = c(colors_4d[2],colors_4d[1]),
                      labels = c("intro.estab_native" = "Successful",
                                 "intro.noestab_native" = "Failed"),
                      name = ' ')+ 
  scale_fill_manual(values = c('#55a868ff', '#ccb974ff',  'grey55','grey85'),
                    labels = c("i_exclusion_j_probs" = "Exotics Win",
                               "j_exclusion_i_probs" = "Natives Win",
                               "coexists_probs" = "Coexistence",
                               "priority_probs" = "Priority Effect"),
                    name = 'Outcome')+
  theme_for_coe_plot()+
  guides(
    fill = guide_legend(order = 1),  # Put shape legend first
    color = guide_legend(order = 2)   # Put color legend second
  ) +
  theme(legend.position = c(0.4, 0.7),
        legend.direction = "vertical", legend.box = "horizontal"
        #legend.box.background = element_rect(colour = 'black',
        #                                   linewidth = 2),
        #legend.title = element_text(margin = margin(t = 140)),  # create more white space above the legend titles
        #legend.box.margin = margin(b = -80)
  ) 
Fig.1_compare_estab_original_pie


inter_all_fordomin = inter_all_c %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))
inter_all_fordomin_msd = inter_all_fordomin %>% 
  group_by(stage_ij_domin) %>% summarise_at(vars(nd, lgfd),
                                            c(mean), na.rm = T)
inter_all_fordomin_outcome = inter_all_fordomin %>% 
  group_by(stage_ij_domin) %>% summarise_at(vars(coexists,
                                                 priority,
                                                 i_exclusion_j,
                                                 j_exclusion_i),
                                            c(probs = 
                                                function(x){
                                                  y = length(x[which(x == 1)])/length(x)}))
inter_all_fordomin_outcomes = cbind(inter_all_fordomin_msd,
                                    inter_all_fordomin_outcome[,2:ncol(inter_all_fordomin_outcome)])

Fig.1_compare_domin_original_pie = ggplot()+
  geom_rect(aes( xmin = -1, xmax = 1, ymax = 1, ymin =0), fill = '#ccb974ff',
            alpha = 0.4)+
  geom_rect(aes( xmin = -1, xmax = 1, ymax = 0, ymin =-1), fill = '#55a868ff',
            alpha = 0.4)+
  geom_path(data = OneLine,
            aes(x, y), col = 'red', size = 1.5)+
  geom_line(data = BoundaryValues,
            aes(x, y2), col = 'red', size = 1.5)+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_ribbon(data = BoundaryValues,
              aes(x = x, ymin = y2, ymax = y1), inherit.aes = FALSE, fill = 'grey75')+
  geom_ribbon(data = BoundaryValues_P,
              aes(x = x, ymin = y1, ymax = y2), inherit.aes = FALSE, fill = 'grey90')+
  geom_scatterpie(aes(x = nd, y = lgfd, group = stage_ij_domin), 
                  data = inter_all_fordomin_outcomes,
                  pie_scale = 60,
                  cols= c('j_exclusion_i_probs', 'i_exclusion_j_probs',
                          'coexists_probs' ,'priority_probs'))+
  #coord_fixed(xlim =c(-0.05, 0.35), ylim= c(0, 0.4))+
  coord_fixed(xlim =c(-0.05, 0.07), ylim= c(-0.02, 0.1))+
  geom_point(aes(x = nd, y = lgfd, col = stage_ij_domin), size = 6,
             inter_all_fordomin_outcomes)+ 
  geom_path(aes(nd, lgfd), arrow = arrow(length = unit(0.3, "cm"), ends = 'first'), size =1.5,
            data = inter_all_fordomin_outcomes)+
  labs(x = 'Niche difference (ND)', 
       y = '  ') +
  annotate(geom="text", x=-0.03, y=0.095, family = 'Arial',
           label='Dominance', fontface = 'bold')+
  scale_colour_manual(values = c(colors_4d[2],colors_4d[1]),
                      labels = c("estab.domin_native" = "Successful",
                                 "estab.nodomin_native" = "Failed"),
                      name = ' ')+ 
  scale_fill_manual(values = c('#55a868ff', '#ccb974ff',  'grey55','grey85'),
                    labels = c("i_exclusion_j_probs" = "Exotics Win",
                               "j_exclusion_i_probs" = "Natives Win",
                               "coexists_probs" = "Coexistence",
                               "priority_probs" = "Priority Effect"),
                    name = 'Outcome')+
  theme_for_coe_plot()+
  #theme_classic()+
  guides(
    fill = guide_legend(order = 1),  # Put shape legend first
    color = guide_legend(order = 2)   # Put color legend second
  ) +
  theme(#legend.position = c(0.4, 0.7),
    legend.position = 'none',
    legend.direction = "vertical", legend.box = "horizontal"
    #legend.box.background = element_rect(colour = 'black',
    #                                   linewidth = 2),
    #legend.title = element_text(margin = margin(t = 60)),  # create more white space above the legend titles
    #legend.box.margin = margin(b = -140)
  )
Fig.1_compare_domin_original_pie


#### Stage split two parts: estab noestab domin nodomin ####
##### just mean and sd of raw data #####
inter_all_forexotics = inter_all_c %>%
  filter(stage_i != 'native', stage_j == 'native')
length(unique(inter_all_forexotics$species_i))
length(unique(inter_all_forexotics$species_j))

inter_all_forestab = inter_all_c %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))

BoundaryValues=data.frame(x = seq(0,0.99,l=100), 
                           y1 =  log10(1/ (1-seq(0,0.99,l=100))),
                           y2 = log10((1-seq(0,0.99,l=100))))      


BoundaryValues_P=data.frame(x = seq(0,-0.99,l=100), 
                             y1 =   log10(1/ (1-seq(0,-0.99,l=100))),
                             y2 = log10((1-seq(0,-0.99,l=100))))      
OneLine= data.frame(x= c(rev(BoundaryValues$x),BoundaryValues$x),
                     y= c(rev(BoundaryValues$y2),BoundaryValues$y1 ))

length(unique(inter_all_forestab$sp_pair))
length(unique(inter_all_forestab$species_i))
length(unique(inter_all_forestab$species_j))

inter_all_forestab %>% 
  group_by(stage_ij_estab) %>% summarise_at(vars(nd, lgfd, fd),
                                            c(length))

inter_all_forestab_mse = inter_all_forestab %>% 
  group_by(stage_ij_estab) %>% summarise_at(vars(nd, lgfd, fd),
                                            c(mean, function(x){
                                              y = sd(x)
                                            }, function(x){
                                              y = sd(x)/sqrt(length(x))
                                            }))

colnames(inter_all_forestab_mse) = c('stage', 'nd_mean', 'lgfd_mean','fd_mean',
                                     'nd_sd', 'lgfd_sd', 'fd_sd',
                                     'nd_se', 'lgfd_se', 'fd_se')

windowsFonts(Arial=windowsFont("Arial"))

Fig.1_compare_estab_original = ggplot(NULL) +
  geom_rect(aes( xmin = -1, xmax = 1, ymax = 1, ymin =0), fill = 'white',
            alpha = 0.4)+
  geom_rect(aes( xmin = -1, xmax = 1, ymax = 0, ymin =-1), fill = 'white',
            alpha = 0.4)+
  geom_ribbon(data = BoundaryValues,
              aes(x = x, ymin = y2, ymax = y1), inherit.aes = FALSE, fill = 'grey75')+
  geom_ribbon(data = BoundaryValues_P,
              aes(x = x, ymin = y1, ymax = y2), inherit.aes = FALSE, fill = 'grey90')+
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.3) +
  #scale_y_log10(name = expression(paste("Fitness difference: ",kappa[j]/kappa[i])) )+
  #coord_cartesian(xlim=c(-0.5, 0.5),ylim=c(0.5, 1.5))+ # just control the lab, not the data
  geom_point(data = (inter_all_forestab %>% filter(nd == max(inter_all_forestab$nd))),
             aes(x = nd, y = fd),
             color = 'white',
             size = 1)+
  geom_pointrange(data = inter_all_forestab_mse,
                  mapping = aes(x = nd_mean,y = lgfd_mean,
                                ymin = lgfd_mean-lgfd_sd,
                                ymax = lgfd_mean+lgfd_sd,
                                color = stage,
                                shape = stage),
                  alpha = 1.2
                  ,fatten = 12,linewidth = 1.5
  )+
  geom_pointrange(data = inter_all_forestab_mse,
                  mapping = aes(x = nd_mean,y = lgfd_mean,
                                xmin = nd_mean-nd_sd,
                                xmax = nd_mean+nd_sd,
                                color = stage,
                                shape = stage),
                  alpha = 1.2
                  ,fatten = 12,linewidth = 1.5
  )+
  scale_color_manual(values = c(colors_2d[1],
                                colors_2d[2]),
                     labels = c("intro.estab_native" = "Successful",
                                "intro.noestab_native" = "Failed"),
                     name = ' ')+
  scale_shape_manual(values = c(16,1),
                     labels = c("intro.estab_native" = "Successful",
                                "intro.noestab_native" = "Failed"),
                     name = ' ')+
  scale_alpha_manual(values = c(0.3, 0.05),
                     labels = c("intro.estab_native" = "Successful",
                                "intro.noestab_native" = "Failed"),
                     name = ' ',
                     guide = NULL)+
  coord_cartesian(xlim=c(-0.06, 0.1),ylim=c(-0.3, 0.35))+ # just control the axis, not the data
  #coord_cartesian(xlim=c(-0.3, 0.35),ylim=c(-0.3, 0.35))+ # just control the axis, not the data
  annotate(geom="text", x=-0.032, y=0.35, family = 'Arial',
           label='Establishment', size = 6)+
  labs(x = 'Niche difference (ND)', 
       y = 'Relative fitness difference (RFD)') +
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(linewidth = 0.5),
                             family = 'Arial'),
         shape=guide_legend(title=NULL,
                            family = 'Arial')) +
  theme_for_coe_plot()

Fig.1_compare_estab_original

## Dominance
inter_all_fordomin = inter_all_c %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))
unique(inter_all_fordomin$stage_i)
unique(inter_all_fordomin$stage_j)

inter_all_fordomin %>% 
  group_by(stage_ij_domin) %>% summarise_at(vars(nd, lgfd, fd),
                                            c(length))

inter_all_fordomin_msd = inter_all_fordomin %>% 
  group_by(stage_ij_domin) %>% summarise_at(vars(nd, lgfd),
                                            c(mean, sd), na.rm = T)

colnames(inter_all_fordomin_msd) = c('stage', 'nd_mean',
                                     'lgfd_mean', 'nd_sd',
                                     'lgfd_sd')


Fig.1_compare_domin_original = ggplot(NULL) +
  geom_rect(aes( xmin = -1, xmax = 1, ymax = 1, ymin =0), fill = 'white',
            alpha = 0.4)+
  geom_rect(aes( xmin = -1, xmax = 1, ymax = 0, ymin =-1), fill = 'white',
            alpha = 0.4)+
  geom_ribbon(data = BoundaryValues,
              aes(x = x, ymin = y2, ymax = y1), inherit.aes = FALSE, fill = 'grey75')+
  geom_ribbon(data = BoundaryValues_P,
              aes(x = x, ymin = y1, ymax = y2), inherit.aes = FALSE, fill = 'grey90')+
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.3) +
  geom_point(data = (inter_all_fordomin %>% filter(nd == max(inter_all_fordomin$nd))),
             aes(x = nd, y = lgfd), color = 'white',
             alpha = 0.2, size = 1)+
  geom_pointrange(data = inter_all_fordomin_msd,
                  mapping = aes(x = nd_mean,y = lgfd_mean,
                                ymin = lgfd_mean-lgfd_sd,
                                ymax = lgfd_mean+lgfd_sd,
                                color = stage,
                                shape = stage)
                  ,fatten = 12,linewidth = 1.5
  )+
  geom_pointrange(data = inter_all_fordomin_msd,
                  mapping = aes(x = nd_mean, y = lgfd_mean,
                                xmin = nd_mean-nd_sd,
                                xmax = nd_mean+nd_sd,
                                color = stage,
                                shape = stage)
                  ,fatten = 12,linewidth = 1.5
  )+
  scale_color_manual(values = c(colors_2d[1],
                                colors_2d[2]),
                     labels = c("estab.domin_native" = "Successful",
                                "estab.nodomin_native" = "Failed"),
                     name = ' ')+
  scale_shape_manual(values = c(16,1),
                     labels = c("estab.domin_native" = "Successful",
                                "estab.nodomin_native" = "Failed"),
                     name = ' ')+
  coord_cartesian(xlim=c(-0.06, 0.1),ylim=c(-0.3, 0.35))+ # just control the axis, not the data
  #coord_cartesian(xlim=c(-0.3, 0.35),ylim=c(-0.3, 0.35))+ # just control the axis, not the data
  annotate(geom="text", x=-0.037, y=0.35, family = 'Arial',
           label='Dominance', size = 6)+
  labs(x = 'Niche difference (ND)', 
       y = '  ') +
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(linewidth = 0.3)),
         shape=guide_legend(title=NULL,
                            override.aes = list(linewidth = 0.3))) +
  theme_for_coe_plot() + 
  theme(legend.position = 'none')
Fig.1_compare_domin_original


##### including parameter uncertainty: just mean and sd of raw data #####
inter_all_forexotics = inter_all_c %>%
  filter(stage_i != 'native', stage_j == 'native')
length(unique(inter_all_forexotics$species_i))
length(unique(inter_all_forexotics$species_j))

inter_all_forestab = inter_all_c %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))

spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)

length(unique(inter_all_forestab$sp_pair))
length(unique(inter_all_forestab$species_i))
length(unique(inter_all_forestab$species_j))

inter_all_forestab_msd_2 = inter_all_forestab %>% 
  group_by(stage_ij_estab) %>% summarise_at(vars(nd_mean, lgfd_mean),
                                            c(mean, sd), na.rm = T)

colnames(inter_all_forestab_msd_2) = c('stage', 'nd_mean_2',
                                       'lgfd_mean_2', 'nd_sd',
                                       'lgfd_sd')

Fig.1_compare_estab_original_2 = ggplot(NULL) +
  geom_point(data = inter_all_forestab, aes(x = nd_mean, y = lgfd_mean,
                                            color = stage_ij_estab,
                                            shape = stage_ij_estab,
                                            alpha = stage_ij_estab),
             size = 1)+
  geom_pointrange(data = inter_all_forestab_msd_2,
                  mapping = aes(x = nd_mean_2,y = lgfd_mean_2,
                                ymin = lgfd_mean_2-lgfd_sd,
                                ymax = lgfd_mean_2+lgfd_sd,
                                color = stage,
                                shape = stage),
                  alpha = 1.2
                  ,fatten = 8,linewidth = 0.8
  )+
  geom_pointrange(data = inter_all_forestab_msd_2,
                  mapping = aes(x = nd_mean_2,y = lgfd_mean_2,
                                xmin = nd_mean_2-nd_sd,
                                xmax = nd_mean_2+nd_sd,
                                color = stage,
                                shape = stage),
                  alpha = 1.2
                  ,fatten = 8,linewidth = 0.8
  )+
  scale_color_manual(values = c(colors_2d[1],
                                colors_2d[2]),
                     labels = c("intro.estab_native" = "Successful",
                                "intro.noestab_native" = "Failed"),
                     name = ' ')+
  scale_shape_manual(values = c(16,1),
                     labels = c("intro.estab_native" = "Successful",
                                "intro.noestab_native" = "Failed"),
                     name = ' ')+
  scale_alpha_manual(values = c(0.05, 0.02),
                     labels = c("intro.estab_native" = "Successful",
                                "intro.noestab_native" = "Failed"),
                     name = ' ',
                     guide = NULL)+
  coord_cartesian(xlim=c(-2, 1),ylim=c(-1, 1))+ # just control the lab, not the data
  #scale_x_continuous(limits=c(-1, 1)) +
  #scale_y_continuous(limits=c(-0.9, 0.9)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.3) +
  annotate(geom="text", x=-0.6, y=0.8, family = 'Arial',
           label='Establishment')+
  labs(x = 'Niche difference (ND)', 
       y = 'Relative fitness difference\n(RFD)') +
  stat_function(fun=spss.f, colour="black", linewidth = 0.3)+
  stat_function(fun=spss.ff, colour="black", linewidth = 0.3)+
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(linewidth = 0.3),
                             family = 'Arial'),
         shape=guide_legend(title=NULL,
                            family = 'Arial')) +
  theme_for_coe_plot()
Fig.1_compare_estab_original_2

inter_all_fordomin = inter_all_c %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))

inter_all_fordomin_msd = inter_all_fordomin %>% 
  group_by(stage_ij_domin) %>% summarise_at(vars(nd, lgfd),
                                            c(mean, sd), na.rm = T)

colnames(inter_all_fordomin_msd) = c('stage', 'nd_mean',
                                     'lgfd_mean', 'nd_sd',
                                     'lgfd_sd')

Fig.1_compare_domin_original = ggplot(NULL) +
  geom_point(data = inter_all_fordomin, aes(x = nd, y = lgfd,
                                            color = stage_ij_domin,
                                            shape = stage_ij_domin),
             alpha = 0.05, size = 1)+
  geom_pointrange(data = inter_all_fordomin_msd,
                  mapping = aes(x = nd_mean,y = lgfd_mean,
                                ymin = lgfd_mean-lgfd_sd,
                                ymax = lgfd_mean+lgfd_sd,
                                color = stage,
                                shape = stage)
                  ,fatten = 8,linewidth = 0.8
  )+
  geom_pointrange(data = inter_all_fordomin_msd,
                  mapping = aes(x = nd_mean, y = lgfd_mean,
                                xmin = nd_mean-nd_sd,
                                xmax = nd_mean+nd_sd,
                                color = stage,
                                shape = stage)
                  ,fatten = 8,linewidth = 0.8
  )+
  scale_color_manual(values = c(colors_2d[1],
                                colors_2d[2]),
                     labels = c("estab.domin_native" = "Successful",
                                "estab.nodomin_native" = "Failed"),
                     name = ' ')+
  scale_shape_manual(values = c(16,1),
                     labels = c("estab.domin_native" = "Successful",
                                "estab.nodomin_native" = "Failed"),
                     name = ' ')+
  #scale_color_brewer(palette="Set2")+
  coord_cartesian(xlim=c(-1, 1),ylim=c(-0.9, 0.9))+ # just control the lab, not the data
  #scale_x_continuous(limits=c(-1, 1)) +
  #scale_y_continuous(limits=c(-0.9, 0.9)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.3) +
  annotate(geom="text", x=-0.6, y=0.8, family = 'Arial',
           label='Dominance')+
  labs(x = 'Niche difference (ND)', 
       y = '  ') +
  stat_function(fun=spss.f, colour="black", linewidth = 0.3)+
  stat_function(fun=spss.ff, colour="black", linewidth = 0.3)+
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(linewidth = 0.3)),
         shape=guide_legend(title=NULL,
                            override.aes = list(linewidth = 0.3))) +
  theme_for_coe_plot() 
#+ theme(legend.position = 'none')
Fig.1_compare_domin_original


##### Normal methods #####
unique(inter_all_c$stage_ij_estab)
inter_all_forestab = inter_all_c %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))


mod_nd_compare_1 = lmer(nd ~ stage_ij_estab+(1|field/plot),
                        data = inter_all_forestab, 
                        REML = F)
mod_nd_compare_lm = lm(nd ~ stage_ij_estab,
                       data = inter_all_forestab)
anova(mod_nd_compare_1, mod_nd_compare_lm)
mod_nd_compare_estab_1 = lmer(nd ~ stage_ij_estab+(1|field/plot),
                              data = inter_all_forestab, 
                              REML = T)
mod_nd_compare_estab_2 = lmer(nd ~ stage_ij_estab+(1|field/plot)+(1|sp_pair),
                              data = inter_all_forestab, 
                              REML = T)
mod_nd_compare_estab = lmer(nd ~ stage_ij_estab+(1|field/plot)+
                              (1|species_i) +(1|species_j) +
                              (1|sp_pair),
                            data = inter_all_forestab, 
                            REML = T)
anova(mod_nd_compare_estab_1, mod_nd_compare_estab) ### select mod_nd_compare_estab
anova(mod_nd_compare_estab_2, mod_nd_compare_estab) ### select mod_nd_compare_estab

mod_fd_compare_estab_1 = lmer(lgfd ~ stage_ij_estab+(1|field/plot),
                              data = inter_all_forestab,
                              REML = T)
mod_fd_compare_estab = lmer(lgfd ~ stage_ij_estab+(1|field/plot)+
                              (1|species_i) + (1|species_j) +
                              (1|sp_pair),
                            data = inter_all_forestab,
                            REML = T)
mod_fd_compare_estab_2 = lmer(lgfd ~ stage_ij_estab+(1|field/plot)+
                                (1|sp_pair),
                              data = inter_all_forestab,
                              REML = T)

anova(mod_fd_compare_estab_1, mod_fd_compare_estab) ### select mod_fd_compare_estab
anova(mod_fd_compare_estab_2, mod_fd_compare_estab) ### select mod_fd_compare_estab

summary(mod_nd_compare_estab)
summary(mod_fd_compare_estab)

pre_mod_nd_compare_estab = as.data.frame(ggpredict(mod_nd_compare_estab,
                                                   terms = 'stage_ij_estab',
                                                   type = 're'))

pre_mod_fd_compare_estab = as.data.frame(ggpredict(mod_fd_compare_estab,
                                                   terms = 'stage_ij_estab',
                                                   type = 're'))

pre_ndfd_compare_estab = data.frame(nd = pre_mod_nd_compare_estab$predicted,
                                    lgfd = pre_mod_fd_compare_estab$predicted,
                                    stage = pre_mod_nd_compare_estab$x,
                                    nd_sd = pre_mod_nd_compare_estab$std.error,
                                    nd_low = pre_mod_nd_compare_estab$conf.low,
                                    nd_high = pre_mod_nd_compare_estab$conf.high,
                                    fd_sd = pre_mod_fd_compare_estab$std.error,
                                    fd_low = pre_mod_fd_compare_estab$conf.low,
                                    fd_high = pre_mod_fd_compare_estab$conf.high)

spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)
# Define a collection of palettes to alter the default based on number of levels to encode
# Template function for creating densities grouped by a variable

Fig.1_compare_estab = ggplot(NULL) +
  geom_point(data = inter_all_forestab, aes(x = nd, y = lgfd,
                                            color = stage_ij_estab),
             alpha = 0.1, size = 30)+
  geom_point(data = pre_ndfd_compare_estab,
             mapping = aes(x = nd,y = lgfd),
             color = 'black',
             size = 50)+
  geom_pointrange(data = pre_ndfd_compare_estab,
                  mapping = aes(x = nd,y = lgfd, ymin = fd_low,
                                ymax = fd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  geom_pointrange(data = pre_ndfd_compare_estab,
                  mapping = aes(y = lgfd,x = nd, xmin = nd_low,
                                xmax = nd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  scale_color_discrete(type = c('#efcc01',
                                '#025302'
  ))+
  #scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(-1, 1.5)) +
  scale_y_continuous(limits=c(-0.8, 0.8)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=3) +
  xlab("Niche difference") +
  ylab(expression(paste(log[10],
                        "(Fitness difference)",
                        sep=' ')))+
  stat_function(fun=spss.f, colour="black", size = 5)+
  stat_function(fun=spss.ff, colour="black", size = 5)+
  #scale_x_continuous(breaks=c(-1, 1)) +
  annotate(geom="text", x=0.2, y=0.7, label=str_wrap("Invasion and exclude residents",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=0.25, y=-0.7, label=str_wrap("Residents repel invasion",
                                                       width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=0.7, label=str_wrap("Invasion and Coexistence",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=-0.7, label=str_wrap("Invasion and Coexistence",
                                                      width = 20), size = 50) +
  #annotate(geom="text", x=-0.7, y=0.1, label="Priority effect", size=20) +
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(size = 8))) +
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 0.5),
        panel.grid = element_blank(),
        legend.position = c(0.18, 0.9),
        legend.text = element_text(size=120),
        #legend.title = element_text(size=80),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000',face = 'bold'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000',face = 'bold'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000'))

Fig.1_compare_estab

#### domin nodomin
unique(inter_all_c$stage_ij_domin)
inter_all_fordomin = inter_all_c %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))

mod_nd_compare_domin_1 = lmer(nd ~ stage_ij_domin+(1|field/plot),
                              data = inter_all_fordomin, 
                              REML = F)
mod_nd_compare_domin_2 = lmer(nd ~ stage_ij_domin+(1|field/plot)+(1|sp_pair),
                              data = inter_all_fordomin, 
                              REML = T)
mod_nd_compare_domin = lmer(nd ~ stage_ij_domin+(1|field/plot)+
                              (1|species_i) + (1|species_j) +
                              (1|sp_pair),
                            data = inter_all_fordomin, 
                            REML = T)
anova(mod_nd_compare_domin, mod_nd_compare_domin_1) ## select mod_nd_compare_domin
anova(mod_nd_compare_domin_2, mod_nd_compare_domin) ## select mod_nd_compare_domin

mod_fd_compare_domin_1 = lmer(lgfd ~ stage_ij_domin+(1|field/plot),
                              data = inter_all_fordomin,
                              REML = T)
mod_fd_compare_domin_2 = lmer(lgfd ~ stage_ij_domin+(1|field/plot)+(1|sp_pair),
                              data = inter_all_fordomin,
                              REML = T)
mod_fd_compare_domin = lmer(lgfd ~ stage_ij_domin+(1|field/plot)+
                              (1|species_i) + (1|species_j) +
                              (1|sp_pair),
                            data = inter_all_fordomin,
                            REML = T)
anova(mod_fd_compare_domin,mod_fd_compare_domin_1) ## select mod_fd_compare_domin
summary(mod_nd_compare_domin)
summary(mod_fd_compare_domin)


pre_mod_nd_compare_domin = as.data.frame(ggpredict(mod_nd_compare_domin,
                                                   terms = 'stage_ij_domin',
                                                   type = 're'))

pre_mod_fd_compare_domin = as.data.frame(ggpredict(mod_fd_compare_domin,
                                                   terms = 'stage_ij_domin',
                                                   type = 're'))

pre_ndfd_compare_domin = data.frame(nd = pre_mod_nd_compare_domin$predicted,
                                    lgfd = pre_mod_fd_compare_domin$predicted,
                                    stage = pre_mod_nd_compare_domin$x,
                                    nd_sd = pre_mod_nd_compare_domin$std.error,
                                    nd_low = pre_mod_nd_compare_domin$conf.low,
                                    nd_high = pre_mod_nd_compare_domin$conf.high,
                                    fd_sd = pre_mod_fd_compare_domin$std.error,
                                    fd_low = pre_mod_fd_compare_domin$conf.low,
                                    fd_high = pre_mod_fd_compare_domin$conf.high)

spss.f = function(x) log10(1-x)
x = seq(-1, 1, 0.001)
spss.ff = function(x) -log10(1-x)

Fig.1_compare_domin = ggplot(NULL) +
  geom_point(data = inter_all_fordomin, aes(x = nd, y = lgfd,
                                            color = stage_ij_domin
  ), alpha = 0.1,
  size = 30)+
  geom_point(data = pre_ndfd_compare_domin,
             mapping = aes(x = nd,y = lgfd),
             color = 'black',
             size = 50)+
  geom_pointrange(data = pre_ndfd_compare_domin,
                  mapping = aes(x = nd,y = lgfd, ymin = fd_low,
                                ymax = fd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  geom_pointrange(data = pre_ndfd_compare_domin,
                  mapping = aes(y = lgfd,x = nd, xmin = nd_low,
                                xmax = nd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  scale_color_discrete(type = c('#f40407',
                                '#0266be'))+
  #scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(-1, 1.5)) +
  scale_y_continuous(limits=c(-0.8, 0.8)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=3) +
  xlab("Niche difference") +
  ylab(expression(paste(log[10],
                        "(Fitness difference)",
                        sep=' ')))+
  stat_function(fun=spss.f, colour="black", size = 5)+
  stat_function(fun=spss.ff, colour="black", size = 5)+
  #scale_x_continuous(breaks=c(-1, 1)) +
  annotate(geom="text", x=0.3, y=0.7, label=str_wrap("Invasion and exclude residents",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=0.3, y=-0.7, label=str_wrap("Residents repel invasion",
                                                      width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=0.7, label=str_wrap("Invasion and Coexistence",
                                                     width = 20), size = 50) +
  annotate(geom="text", x=1.2, y=-0.7, label=str_wrap("Invasion and Coexistence",
                                                      width = 20), size = 50) +
  #annotate(geom="text", x=-0.7, y=0.1, label="Priority effect", size=20) +
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(size = 8))) +
  theme(panel.background = element_rect(fill = 'white',color = 'black',linewidth = 0.5),
        panel.grid = element_blank(),
        legend.position = c(0.18, 0.9),
        legend.text = element_text(size = 120),
        plot.margin = margin(10,10,10,10),
        plot.background = element_blank(),
        text = element_text(size = 150),
        axis.ticks = element_line(linewidth = 3),
        axis.ticks.length = unit(-5,'lines'),
        axis.title.y = element_text(margin = margin(0,3,0,0),color = '#000000',face = 'bold'),
        axis.title.x = element_text(margin = margin(5,0,0,0),color = '#000000',face = 'bold'),
        axis.text.y = element_text(margin = margin(0,3,0,0),color = '#000000'),
        axis.text.x = element_text(margin = margin(5,0,0,0),color = '#000000'))
Fig.1_compare_domin

#ggsave('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.1_compare_domin.svg',
#       plot=Fig.1_compare_domin, 
#       width=105, height=80, dpi=600, units='cm',
#       limitsize=F)

gap = ggplot(NULL)+theme_void()
library(ggpubr)
Fig.1_split_stage = ggarrange(gap,Fig.1_compare_estab,gap,
                              Fig.1_compare_domin, nrow = 1, ncol = 4,
                              labels = c('','a)','','b)'), hjust = 1,
                              vjust = 1.5, widths = c(0.1, 1, 0.1, 1),
                              font.label = list(size = 150))

ggsave(plot = Fig.1_split_stage,
       'results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.1_split_stage.svg',
       width = 320,height = 100, dpi = 300, units = 'cm',
       limitsize = F)

Fig.1_split_stage_2 = ggarrange(gap,Fig.1_compare_estab,gap,
                                Fig.1_compare_domin, nrow = 4, ncol = 1,
                                labels = c('','a)','','b)'), hjust = -0.5,
                                vjust = 0, heights = c(0.1, 1, 0.1, 1),
                                font.label = list(size = 150))

ggsave(plot = Fig.1_split_stage_2,
       'results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.1_split_stage_2.svg',
       width = 165,height = 340, dpi = 300, units = 'cm',
       limitsize = F)

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






##### Calculating species level ND & RFD using N-F mapping #####
# some functions are adapted from Spaak_et_al (2021) Oikos
#library(reticulate)
#source('code/function/invasion_graph_main_functions.R')
#source_python("code/function/scheme_plot.py")
#source_python("code/function/numerical_NFD.py")

inter_all_c$ra_m_real_t_i = as.numeric(inter_all_c$ra_m_real_t_i)
inter_all_c$ra_m_real_t_j = as.numeric(inter_all_c$ra_m_real_t_j)
inter_all_c_l_2 = split(inter_all_c, inter_all_c$f_p)
dat_all_2 = split(dat_all_alltime, dat_all_alltime$f_p)

dat_equ_check = data.frame()

for (i in 1:length(inter_all_c_l_2)) {
  #i = 1
  trans_plot = inter_all_c_l_2[[i]]
  inv_suc = trans_plot %>% filter(stage_i != 'native' & stage_j == 'native')
  inv_sp = unique(inv_suc$species_i)
  native_sp = unique(inv_suc$species_j)
  inv_suc_all = dat_all_2[[i]]
  
  #for (j in 1:length(inv_sp)) {
    #j = 1
    inv_suc_sp = inv_suc %>% filter(species_i == inv_sp[j])
    inv_suc_sp_subplot = trans_plot %>% filter(species_i %in% c(inv_sp[j],
                                                              native_sp) &
                                                 species_j %in% c(inv_sp[j],
                                                                native_sp))
    
    inv_suc_sp_subplot = arrange(inv_suc_sp_subplot, inv_suc_sp_subplot$species_i)
    spe = unique(inv_suc_sp_subplot$species_i)
    aii = unique(inv_suc_sp_subplot$a_lv_ii)
    a = dcast(inv_suc_sp_subplot, species_i ~ species_j, value.var = 'a_lv_ij',
              fun.aggregate = sum, fill = 0)
    
    a = a[,-1]
    a = as.matrix(a)
    row.names(a) = colnames(a)
    focal_sp = which(colnames(a) == inv_sp[j])
    non_foc = which((colnames(a) %in% native_sp))
    
    for (j in 1:length(aii)) {
      a[colnames(a)[j], colnames(a)[j]] = aii[j]
    }
    a = a[non_foc, non_foc]
    non_zero_rows <- apply(a, 1, function(row) any(row != 0))
    non_zero_cols <- apply(a, 2, function(col) any(col != 0))
    
    xtemp=try(solve(a = a[non_zero_rows, non_zero_cols],
                       b = -1 * rep(1, length(non_foc))))
    if (min(xtemp) > 0 & class(xtemp) == 'numeric'){
      is_equ = 1
    } else if(class(xtemp) != 'numeric') {
      is_equ = 'error'
    } else if(min(xtemp) < 0) {
      is_equ = 0
    }

    dat_equ_check_1 = data.frame(f_p = unique(trans_plot$f_p),
                              plot = unique(trans_plot$f_p),
                              field = unique(trans_plot$field),
                              is_equ = is_equ)
                              #species = inv_sp[j],) 
                              
    dat_equ_check = rbind(dat_equ_check, dat_equ_check_1)
  #}
}
length(which(dat_equ_check$is_equ == 0))/480
# 0.9729167 native communities of the plots cannot coexist


#### Single regression results for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for all species ####
##### establishment predictive curves for md.ab #####
setwd("D:/R projects/BSS")
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
(estab.mnd.a_single.logistic=ggplot(data=lincombs.data.estab.mnd.a_single,aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd.a, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",x=min(dat_suc_sp$mnd.a)+(max(dat_suc_sp$mnd.a)-min(dat_suc_sp$mnd.a))*0.3,
             y=c(0.90,0.80),
             label=c("italic(β)['ND'[ab]] == '6.15'",
                     "'95%CI' == '[4.60, 7.70]'"),
             parse=T,size=4)+
    theme_regular()
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
(estab.mlgfd.a_single.logistic=ggplot(data=lincombs.data.estab.mlgfd.a_single,aes(x=mlgfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=2,
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
    annotate(geom="text",x=(max(dat_suc_sp$mlgfd.a)+min(dat_suc_sp$mlgfd.a))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['RFD'[ab]] == '2.93'",
                     "'95%CI' == '[2.49, 3.37]'"),
             parse=T,size=4)+
    theme_regular()
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
(estab.mpd.a_single.logistic=ggplot(data=lincombs.data.estab.mpd.a_single,aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=2,
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
    annotate(geom="text",x=(max(dat_suc_sp$mpd.a_all)+min(dat_suc_sp$mpd.a_all))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['MPD'[ab]] == '-0.0180'",
                     "'95%CI' == '[-0.0196, -0.0165]'"),
             parse=T,size=4)+
    theme_regular()
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
(estab.mconti_func_d.a_single.logistic=ggplot(data=lincombs.data.estab.mconti_func_d.a_single,
                                              aes(x=mconti_func_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=2,
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
             x=(max(dat_suc_sp$mconti_func_d.a_all)+min(dat_suc_sp$mconti_func_d.a_all))*0.5,
             y=c(0.9,0.80),
             label=c("italic(β)['MFD'[ab]] == '-0.6965'",
                     "'95%CI' == '[-0.8422, -0.5507]'"),
             parse=T,size=4)+
    theme_regular()
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
(domin.mnd.a_single.logistic=ggplot(data=lincombs.data.domin.mnd.a_single,aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_x_continuous(breaks = seq(-0.35, 0.7, 0.35))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd.a, y=domin),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Niche difference (ND)', y='Dominance probability')+
    annotate(geom="text",
             x=min(dat_dom_sp$mnd.a)+(max(dat_dom_sp$mnd.a)-min(dat_dom_sp$mnd.a))*0.3,
             y=c(0.90,0.80),
             label=c("italic(β)['ND'[ab]] == '5.08'",
                     "'95%CI' == '[2.35, 7.81]'"),
             parse=T,size=4)+
    theme_regular()
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
(domin.mlgfd.a_single.logistic=ggplot(data=lincombs.data.domin.mlgfd.a_single,aes(x=mlgfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mlgfd.a, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Relative fitness difference (RFD)', y='')+
    annotate(geom="text",
             x=(max(dat_dom_sp$mlgfd.a)+min(dat_dom_sp$mlgfd.a))*0.5,
             y=c(0.2,0.1),
             label=c("italic(β)['RFD'[ab]] == '1.32'",
                     "'95%CI' == '[0.70, 1.94]'"),
             parse=T,size=4)+
    theme_regular()
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
(domin.mpd.a_single.logistic=ggplot(data=lincombs.data.domin.mpd.a_single,aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=2,
              linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mpd.a_all, y=domin),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Phylogenetic difference (PD)', y='Dominance probability')+
    annotate(geom="text",x=(max(dat_dom_sp$mpd.a_all)+min(dat_dom_sp$mpd.a_all))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['MPD'[ab]] == '-0.0028'",
                     "'95%CI' == '[-0.0058, 0.0002]'"),
             parse=T,size=4)+
    theme_regular()
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
    labs(x='Functional difference (FD)', y='  ')+
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
(domin.mconti_func_d.a_single.logistic=ggplot(data=lincombs.data.domin.mconti_func_d.a_single,
                                              aes(x=mconti_func_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=2,
              linetype= 2)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mconti_func_d.a_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Functional difference (FD)', y='  ')+
    annotate(geom="text",x=(max(dat_dom_sp$mconti_func_d.a_all)+min(dat_dom_sp$mconti_func_d.a_all))*0.5,
             y=c(0.9,0.80),
             label=c("italic(β)['MFD'[ab]] == '-0.0776'",
                     "'95%CI' == '[-0.4201, 0.2635]'"),
             parse=T,size=4)+
    theme_regular()
)
#mconti_func_d.a_all -0.0776    -0.4201     0.2635


##### establishment predictive curves for md #####
setwd("D:/R projects/BSS")
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
(estab.mnd_single.logistic=ggplot(data=lincombs.data.estab.mnd_single,aes(x=mnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",
             x=min(dat_suc_sp$mnd)+(max(dat_suc_sp$mnd)-min(dat_suc_sp$mnd))*0.3,
             y=c(0.90,0.80),
             label=c("italic(β)['ND'] == '12.87'",
                     "'95%CI' == '[9.98, 15.76]'"),
             parse=T,size=3.5)+
    theme_regular()
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
(estab.mlgfd_single.logistic=ggplot(data=lincombs.data.estab.mlgfd_single,aes(x=mlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=2,
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
    annotate(geom="text",
             x=min(dat_suc_sp$mlgfd)+(max(dat_suc_sp$mlgfd)-min(dat_suc_sp$mlgfd))*0.7,
             y=c(0.30,0.20),
             label=c("italic(β)['RFD'] == '3.79'",
                     "'95%CI' == '[3.27, 4.32]'"),
             parse=T,size=3.5)+
    theme_regular()
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
(estab.mpd_single.logistic=ggplot(data=lincombs.data.estab.mpd_single,aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=2,
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
    annotate(geom="text",x=min(dat_suc_sp$mpd_all)+(max(dat_suc_sp$mpd_all)-min(dat_suc_sp$mpd_all))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['MPD'] == '-0.0261'",
                     "'95%CI' == '[-0.0284, -0.0238]'"),
             parse=T,size=3.5)+
    theme_regular()
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
(estab.mconti_func_d_single.logistic=ggplot(data=lincombs.data.estab.mconti_func_d_single,
                                            aes(x=mconti_func_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=2,
              linetype= 2)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d_all, y=estab),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=min(dat_suc_sp$mconti_func_d_all)+
               (max(dat_suc_sp$mconti_func_d_all)-min(dat_suc_sp$mconti_func_d_all))*0.5,
             y=c(0.9,0.80),
             label=c("italic(β)['MFD'] == '-0.1117'",
                     "'95%CI' == '[-0.3285, 0.1051]'"),
             parse=T,size=3.5)+
    theme_regular()
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
(domin.mnd_single.logistic=ggplot(data=lincombs.data.domin.mnd_single,aes(x=mnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd, y=domin),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Niche difference (ND)', y='Dominance probability')+
    annotate(geom="text",x=min(dat_dom_sp$mnd)+
               (max(dat_dom_sp$mnd)-min(dat_dom_sp$mnd))*0.3,
             y=c(0.90,0.80),
             label=c("italic(β)['ND'] == '11.65'",
                     "'95%CI' == '[6.57, 16.72]'"),
             parse=T,size=3.5)+
    theme_regular()
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
(domin.mlgfd_single.logistic=ggplot(data=lincombs.data.domin.mlgfd_single,aes(x=mlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mlgfd, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Relative fitness difference (RFD)', y='')+
    annotate(geom="text",x=min(dat_dom_sp$mlgfd)+
               (max(dat_dom_sp$mlgfd)-min(dat_dom_sp$mlgfd))*0.7,
             y=c(0.30,0.20),
             label=c("italic(β)['RFD'] == '1.88'",
                     "'95%CI' == '[1.10, 2.67]'"),
             parse=T,size=3.5)+
    theme_regular()
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
(domin.mpd_single.logistic=ggplot(data=lincombs.data.domin.mpd_single,aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=2,
              linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mpd_all, y=domin),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Phylogenetic difference (PD)', y='Dominance probability')+
    annotate(geom="text",x=min(dat_dom_sp$mpd_all)+
               (max(dat_dom_sp$mpd_all)-min(dat_dom_sp$mpd_all))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['MPD'] == '-0.0015'",
                     "'95%CI' == '[-0.0062, 0.0032]'"),
             parse=T,size=3.5)+
    theme_regular()
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
(domin.mconti_func_d_single.logistic=ggplot(data=lincombs.data.domin.mconti_func_d_single,
                                            aes(x=mconti_func_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=2, linetype = 2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_x_continuous(breaks = c(0.2, 0.3, 0.4))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mconti_func_d_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Functional difference (FD)', y='  ')+
    annotate(geom="text",x=min(dat_dom_sp$mconti_func_d_all)+
               (max(dat_dom_sp$mconti_func_d_all)-min(dat_dom_sp$mconti_func_d_all))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['MFD'] == '-0.2371'",
                     "'95%CI' == '[-0.7147, 0.2411]'"),
             parse=T,size=3.5)+
    theme_regular()
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

##### establishment predictive curves for mean nearest differences #####
setwd("D:/R projects/BSS")
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


### get establishment predict data for mnnd
lincombs.data.estab.mnnd_single = data.frame(mnnd=seq(min(dat_suc_sp$mnnd),
                                                      max(dat_suc_sp$mnnd),
                                                      length=100))

lincombs.matrix.estab.mnnd_single=model.matrix(~mnnd,
                                               data=lincombs.data.estab.mnnd_single)
lincombs.matrix.estab.mnnd_single=as.data.frame(lincombs.matrix.estab.mnnd_single)
lincombs.estab.mnnd_single=inla.make.lincombs(lincombs.matrix.estab.mnnd_single)

inla.model_lincombs.estab.mnnd_single = pglmm(estab ~ mnnd+#(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "binomial", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.estab.mnnd_single,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)

lincombs.posterior.estab.mnnd_single = inla.model_lincombs.estab.mnnd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnnd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnnd_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mnnd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnnd_single$lower=unlist(lapply(lincombs.posterior.estab.mnnd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnnd_single$upper=unlist(lapply(lincombs.posterior.estab.mnnd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnnd_single
save(lincombs.data.estab.mnnd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnnd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnnd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnnd_single.rdata")
(estab.mnnd_single.logistic=ggplot(data=lincombs.data.estab.mnnd_single,aes(x=mnnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnnd, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    scale_x_continuous(breaks = round(seq(-0.60, 0.2, length.out = 5), 1))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",x=min(dat_suc_sp$mnnd)+
               (max(dat_suc_sp$mnnd)-min(dat_suc_sp$mnnd))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['ND'[nearest]] == '7.14'",
                     "'95%CI' == '[4.81, 9.46]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mnnd         7.14       4.81       9.46



### get establishment predict data for mnlgfd
lincombs.data.estab.mnlgfd_single = data.frame(mnlgfd=seq(min(dat_suc_sp$mnlgfd),
                                                          max(dat_suc_sp$mnlgfd),
                                                          length=100))

lincombs.matrix.estab.mnlgfd_single=model.matrix(~mnlgfd,
                                                 data=lincombs.data.estab.mnlgfd_single)
lincombs.matrix.estab.mnlgfd_single=as.data.frame(lincombs.matrix.estab.mnlgfd_single)
lincombs.estab.mnlgfd_single=inla.make.lincombs(lincombs.matrix.estab.mnlgfd_single)

inla.model_lincombs.estab.mnlgfd_single = pglmm(estab ~ mnlgfd+#(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_suc_sp,
                                                family = "binomial", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.estab.mnlgfd_single,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

lincombs.posterior.estab.mnlgfd_single = inla.model_lincombs.estab.mnlgfd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnlgfd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnlgfd_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mnlgfd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnlgfd_single$lower=unlist(lapply(lincombs.posterior.estab.mnlgfd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnlgfd_single$upper=unlist(lapply(lincombs.posterior.estab.mnlgfd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnlgfd_single
save(lincombs.data.estab.mnlgfd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnlgfd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnlgfd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnlgfd_single.rdata")
(estab.mnlgfd_single.logistic=ggplot(data=lincombs.data.estab.mnlgfd_single,aes(x=mnlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnlgfd, y=estab),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y=' ')+
    annotate(geom="text",x=min(dat_suc_sp$mnlgfd)+
               (max(dat_suc_sp$mnlgfd)-min(dat_suc_sp$mnlgfd))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['RFD'[nearest]] == '1.06'",
                     "'95%CI' == '[0.79, 1.32]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mnlgfd       1.06       0.79       1.32


### get establishment predict data for mntd_all
lincombs.data.estab.mntd_single = data.frame(mntd_all=seq(min(dat_suc_sp$mntd_all),
                                                          max(dat_suc_sp$mntd_all),
                                                          length=100))

lincombs.matrix.estab.mntd_single=model.matrix(~mntd_all,
                                               data=lincombs.data.estab.mntd_single)
lincombs.matrix.estab.mntd_single=as.data.frame(lincombs.matrix.estab.mntd_single)
lincombs.estab.mntd_single=inla.make.lincombs(lincombs.matrix.estab.mntd_single)

inla.model_lincombs.estab.mntd_single = pglmm(estab ~ mntd_all+#(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "binomial", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.estab.mntd_single,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)

lincombs.posterior.estab.mntd_single = inla.model_lincombs.estab.mntd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mntd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mntd_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mntd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mntd_single$lower=unlist(lapply(lincombs.posterior.estab.mntd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mntd_single$upper=unlist(lapply(lincombs.posterior.estab.mntd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mntd_single
save(lincombs.data.estab.mntd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mntd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mntd_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mntd_single.rdata")
(estab.mntd_single.logistic=ggplot(data=lincombs.data.estab.mntd_single,aes(x=mntd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mntd_all, y=estab),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",x=min(dat_suc_sp$mntd_all)+
               (max(dat_suc_sp$mntd_all)-min(dat_suc_sp$mntd_all))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['NPD'] == '-0.0022'",
                     "'95%CI' == '[-0.0031, -0.0013]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mntd_all    -0.0022    -0.0031    -0.0013



### get establishment predict data for mnconti_func_d_all
lincombs.data.estab.mnconti_func_d_single = data.frame(mnconti_func_d_all = seq(min(dat_suc_sp$mnconti_func_d_all),
                                                                                max(dat_suc_sp$mnconti_func_d_all),
                                                                                length=100))

lincombs.matrix.estab.mnconti_func_d_single=model.matrix(~mnconti_func_d_all,
                                                         data=lincombs.data.estab.mnconti_func_d_single)
lincombs.matrix.estab.mnconti_func_d_single=as.data.frame(lincombs.matrix.estab.mnconti_func_d_single)
lincombs.estab.mnconti_func_d_single=inla.make.lincombs(lincombs.matrix.estab.mnconti_func_d_single)

inla.model_lincombs.estab.mnconti_func_d_single = pglmm(estab ~ mnconti_func_d_all+#(1|species) + 
                                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                                        family = "binomial", cov_ranef = list(species = tree),
                                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                    config = TRUE),
                                                                             quantiles=c(0.025,0.5,0.975),
                                                                             lincomb=lincombs.estab.mnconti_func_d_single,
                                                                             control.predictor=list(compute=T)),
                                                        bayes = T)

lincombs.posterior.estab.mnconti_func_d_single = inla.model_lincombs.estab.mnconti_func_d_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnconti_func_d_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnconti_func_d_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mnconti_func_d_single,
                                                                        function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnconti_func_d_single$lower=unlist(lapply(lincombs.posterior.estab.mnconti_func_d_single,
                                                              function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnconti_func_d_single$upper=unlist(lapply(lincombs.posterior.estab.mnconti_func_d_single,
                                                              function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnconti_func_d_single
save(lincombs.data.estab.mnconti_func_d_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnconti_func_d_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnconti_func_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnconti_func_d_single.rdata")
(estab.mnconti_func_d_single.logistic=ggplot(data=lincombs.data.estab.mnconti_func_d_single,
                                             aes(x=mnconti_func_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=2,
              linetype= 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnconti_func_d_all, y=estab),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=min(dat_suc_sp$mnconti_func_d_all)+
               (max(dat_suc_sp$mnconti_func_d_all)-min(dat_suc_sp$mnconti_func_d_all))*0.5,
             y=c(0.9,0.80),
             label=c("italic(β)['NFD'] == '1.3311'",
                     "'95%CI' == '[1.1207, 1.5416]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mnconti_func_d_all  1.3311     1.1207     1.5416


### get establishment predict data for mnfunc_d_all
lincombs.data.estab.mnfunc_d_single = data.frame(mnfunc_d_all = seq(min(dat_suc_sp$mnfunc_d_all),
                                                                    max(dat_suc_sp$mnfunc_d_all),
                                                                    length=100))

lincombs.matrix.estab.mnfunc_d_single=model.matrix(~mnfunc_d_all,
                                                   data=lincombs.data.estab.mnfunc_d_single)
lincombs.matrix.estab.mnfunc_d_single=as.data.frame(lincombs.matrix.estab.mnfunc_d_single)
lincombs.estab.mnfunc_d_single=inla.make.lincombs(lincombs.matrix.estab.mnfunc_d_single)

inla.model_lincombs.estab.mnfunc_d_single = pglmm(estab ~ mnfunc_d_all+#(1|species) + 
                                                    (1|f_p) + (1|field), data = dat_suc_sp,
                                                  family = "binomial", cov_ranef = list(species = tree),
                                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                              config = TRUE),
                                                                       quantiles=c(0.025,0.5,0.975),
                                                                       lincomb=lincombs.estab.mnfunc_d_single,
                                                                       control.predictor=list(compute=T)),
                                                  bayes = T)

lincombs.posterior.estab.mnfunc_d_single = inla.model_lincombs.estab.mnfunc_d_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnfunc_d_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnfunc_d_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mnfunc_d_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnfunc_d_single$lower=unlist(lapply(lincombs.posterior.estab.mnfunc_d_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnfunc_d_single$upper=unlist(lapply(lincombs.posterior.estab.mnfunc_d_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnfunc_d_single
save(lincombs.data.estab.mnfunc_d_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnfunc_d_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnfunc_d_single.rdata")
(estab.mnfunc_d_single.logistic=ggplot(data=lincombs.data.estab.mnfunc_d_single,
                                       aes(x=mnfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1,
              linetype= 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnfunc_d_all, y=estab),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=c(0.08,0.08),y=c(0.9,0.80),
             label=c("italic(β)['NFD'] == '6.1365'",
                     "'95%CI' == '[4.9817, 7.2916]'"),
             parse=T,size=3.5)+
    theme_regular()
)


##### dominance predictive curves for mean nearest differences #####
### get dominance predict data for mnnd
lincombs.data.domin.mnnd_single = data.frame(mnnd=seq(min(dat_dom_sp$mnnd),
                                                      max(dat_dom_sp$mnnd),
                                                      length=100))

lincombs.matrix.domin.mnnd_single=model.matrix(~mnnd,
                                               data=lincombs.data.domin.mnnd_single)
lincombs.matrix.domin.mnnd_single=as.data.frame(lincombs.matrix.domin.mnnd_single)
lincombs.domin.mnnd_single=inla.make.lincombs(lincombs.matrix.domin.mnnd_single)

inla.model_lincombs.domin.mnnd_single = pglmm(domin ~ mnnd+#(1|species) + 
                                                (1|f_p) + (1|field), data = dat_dom_sp,
                                              family = "binomial", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.domin.mnnd_single,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)

lincombs.posterior.domin.mnnd_single = inla.model_lincombs.domin.mnnd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnnd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnnd_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mnnd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnnd_single$lower=unlist(lapply(lincombs.posterior.domin.mnnd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnnd_single$upper=unlist(lapply(lincombs.posterior.domin.mnnd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnnd_single
save(lincombs.data.domin.mnnd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnnd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnnd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnnd_single.rdata")
(domin.mnnd_single.logistic=ggplot(data=lincombs.data.domin.mnnd_single,aes(x=mnnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnnd, y=domin),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Niche difference (ND)', y='Dominance probability')+
    annotate(geom="text",x=min(dat_dom_sp$mnnd)+
               (max(dat_dom_sp$mnnd)-min(dat_dom_sp$mnnd))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['ND'[nearest]] == '5.11'",
                     "'95%CI' == '[1.51, 8.71]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mnnd        5.11       1.51       8.71



### get dominance predict data for mnlgfd
lincombs.data.domin.mnlgfd_single = data.frame(mnlgfd=seq(min(dat_dom_sp$mnlgfd),
                                                          max(dat_dom_sp$mnlgfd),
                                                          length=100))

lincombs.matrix.domin.mnlgfd_single=model.matrix(~mnlgfd,
                                                 data=lincombs.data.domin.mnlgfd_single)
lincombs.matrix.domin.mnlgfd_single=as.data.frame(lincombs.matrix.domin.mnlgfd_single)
lincombs.domin.mnlgfd_single=inla.make.lincombs(lincombs.matrix.domin.mnlgfd_single)

inla.model_lincombs.domin.mnlgfd_single = pglmm(domin ~ mnlgfd+#(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_dom_sp,
                                                family = "binomial", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.domin.mnlgfd_single,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

lincombs.posterior.domin.mnlgfd_single = inla.model_lincombs.domin.mnlgfd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnlgfd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnlgfd_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mnlgfd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnlgfd_single$lower=unlist(lapply(lincombs.posterior.domin.mnlgfd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnlgfd_single$upper=unlist(lapply(lincombs.posterior.domin.mnlgfd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnlgfd_single
save(lincombs.data.domin.mnlgfd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnlgfd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnlgfd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnlgfd_single.rdata")
(domin.mnlgfd_single.logistic=ggplot(data=lincombs.data.domin.mnlgfd_single,aes(x=mnlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnlgfd, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Relative fitness difference (RFD)', y='')+
    annotate(geom="text",x=min(dat_dom_sp$mnlgfd)+
               (max(dat_dom_sp$mnlgfd)-min(dat_dom_sp$mnlgfd))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['RFD'[nearest]] == '0.41'",
                     "'95%CI' == '[0.07, 0.75]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mnlgfd       0.41       0.07       0.75



### get dominance predict data for mntd_all
lincombs.data.domin.mntd_single = data.frame(mntd_all=seq(min(dat_dom_sp$mntd_all),
                                                          max(dat_dom_sp$mntd_all),
                                                          length=100))

lincombs.matrix.domin.mntd_single=model.matrix(~mntd_all,
                                               data=lincombs.data.domin.mntd_single)
lincombs.matrix.domin.mntd_single=as.data.frame(lincombs.matrix.domin.mntd_single)
lincombs.domin.mntd_single=inla.make.lincombs(lincombs.matrix.domin.mntd_single)

inla.model_lincombs.domin.mntd_single = pglmm(domin ~ mntd_all+#(1|species) + 
                                                (1|f_p) + (1|field), data = dat_dom_sp,
                                              family = "binomial", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.domin.mntd_single,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)

lincombs.posterior.domin.mntd_single = inla.model_lincombs.domin.mntd_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mntd_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mntd_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mntd_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mntd_single$lower=unlist(lapply(lincombs.posterior.domin.mntd_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mntd_single$upper=unlist(lapply(lincombs.posterior.domin.mntd_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mntd_single
save(lincombs.data.domin.mntd_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mntd_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mntd_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mntd_single.rdata")
(domin.mntd_single.logistic=ggplot(data=lincombs.data.domin.mntd_single,aes(x=mntd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=2,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mntd_all, y=domin),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Phylogenetic difference (PD)', y='Dominance probability')+
    annotate(geom="text",x=min(dat_dom_sp$mntd_all)+
               (max(dat_dom_sp$mntd_all)-min(dat_dom_sp$mntd_all))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['NPD'] == '-0.0054'",
                     "'95%CI' == '[-0.0070, -0.0038]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mntd_all    -0.0054    -0.0070    -0.0038



### get dominance predict data for mnconti_func_d_all
lincombs.data.domin.mnconti_func_d_single = data.frame(mnconti_func_d_all = seq(min(dat_dom_sp$mnconti_func_d_all),
                                                                                max(dat_dom_sp$mnconti_func_d_all),length=100))

lincombs.matrix.domin.mnconti_func_d_single=model.matrix(~mnconti_func_d_all,
                                                         data=lincombs.data.domin.mnconti_func_d_single)
lincombs.matrix.domin.mnconti_func_d_single=as.data.frame(lincombs.matrix.domin.mnconti_func_d_single)
lincombs.domin.mnconti_func_d_single=inla.make.lincombs(lincombs.matrix.domin.mnconti_func_d_single)

inla.model_lincombs.domin.mnconti_func_d_single = pglmm(domin ~ mnconti_func_d_all+#(1|species) + 
                                                          (1|f_p) + (1|field), data = dat_dom_sp,
                                                        family = "binomial", cov_ranef = list(species = tree),
                                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                    config = TRUE),
                                                                             quantiles=c(0.025,0.5,0.975),
                                                                             lincomb=lincombs.domin.mnconti_func_d_single,
                                                                             control.predictor=list(compute=T)),
                                                        bayes = T)

lincombs.posterior.domin.mnconti_func_d_single = inla.model_lincombs.domin.mnconti_func_d_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnconti_func_d_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnconti_func_d_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mnconti_func_d_single,
                                                                        function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnconti_func_d_single$lower=unlist(lapply(lincombs.posterior.domin.mnconti_func_d_single,
                                                              function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnconti_func_d_single$upper=unlist(lapply(lincombs.posterior.domin.mnconti_func_d_single,
                                                              function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnconti_func_d_single
save(lincombs.data.domin.mnconti_func_d_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnconti_func_d_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnconti_func_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnconti_func_d_single.rdata")
(domin.mnconti_func_d_single.logistic=ggplot(data=lincombs.data.domin.mnconti_func_d_single,
                                             aes(x=mnconti_func_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=2)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #scale_x_continuous(breaks = c(0, 0.1, 0.2))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnconti_func_d_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Functional difference (FD)', y='  ')+
    annotate(geom="text",x=min(dat_dom_sp$mnconti_func_d_all)+
               (max(dat_dom_sp$mnconti_func_d_all)-min(dat_dom_sp$mnconti_func_d_all))*0.5,
             y=c(0.90,0.80),
             label=c("italic(β)['NFD'] == '0.6679'",
                     "'95%CI' == '[0.3004, 1.0355]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mnconti_func_d_all  0.6679     0.3004     1.0355




### get dominance predict data for mnfunc_d_all
lincombs.data.domin.mnfunc_d_single = data.frame(mnfunc_d_all = seq(min(dat_dom_sp$mnfunc_d_all),
                                                                    max(dat_dom_sp$mnfunc_d_all),length=100))

lincombs.matrix.domin.mnfunc_d_single=model.matrix(~mnfunc_d_all,
                                                   data=lincombs.data.domin.mnfunc_d_single)
lincombs.matrix.domin.mnfunc_d_single=as.data.frame(lincombs.matrix.domin.mnfunc_d_single)
lincombs.domin.mnfunc_d_single=inla.make.lincombs(lincombs.matrix.domin.mnfunc_d_single)

inla.model_lincombs.domin.mnfunc_d_single = pglmm(domin ~ mnfunc_d_all+#(1|species) + 
                                                    (1|f_p) + (1|field), data = dat_dom_sp,
                                                  family = "binomial", cov_ranef = list(species = tree),
                                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                              config = TRUE),
                                                                       quantiles=c(0.025,0.5,0.975),
                                                                       lincomb=lincombs.domin.mnfunc_d_single,
                                                                       control.predictor=list(compute=T)),
                                                  bayes = T)

lincombs.posterior.domin.mnfunc_d_single = inla.model_lincombs.domin.mnfunc_d_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnfunc_d_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnfunc_d_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mnfunc_d_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnfunc_d_single$lower=unlist(lapply(lincombs.posterior.domin.mnfunc_d_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnfunc_d_single$upper=unlist(lapply(lincombs.posterior.domin.mnfunc_d_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnfunc_d_single
save(lincombs.data.domin.mnfunc_d_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnfunc_d_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnfunc_d_single.rdata")
(domin.mnfunc_d_single.logistic=ggplot(data=lincombs.data.domin.mnfunc_d_single,
                                       aes(x=mnfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_x_continuous(breaks = c(0, 0.1, 0.2))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnfunc_d_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Functional difference (FD)', y='  ')+
    annotate(geom="text",x=c(0.08,0.08),y=c(0.90,0.80),
             label=c("italic(β)['NFD'] == '7.9398'",
                     "'95%CI' == '[5.3256, 10.5570]'"),
             parse=T,size=3.5)+
    theme_regular()
)



##### predictive curves for intrinsic growth rate #####
setwd("D:/R projects/BSS")
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


### get establishment predict data for demo_rate
lincombs.data.estab.demo_rate_single = data.frame(demo_rate=seq(min(dat_suc_sp$demo_rate),
                                                                max(dat_suc_sp$demo_rate),
                                                                length=100))

lincombs.matrix.estab.demo_rate_single=model.matrix(~demo_rate,
                                                    data=lincombs.data.estab.demo_rate_single)
lincombs.matrix.estab.demo_rate_single=as.data.frame(lincombs.matrix.estab.demo_rate_single)
lincombs.estab.demo_rate_single=inla.make.lincombs(lincombs.matrix.estab.demo_rate_single)

inla.model_lincombs.estab.demo_rate_single = pglmm(estab ~ demo_rate+#(1|species) + 
                                                     (1|f_p) + (1|field), data = dat_suc_sp,
                                                   family = "binomial", cov_ranef = list(species = tree),
                                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                               config = TRUE),
                                                                        quantiles=c(0.025,0.5,0.975),
                                                                        lincomb=lincombs.estab.demo_rate_single,
                                                                        control.predictor=list(compute=T)),
                                                   bayes = T)

lincombs.posterior.estab.demo_rate_single = inla.model_lincombs.estab.demo_rate_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.demo_rate_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.demo_rate_single$predicted.value=unlist(lapply(lincombs.posterior.estab.demo_rate_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.demo_rate_single$lower=unlist(lapply(lincombs.posterior.estab.demo_rate_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.demo_rate_single$upper=unlist(lapply(lincombs.posterior.estab.demo_rate_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.demo_rate_single
save(lincombs.data.estab.demo_rate_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.demo_rate_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ demo_rate
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.demo_rate_single.rdata")
(estab.demo_rate_single.logistic=ggplot(data=lincombs.data.estab.demo_rate_single,aes(x=demo_rate, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='black',alpha=0.2)+
    geom_line(color='black',size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=demo_rate, y=estab),
               color = 'grey',
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",x=c(2.3, 2.3),y=c(0.90,0.80),
             label=c("italic(β)['Intrinsic growth rate'] == '1.02'",
                     "'95%CI' == '[0.89, 1.15]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#demo_rate    1.02       0.89       1.15


### get dominance predict data for demo_rate
lincombs.data.domin.demo_rate_single = data.frame(demo_rate=seq(min(dat_dom_sp$demo_rate),
                                                                max(dat_dom_sp$demo_rate),
                                                                length=100))

lincombs.matrix.domin.demo_rate_single=model.matrix(~demo_rate,
                                                    data=lincombs.data.domin.demo_rate_single)
lincombs.matrix.domin.demo_rate_single=as.data.frame(lincombs.matrix.domin.demo_rate_single)
lincombs.domin.demo_rate_single=inla.make.lincombs(lincombs.matrix.domin.demo_rate_single)

inla.model_lincombs.domin.demo_rate_single = pglmm(domin ~ demo_rate+#(1|species) + 
                                                     (1|f_p) + (1|field), data = dat_dom_sp,
                                                   family = "binomial", cov_ranef = list(species = tree),
                                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                               config = TRUE),
                                                                        quantiles=c(0.025,0.5,0.975),
                                                                        lincomb=lincombs.domin.demo_rate_single,
                                                                        control.predictor=list(compute=T)),
                                                   bayes = T)

lincombs.posterior.domin.demo_rate_single = inla.model_lincombs.domin.demo_rate_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.demo_rate_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

# inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.demo_rate_single$predicted.value=unlist(lapply(lincombs.posterior.domin.demo_rate_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.demo_rate_single$lower=unlist(lapply(lincombs.posterior.domin.demo_rate_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.demo_rate_single$upper=unlist(lapply(lincombs.posterior.domin.demo_rate_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.demo_rate_single
save(lincombs.data.domin.demo_rate_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.demo_rate_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominlishment ~ demo_rate
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.demo_rate_single.rdata")
(domin.demo_rate_single.logistic=ggplot(data=lincombs.data.domin.demo_rate_single,aes(x=demo_rate, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='black',alpha=0.2)+
    geom_line(color='black',size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=demo_rate, y=domin),
               color = 'grey',
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Intrinsic growth rate', y='Dominance probability')+
    annotate(geom="text",x=c(2.8, 2.8),y=c(0.90,0.80),
             label=c("italic(β)['Intrinsic growth rate'] == '-0.55'",
                     "'95%CI' == '[-0.80, -0.30]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#demo_rate   -0.55      -0.80      -0.30





#### Fast start for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for all species ####
setwd("D:/R projects/BSS")
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

car::vif(
  glmer(estab ~ mnnd + mnlgfd + mntd_all + mnconti_func_d_all + (1|f_p)+ (1|field),
        family=binomial,data=dat_suc_sps)
)

### all functional trait distance
car::vif(
  glmer(estab ~ mnd + mlgfd + mpd_all + mfunc_d_all + (1|f_p)+ (1|field),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnd.a + mlgfd.a + mpd.a_all + mfunc_d.a_all + (1|f_p)+ (1|field),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnnd + mnlgfd + mntd_all + mnfunc_d_all + (1|species) + (1|f_p)+ (1|field),
        family=binomial,data=dat_suc_sps)
)
###no co-linearity problem, all VIF < 3

### only continuous trait distance 
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

estab_model_mnd_conti_func_d_all  = pglmm(estab~mnnd+mnlgfd+mntd_all+mnconti_func_d_all
                                          #+(1|species) 
                                          + (1|f_p) + (1|field), data = dat_suc_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,
                                                                                      waic=T,
                                                                                      cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)

summary(estab_model_md_conti_func_d_all)
summary(estab_model_md.a_conti_func_d_all)
summary(estab_model_mnd_conti_func_d_all)

### all functional trait distance
estab_model_md_func_d_all = pglmm(estab~mnd+mlgfd+mpd_all+mfunc_d_all+#(1|species) + 
                                    (1|f_p) + (1|field), data = dat_suc_sps,
                                  family = "binomial", cov_ranef = list(species = tree),
                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                              config = TRUE),
                                                       quantiles=c(0.025,0.5,0.975)),
                                  bayes = T)

estab_model_md.a_func_d_all = pglmm(estab~mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all+#(1|species) + 
                                      (1|f_p) + (1|field), data = dat_suc_sps,
                                    family = "binomial", cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975)),
                                    bayes = T)

estab_model_mnd_func_d_all = pglmm(estab~mnnd+mnlgfd+mntd_all+mnfunc_d_all+#(1|species) + 
                                     (1|f_p) + (1|field), data = dat_suc_sps,
                                   family = "binomial", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T)

summary(estab_model_md_func_d_all)
summary(estab_model_md.a_func_d_all)
summary(estab_model_mnd_func_d_all)

### split fd to demo_ratio and comp_ratio
estab_model_md_ratio_func_d_all = pglmm(estab~mnd+mlgfd
                                        +mpd_all+mfunc_d_all+
                                          demo_rate+ #+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sps,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = T),
                                                             quantiles=c(0.025,0.5,0.975)),
                                        bayes = T)

estab_model_md.a_ratio_func_d_all = pglmm(estab~mnd.a+mlgfd.a
                                          +mpd.a_all+mfunc_d.a_all+
                                            demo_rate+ #+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_suc_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = T),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)

estab_model_mnd_ratio_func_d_all = pglmm(estab~mnnd+mnlgfd
                                         +mntd_all+mnfunc_d_all+
                                           demo_rate+ #+(1|species) + 
                                           (1|f_p) + (1|field), data = dat_suc_sps,
                                         family = "binomial", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = T),
                                                              quantiles=c(0.025,0.5,0.975)),
                                         bayes = T)
summary(estab_model_md_ratio_func_d_all)
summary(estab_model_md.a_ratio_func_d_all)
summary(estab_model_mnd_ratio_func_d_all)

# plot 
estab_data.inla.md_func_d.all_intercept1 = estab_model_md_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.md_func_d.all_intercept = estab_data.inla.md_func_d.all_intercept1%>%
  mutate(rowname=c("MND","MRFD","MPD_all","MFD_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

estab_data.inla.md.a_func_d.all_intercept1 = estab_model_md.a_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.md.a_func_d.all_intercept = estab_data.inla.md.a_func_d.all_intercept1%>%
  mutate(rowname=c("MND.ab","MRFD.ab","MPD.ab_all","MFD.ab_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

estab_data.inla.mnd_func_d.all_intercept1 = estab_model_mnd_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.mnd_func_d.all_intercept = estab_data.inla.mnd_func_d.all_intercept1%>%
  mutate(rowname=c("MNND","MNRFD","MNTD_all","MNFD_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

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

estab_data.inla.mnd_conti_func_d.all_intercept1 = estab_model_mnd_conti_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.mnd_conti_func_d.all_intercept = estab_data.inla.mnd_conti_func_d.all_intercept1%>%
  mutate(rowname=c("MNND","MNRFD","MNTD_all","MNFD_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))


### Effect size plot for mean differences
# point + effect size
require(ggplot2)
require(ggpubr)
(estab_nd_rfd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(estab_data.inla.md_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MRFD', 'MND'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=(estab_data.inla.md_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MRFD', 'MND'))),
               mapping = aes(x=mean,y=rowname),
               size = 4.4)+
    geom_errorbar(data=estab_data.inla.md_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MRFD', 'MND')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD', 'MND'),
                     label =c('RFD', 'ND'))+
    scale_x_continuous(limits=c(-0.15,0.5))+
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
    annotate(geom="text", x=0.32, y=2.45, family = 'Arial',
             label='Establishment', size = 6)+
    labs(x = '  ', y = NULL)+
    guides(color="none")+
    theme_regular())

(estab_pd_fd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(estab_data.inla.md_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MFD_all', 'MPD_all'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=(estab_data.inla.md_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MFD_all', 'MPD_all'))),
               mapping = aes(x=mean,y=rowname),
               size = 4.4)+
    geom_errorbar(data=estab_data.inla.md_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MFD_all', 'MPD_all')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD_all','MPD_all'),
                     label=c('MFD','MPD'))+
    scale_x_continuous(limits=c(-1,0.15), breaks = seq(-1, 0, 0.5))+
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
    annotate(geom="text", x=-0.7, y=2.45, family = 'Arial',
             label='Establishment', size = 6)+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

# R2%
(estab_md.all.varied.intercept.R2.plot = 
    ggplot(data=estab_data.inla.md_func_d.all_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MFD_all','MPD_all', 'MRFD', 'MND'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
established_all_md.all = ggarrange(estab_md.all.varied.intercept.plot,
                                   estab_md.all.varied.intercept.R2.plot,
                                   widths=c(2,1.2))


### Effect size plot for abundance weighted mean differences
# point + effect size
(estab_nd_rfd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(estab_data.inla.md.a_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MRFD.ab', 'MND.ab'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=(estab_data.inla.md.a_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MRFD.ab', 'MND.ab'))),
               mapping = aes(x=mean,y=rowname),
               size = 4.4)+
    geom_errorbar(data=estab_data.inla.md.a_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MRFD.ab', 'MND.ab')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD.ab', 'MND.ab'),
                     label = c(expression(RFD['ab']), expression(ND['ab'])))+
    scale_x_continuous(limits=c(-0.15,0.5))+
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
    annotate(geom="text", x=0.295, y=2.45, family = 'Arial',
             label='Establishment',
             size = 6)+
    labs(x = '  ', y = '  ')+
    guides(color="none")+
    theme_regular())

(estab_pd_fd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(estab_data.inla.md.a_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=(estab_data.inla.md.a_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all'))),
               mapping = aes(x=mean,y=rowname),
               size = 4.4)+
    geom_errorbar(data=estab_data.inla.md.a_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all'),
                     label=c(expression(MFD['ab']), expression(MPD['ab'])))+
    scale_x_continuous(limits=c(-1,0.15), breaks = seq(-1, 0, 0.5))+
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
    annotate(geom="text", x=-0.68, y=2.45, family = 'Arial',
             label='Establishment', size = 6)+
    labs(x = '  ', y = ' ')+
    guides(color="none")+
    theme_regular())

(estab_md.a_ratio.all.varied.intercept.plot =
    ggplot(data=estab_data.inla.md.a_ratio_func_d.all_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 3)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all', 'MRFD.ab', 'MND.ab', 'demo_rate'),
                     label = c(expression(MFD['ab']), expression(MPD['ab']),
                               expression(RFD['ab']), expression(ND['ab']),
                               'Intrinsic growth rate'))+
    scale_x_continuous(limits=c(-1.5,0.5))+
    scale_color_manual(values = c('purple',"#80defb", "#5a57fe","#fbd76c", "#f55756"),
                       name = ' ')+
    annotate(geom="text", x=1, y=5.5, family = 'Arial',
             label='Establishment')+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())

(estab_demo_rate.a.all.varied.intercept.plot =
    ggplot(data=(data=estab_data.inla.md.a_ratio_func_d.all_intercept %>% 
                   filter(rowname %in% c())),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 3)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('demo_rate'),
                     label = c('Intrinsic growth rate'))+
    scale_x_continuous(limits=c(-1,1.1))+
    
    annotate(geom="text", x=0.6, y=1.5, family = 'Arial',
             label='Establishment')+
    labs(x = 'Standardized effects', y = '  ')+
    guides(color="none")+
    theme_regular())

# R2%
(estab_md.a.all.varied.intercept.R2.plot = 
    ggplot(data=estab_data.inla.md.a_func_d.all_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all', 'MRFD.ab', 'MND.ab'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
established_all_md.a.all = ggarrange(estab_md.a.all.varied.intercept.plot,
                                     estab_md.a.all.varied.intercept.R2.plot,
                                     widths=c(2,1.2))


### Effect size plot for mean nearest differences
# point + effect size
(estab_nnd_nrfd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(estab_data.inla.mnd_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MNRFD', 'MNND'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=(estab_data.inla.mnd_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MNRFD', 'MNND'))),
               mapping = aes(x=mean,y=rowname),
               size = 4.4)+
    geom_errorbar(data=estab_data.inla.mnd_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MNRFD', 'MNND')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNRFD', 'MNND'),
                     label = c(expression(RFD['nearest']), expression(ND['nearest'])))+
    scale_x_continuous(limits=c(-0.15,0.5))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MNND" = "MNND",
                                  "MNRFD" = "MNRFD"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[1],
                                 colors_4d[2]),
                      labels = c("MNND" = "MNND",
                                 "MNRFD" = "MNRFD"),
                      name = ' ')+
    annotate(geom="text", x=0.31, y=2.45, family = 'Arial',
             label='Establishment', size = 6)+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

(estab_npd_nfd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(estab_data.inla.mnd_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MNFD_all', 'MNTD_all'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=(estab_data.inla.mnd_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MNFD_all', 'MNTD_all'))),
               mapping = aes(x=mean,y=rowname),
               size = 4.4)+
    geom_errorbar(data=estab_data.inla.mnd_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MNFD_all', 'MNTD_all')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNFD_all','MNTD_all'),
                     label = c('NFD','NPD'))+
    scale_x_continuous(limits=c(-1,0.5), breaks = seq(-1, 0.5, 0.5))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MNTD_all" = "MNTD_all",
                                  "MNFD_all" = "MNFD_all"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[4],
                                 colors_4d[3]),
                      labels = c("MNTD_all" = "MNTD_all",
                                 "MNFD_all" = "MNFD_all"),
                      name = ' ')+
    annotate(geom="text", x=-0.61, y=2.45, family = 'Arial',
             label='Establishment', size = 6)+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

# R2%
(estab_mnd.all.varied.intercept.R2.plot = 
    ggplot(data=estab_data.inla.mnd_func_d.all_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MNFD_all','MNTD_all', 'MNRFD', 'MNND'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
established_all_mnd.all = ggarrange(estab_mnd.all.varied.intercept.plot,
                                    estab_mnd.all.varied.intercept.R2.plot,
                                    widths=c(2,1.2))

###### establishment predictive curves for md ####
### get establishment predict data for mnd
lincombs.data.estab.mnd = data.frame(mnd=seq(min(dat_suc_sp$mnd),max(dat_suc_sp$mnd),length=100),
                                     mlgfd=mean(dat_suc_sp$mlgfd),
                                     mpd_all = mean(dat_suc_sp$mpd_all),
                                     mfunc_d_all = mean(dat_suc_sp$mfunc_d_all))

lincombs.matrix.estab.mnd=model.matrix(~mnd+mlgfd+mpd_all+mfunc_d_all,
                                       data=lincombs.data.estab.mnd)
lincombs.matrix.estab.mnd=as.data.frame(lincombs.matrix.estab.mnd)
lincombs.estab.mnd=inla.make.lincombs(lincombs.matrix.estab.mnd)

inla.model_lincombs.estab.mnd = pglmm(estab ~ mnd+mlgfd+mpd_all+mfunc_d_all+(1|species) + 
                                        (1|f_p) + (1|field), data = dat_suc_sp,
                                      family = "binomial", cov_ranef = list(species = tree),
                                      bayes_options = list(control.compute = list(dic=T,
                                                                                  waic=T,
                                                                                  cpo=T,
                                                                                  config = TRUE),
                                                           quantiles=c(0.025,0.5,0.975),
                                                           lincomb=lincombs.estab.mnd,
                                                           control.predictor=list(compute=T)),
                                      bayes = T)

lincombs.posterior.estab.mnd = inla.model_lincombs.estab.mnd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.mnd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnd$predicted.value=unlist(lapply(lincombs.posterior.estab.mnd,
                                                      function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnd$lower=unlist(lapply(lincombs.posterior.estab.mnd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnd$upper=unlist(lapply(lincombs.posterior.estab.mnd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

save(lincombs.data.estab.mnd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment~mnd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnd.rdata")

(estab.mnd.partial.logistic=ggplot(data=lincombs.data.estab.mnd,
                                   aes(x=mnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x='Introduced-native niche difference',
         y="Establishment probability")+
    annotate(geom="text",x=c(0.33,0.33),y=c(0.80,0.70),
             label=c("italic(β)['I-N MND'] == -2.66",
                     "'95%CI' == '[-3.52, -1.81]'"),parse=T,size=3.5)
)


### get establishment predict data for mlgfd
lincombs.data.estab.mlgfd = data.frame(mlgfd=seq(min(dat_suc_sp$mlgfd),max(dat_suc_sp$mlgfd),length=100),
                                       mnd=mean(dat_suc_sp$mnd),
                                       mpd_all = mean(dat_suc_sp$mpd_all),
                                       mfunc_d_all = mean(dat_suc_sp$mfunc_d_all))

lincombs.matrix.estab.mlgfd=model.matrix(~mnd+mlgfd+mpd_all+mfunc_d_all,
                                         data=lincombs.data.estab.mlgfd)
lincombs.matrix.estab.mlgfd=as.data.frame(lincombs.matrix.estab.mlgfd)
lincombs.estab.mlgfd=inla.make.lincombs(lincombs.matrix.estab.mlgfd)

inla.model_lincombs.estab.mlgfd = pglmm(estab ~ mnd+mlgfd+mpd_all+mfunc_d_all+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.estab.mlgfd,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)

lincombs.posterior.estab.mlgfd = inla.model_lincombs.estab.mlgfd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mlgfd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mlgfd$predicted.value=unlist(lapply(lincombs.posterior.estab.mlgfd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mlgfd$lower=unlist(lapply(lincombs.posterior.estab.mlgfd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mlgfd$upper=unlist(lapply(lincombs.posterior.estab.mlgfd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mlgfd
save(lincombs.data.estab.mlgfd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mlgfd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd.rdata")
(estab.mlgfd.partial.logistic=ggplot(data=lincombs.data.estab.mlgfd,aes(x=mlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mlgfd, y=estab),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native fitness difference', y='  ')+
    annotate(geom="text",x=c(0,0),y=c(0.80,0.70),
             label=c("italic(β)['I-N MRFD'] == 2.61",
                     "'95%CI' == '[1.74, 3.49]'"),
             parse=T,size=3.5)
)

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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mpd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mpd_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mpd_all.rdata")
(estab.mpd_all.partial.logistic=ggplot(data=lincombs.data.estab.mpd_all,aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=estab),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native phylogenetic difference', y='Establishment probability')+
    annotate(geom="text",x=c(200,200),y=c(0.80,0.70),
             label=c("italic(β)['I-N MPD'] == '-0.02'",
                     "'95%CI' == '[-0.03, -0.01]'"),
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mfunc_d_all.rdata")
(estab.mfunc_d_all.partial.logistic=ggplot(data=lincombs.data.estab.mfunc_d_all,
                                           aes(x=mfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1,
              linetype= 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d_all, y=estab),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native functional difference', y='  ')+
    annotate(geom="text",x=c(0.3,0.3),y=c(0.80,0.70),
             label=c("italic(β)['I-N MFD'] == 5.53",
                     "'95%CI' == '[-0.07, 11.08]'"),
             parse=T,size=3.5)
)


###### establishment predictive curves for mean nearest difference ####
### get establishment predict data for mnnd
lincombs.data.estab.mnnd = data.frame(mnnd=seq(min(dat_suc_sp$mnnd),max(dat_suc_sp$mnnd),length=100),
                                      mnlgfd=mean(dat_suc_sp$mnlgfd),
                                      mntd_all = mean(dat_suc_sp$mntd_all),
                                      mnfunc_d_all = mean(dat_suc_sp$mnfunc_d_all))

lincombs.matrix.estab.mnnd=model.matrix(~mnnd+mnlgfd+mntd_all+mnfunc_d_all,
                                        data=lincombs.data.estab.mnnd)
lincombs.matrix.estab.mnnd=as.data.frame(lincombs.matrix.estab.mnnd)
lincombs.estab.mnnd=inla.make.lincombs(lincombs.matrix.estab.mnnd)

inla.model_lincombs.estab.mnnd = pglmm(estab ~ mnnd+mnlgfd+mntd_all+mnfunc_d_all+(1|species) + 
                                         (1|f_p) + (1|field), data = dat_suc_sp,
                                       family = "binomial", cov_ranef = list(species = tree),
                                       bayes_options = list(control.compute = list(dic=T,
                                                                                   waic=T,
                                                                                   cpo=T,
                                                                                   config = TRUE),
                                                            quantiles=c(0.025,0.5,0.975),
                                                            lincomb=lincombs.estab.mnnd,
                                                            control.predictor=list(compute=T)),
                                       bayes = T)

lincombs.posterior.estab.mnnd = inla.model_lincombs.estab.mnnd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnnd$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.mnnd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnnd$predicted.value=unlist(lapply(lincombs.posterior.estab.mnnd,
                                                       function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnnd$lower=unlist(lapply(lincombs.posterior.estab.mnnd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnnd$upper=unlist(lapply(lincombs.posterior.estab.mnnd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

save(lincombs.data.estab.mnnd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment~mnnd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnnd.rdata")

(estab.mnnd.partial.logistic=ggplot(data=lincombs.data.estab.mnnd,
                                    aes(x=mnnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnnd, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x='Introduced-native niche difference',
         y="Establishment probability")+
    annotate(geom="text",x=c(0.33,0.33),y=c(0.80,0.70),
             label=c("italic(β)['I-N MNND'] == -0.9087",
                     "'95%CI' == '[-1.2532, -0.5641]'"),parse=T,size=3.5)
)


### get establishment predict data for mnlgfd
lincombs.data.estab.mnlgfd = data.frame(mnlgfd=seq(min(dat_suc_sp$mnlgfd),max(dat_suc_sp$mnlgfd),length=100),
                                        mnnd=mean(dat_suc_sp$mnnd),
                                        mntd_all = mean(dat_suc_sp$mntd_all),
                                        mnfunc_d_all = mean(dat_suc_sp$mnfunc_d_all))

lincombs.matrix.estab.mnlgfd=model.matrix(~mnnd+mnlgfd+mntd_all+mnfunc_d_all,
                                          data=lincombs.data.estab.mnlgfd)
lincombs.matrix.estab.mnlgfd=as.data.frame(lincombs.matrix.estab.mnlgfd)
lincombs.estab.mnlgfd=inla.make.lincombs(lincombs.matrix.estab.mnlgfd)

inla.model_lincombs.estab.mnlgfd = pglmm(estab ~ mnnd+mnlgfd+mntd_all+mnfunc_d_all+(1|species) + 
                                           (1|f_p) + (1|field), data = dat_suc_sp,
                                         family = "binomial", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = TRUE),
                                                              quantiles=c(0.025,0.5,0.975),
                                                              lincomb=lincombs.estab.mnlgfd,
                                                              control.predictor=list(compute=T)),
                                         bayes = T)

lincombs.posterior.estab.mnlgfd = inla.model_lincombs.estab.mnlgfd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnlgfd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnlgfd$predicted.value=unlist(lapply(lincombs.posterior.estab.mnlgfd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnlgfd$lower=unlist(lapply(lincombs.posterior.estab.mnlgfd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnlgfd$upper=unlist(lapply(lincombs.posterior.estab.mnlgfd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnlgfd
save(lincombs.data.estab.mnlgfd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnlgfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnlgfd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnlgfd.rdata")
(estab.mnlgfd.partial.logistic=ggplot(data=lincombs.data.estab.mnlgfd,aes(x=mnlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnlgfd, y=estab),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native fitness difference', y='  ')+
    annotate(geom="text",x=c(0,0),y=c(0.80,0.70),
             label=c("italic(β)['I-N MNRFD'] == 1.82",
                     "'95%CI' == '[1.25, 2.39]'"),
             parse=T,size=3.5)
)

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
inla.model_lincombs.estab.mntd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mntd_all$predicted.value=unlist(lapply(lincombs.posterior.estab.mntd_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mntd_all$lower=unlist(lapply(lincombs.posterior.estab.mntd_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mntd_all$upper=unlist(lapply(lincombs.posterior.estab.mntd_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mntd_all
save(lincombs.data.estab.mntd_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mntd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mntd_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mntd_all.rdata")
(estab.mntd_all.partial.logistic=ggplot(data=lincombs.data.estab.mntd_all,aes(x=mntd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mntd_all, y=estab),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native phylogenetic difference', y='Establishment probability')+
    annotate(geom="text",x=c(200,200),y=c(0.80,0.70),
             label=c("italic(β)['I-N MNTD'] == '-0.0026'",
                     "'95%CI' == '[-0.0050, -0.0001]'"),
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnfunc_d_all.rdata")
(estab.mnfunc_d_all.partial.logistic=ggplot(data=lincombs.data.estab.mnfunc_d_all,
                                            aes(x=mnfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1,
              linetype= 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnfunc_d_all, y=estab),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native functional difference', y='  ')+
    annotate(geom="text",x=c(0.15,0.15),y=c(0.80,0.70),
             label=c("italic(β)['I-N MNFD'] == -1.54",
                     "'95%CI' == '[-4.87, 1.78]'"),
             parse=T,size=3.5)
)



###### establishment predictive curves for md.ab ####
### get establishment predict data for mnd.ab
lincombs.data.estab.mnd.a = data.frame(mnd.a=seq(min(dat_suc_sp$mnd.a),max(dat_suc_sp$mnd.a),length=100),
                                       mlgfd.a=mean(dat_suc_sp$mlgfd.a),
                                       mpd.a_all = mean(dat_suc_sp$mpd.a_all),
                                       mfunc_d.a_all = mean(dat_suc_sp$mfunc_d.a_all))

lincombs.matrix.estab.mnd.a=model.matrix(~mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all,
                                         data=lincombs.data.estab.mnd.a)
lincombs.matrix.estab.mnd.a=as.data.frame(lincombs.matrix.estab.mnd.a)
lincombs.estab.mnd.a=inla.make.lincombs(lincombs.matrix.estab.mnd.a)

inla.model_lincombs.estab.mnd.a = pglmm(estab ~ mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,
                                                                                    waic=T,
                                                                                    cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.estab.mnd.a,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)

lincombs.posterior.estab.mnd.a = inla.model_lincombs.estab.mnd.a$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnd.a$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.mnd.a$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnd.a$predicted.value=unlist(lapply(lincombs.posterior.estab.mnd.a,
                                                        function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnd.a$lower=unlist(lapply(lincombs.posterior.estab.mnd.a,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnd.a$upper=unlist(lapply(lincombs.posterior.estab.mnd.a,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

save(lincombs.data.estab.mnd.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnd.a.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment~mnd.ab
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mnd.a.rdata")
library(scales)
(estab.mnd.a.partial.logistic=ggplot(data=lincombs.data.estab.mnd.a,
                                     aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd.a, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x=' ',
         y="Establishment probability")+
    annotate(geom="text",x=c(0.33,0.33),y=c(0.80,0.70),
             label=c("italic(β)['MND'[ab]] == -2.12",
                     "'95%CI' == '[-2.77, -1.46]'"),parse=T,size=3.5)+    
    theme_regular()
)


### get establishment predict data for mlgfd.ab
lincombs.data.estab.mlgfd.a = data.frame(mlgfd.a=seq(min(dat_suc_sp$mlgfd.a),max(dat_suc_sp$mlgfd.a),length=100),
                                         mnd.a=mean(dat_suc_sp$mnd.a),
                                         mpd.a_all = mean(dat_suc_sp$mpd.a_all),
                                         mfunc_d.a_all = mean(dat_suc_sp$mfunc_d.a_all))

lincombs.matrix.estab.mlgfd.a=model.matrix(~mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all,
                                           data=lincombs.data.estab.mlgfd.a)
lincombs.matrix.estab.mlgfd.a=as.data.frame(lincombs.matrix.estab.mlgfd.a)
lincombs.estab.mlgfd.a=inla.make.lincombs(lincombs.matrix.estab.mlgfd.a)

inla.model_lincombs.estab.mlgfd.a = pglmm(estab ~ mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_suc_sp,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.estab.mlgfd.a,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)

lincombs.posterior.estab.mlgfd.a = inla.model_lincombs.estab.mlgfd.a$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mlgfd.a$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mlgfd.a$predicted.value=unlist(lapply(lincombs.posterior.estab.mlgfd.a,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mlgfd.a$lower=unlist(lapply(lincombs.posterior.estab.mlgfd.a,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mlgfd.a$upper=unlist(lapply(lincombs.posterior.estab.mlgfd.a,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mlgfd.a
save(lincombs.data.estab.mlgfd.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd.a.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mlgfd.ab
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mlgfd.a.rdata")
(estab.mlgfd.a.partial.logistic=ggplot(data=lincombs.data.estab.mlgfd.a,aes(x=mlgfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mlgfd.a, y=estab),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=c(0,0),y=c(0.80,0.70),
             label=c("italic(β)['MRFD'[ab]] == 2.71", "'95%CI' == '[1.99, 3.42]'"),
             parse=T,size=3.5)+
    theme_regular()
)

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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mpd.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mpd.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mpd.a_all.rdata")
(estab.mpd.a_all.partial.logistic=ggplot(data=lincombs.data.estab.mpd.a_all,aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 3)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd.a_all, y=estab),
               color = colors_4d[3],
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mfunc_d.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfunc_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.estab.mfunc_d.a_all.rdata")
(estab.mfunc_d.a_all.partial.logistic=ggplot(data=lincombs.data.estab.mfunc_d.a_all,
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
    annotate(geom="text",x=c(0.3,0.3),y=c(0.80,0.70),
             label=c("italic(β)['MFD'[ab]] == -4.35", "'95%CI' == '[-6.40, -2.30]'"),
             parse=T,size=3.5)+
    theme_regular()
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

car::vif(
  glmer(domin ~ mnnd + mnlgfd + mntd_all + mnconti_func_d_all + 
          (1|f_p)+ (1|field),
        family=binomial,data=dat_dom_sps)
)
### no colinearity problem

### only continuous trait distance 
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
summary(domin_model_mnd_conti_func_d_all)

### all functional trait distance
domin_model_md_func_d_all = pglmm(domin~mnd+mlgfd
                                  +mpd_all+mfunc_d_all
                                  #+(1|species) 
                                  + (1|f_p) 
                                  + (1|field)
                                  , data = dat_dom_sps,
                                  family = "binomial", cov_ranef = list(species = tree),
                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                              config = TRUE),
                                                       quantiles=c(0.025,0.5,0.975)),
                                  bayes = T)

domin_model_md.a_func_d_all = pglmm(domin~mnd.a+mlgfd.a
                                    +mpd.a_all+mfunc_d.a_all
                                    #+(1|species)
                                    + (1|f_p) 
                                    + (1|field)
                                    , data = dat_dom_sps,
                                    family = "binomial", cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975)),
                                    bayes = T)

domin_model_mnd_func_d_all = pglmm(domin~mnnd+mnlgfd
                                   +mntd_all+mnfunc_d_all
                                   #+(1|species) 
                                   + (1|f_p) 
                                   + (1|field)
                                   , data = dat_dom_sps,
                                   family = "binomial", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T)
summary(domin_model_md_func_d_all)
summary(domin_model_md.a_func_d_all)
summary(domin_model_mnd_func_d_all)


### split fd to demo_ratio and comp_ratio
domin_model_md_ratio_func_d_all = pglmm(domin~mnd+mlgfd
                                        +mpd_all+mfunc_d_all+
                                          demo_rate+ #+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_dom_sps,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = T),
                                                             quantiles=c(0.025,0.5,0.975)),
                                        bayes = T)

domin_model_md.a_ratio_func_d_all = pglmm(domin~mnd.a+mlgfd.a
                                          +mpd.a_all+mfunc_d.a_all+
                                            demo_rate+ #(1|species) + 
                                            (1|f_p) + (1|field), data = dat_dom_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = T),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)

domin_model_mnd_ratio_func_d_all = pglmm(domin~mnnd+mnlgfd
                                         +mntd_all+mnfunc_d_all+
                                           demo_rate+ #+(1|species) + 
                                           (1|f_p) + (1|field), data = dat_dom_sps,
                                         family = "binomial", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = T),
                                                              quantiles=c(0.025,0.5,0.975)),
                                         bayes = T)
summary(domin_model_md_ratio_func_d_all)
summary(domin_model_md.a_ratio_func_d_all)
summary(domin_model_mnd_ratio_func_d_all)

# plot 
domin_data.inla.md_func_d.all_intercept1 = domin_model_md_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.md_func_d.all_intercept = domin_data.inla.md_func_d.all_intercept1%>%
  mutate(rowname=c("MND","MRFD","MPD_all","MFD_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

domin_data.inla.md.a_func_d.all_intercept1 = domin_model_md.a_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.md.a_func_d.all_intercept = domin_data.inla.md.a_func_d.all_intercept1%>%
  mutate(rowname=c("MND.ab","MRFD.ab","MPD.ab_all","MFD.ab_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

domin_data.inla.mnd_func_d.all_intercept1 = domin_model_mnd_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.mnd_func_d.all_intercept = domin_data.inla.mnd_func_d.all_intercept1%>%
  mutate(rowname=c("MNND","MNRFD","MNTD_all","MNFD_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

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

domin_data.inla.mnd_conti_func_d.all_intercept1 = domin_model_mnd_conti_func_d_all$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.mnd_conti_func_d.all_intercept = domin_data.inla.mnd_conti_func_d.all_intercept1%>%
  mutate(rowname=c("MNND","MNRFD","MNTD_all","MNFD_all"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))%>% 
  mutate(significant = ifelse(lower / abs(lower) == upper / abs(upper), 'yes', 'no'))

### Effect size plot for mean differences
# point + effect size
(domin_nd_rfd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(domin_data.inla.md_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MRFD', 'MND'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
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
    geom_point(data=(domin_data.inla.md_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MRFD', 'MND'))),
               mapping = aes(x=mean,y=rowname, color = significant),
               size = 4.4)+
    geom_errorbar(data=domin_data.inla.md_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MRFD', 'MND')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper, color = significant),
                  width = 0.12, linewidth = 1.2)+
    scale_color_manual(values = c('black',
                                  'grey'),
                       name = ' ')+
    ggnewscale::new_scale_color()+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD', 'MND'),
                     label =c('RFD', 'ND'))+
    scale_x_continuous(limits=c(-0.15,0.5))+
    
    
    annotate(geom="text", x=0.35, y=2.45, family = 'Arial',
             label='Dominance', size = 6)+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

(domin_pd_fd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(domin_data.inla.md_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MFD_all', 'MPD_all'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
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
    geom_point(data=(domin_data.inla.md_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MFD_all', 'MPD_all'))),
               mapping = aes(x=mean,y=rowname,color = significant),
               size = 4.4)+
    geom_errorbar(data=domin_data.inla.md_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MFD_all', 'MPD_all')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper,color = significant),
                  width = 0.12, linewidth = 1.2)+
    scale_color_manual(values = c('grey',
                                  'black'),
                       name = ' ')+
    ggnewscale::new_scale_color()+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD_all','MPD_all'),
                     label =c('MFD','MPD') )+
    scale_x_continuous(limits=c(-1,0.15), breaks = seq(-1, 0, 0.5))+
    annotate(geom="text", x=-0.76, y=2.45, family = 'Arial',
             label='Dominance', size = 6)+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

# R2%
(domin_md.all.varied.intercept.R2.plot = 
    ggplot(data=domin_data.inla.md_func_d.all_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MFD_all','MPD_all', 'MRFD', 'MND'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
dominance_all_md.all = ggarrange(domin_md.all.varied.intercept.plot,
                                 domin_md.all.varied.intercept.R2.plot,
                                 widths=c(2,1.2))


### Effect size plot for abundance weighted mean differences
# point + effect size
(domin_nd_rfd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(domin_data.inla.md.a_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MRFD.ab', 'MND.ab'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=(domin_data.inla.md.a_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MRFD.ab', 'MND.ab'))),
               mapping = aes(x=mean,y=rowname),
               size = 4.4)+
    geom_errorbar(data=domin_data.inla.md.a_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MRFD.ab', 'MND.ab')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD.ab', 'MND.ab'),
                     label = c(expression(RFD['ab']), expression(ND['ab'])))+
    scale_x_continuous(limits=c(-0.15,0.5))+
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
    annotate(geom="text", x=0.34, y=2.45, family = 'Arial',
             label='Dominance', size = 6)+
    labs(x = 'Standardized effects', y = '  ')+
    guides(color="none")+
    theme_regular())

(domin_pd_fd.a.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(domin_data.inla.md.a_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    geom_point(data=(domin_data.inla.md.a_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all'))),
               mapping = aes(x=mean,y=rowname), col = 'grey',
               size = 4.4)+
    geom_errorbar(data=domin_data.inla.md.a_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MFD.ab_all', 'MPD.ab_all')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper), col = 'grey',
                  width = 0.12, linewidth = 1.2)+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all'),
                     label = c(expression(MFD['ab']), expression(MPD['ab'])))+
    scale_x_continuous(limits=c(-1,0.15), breaks = seq(-1, 0, 0.5))+
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
    annotate(geom="text", x=-0.75, y=2.45, family = 'Arial',
             label='Dominance', size = 6)+
    labs(x = 'Standardized effects', y = '  ')+
    guides(color="none")+
    theme_regular())

(domin_md.a_ratio.all.varied.intercept.plot =
    ggplot(data=domin_data.inla.md.a_ratio_func_d.all_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 3)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all', 'MRFD.ab', 'MND.ab', 'demo_rate'),
                     label = c(expression(MFD['ab']), expression(MPD['ab']),
                               expression(RFD['ab']), expression(ND['ab']),
                               'Intrinsic growth rate'))+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_manual(values = c('purple',"#80defb", "#5a57fe","#fbd76c", "#f55756"),
                       name = ' ')+
    annotate(geom="text", x=1, y=5.5, family = 'Arial',
             label='Dominance')+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())

(domin_demo_rate.a.all.varied.intercept.plot =
    ggplot(data=(data=domin_data.inla.md.a_ratio_func_d.all_intercept %>% 
                   filter(rowname %in% c('demo_rate'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 3)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('demo_rate'),
                     label = c('Intrinsic growth rate'))+
    scale_x_continuous(limits=c(-1,1))+
    scale_color_manual(values = c('purple'),
                       name = ' ')+
    annotate(geom="text", x=0.6, y=1.5, family = 'Arial',
             label='Dominance')+
    labs(x = 'Standardized effects', y = '  ')+
    guides(color="none")+
    theme_regular())


# R2%
(domin_md.a.all.varied.intercept.R2.plot = 
    ggplot(data=domin_data.inla.md.a_func_d.all_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all', 'MRFD.ab', 'MND.ab'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
dominance_all_md.a.all = ggarrange(domin_md.a.all.varied.intercept.plot,
                                   domin_md.a.all.varied.intercept.R2.plot,
                                   widths=c(2,1.2))


### Effect size plot for mean nearest differences
# point + effect size
(domin_nnd_nrfd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(domin_data.inla.mnd_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MNRFD', 'MNND'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c('MNND' = 'MNND',
                                  "MNRFD" = "MNRFD"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[1],
                                 colors_4d[2]),
                      labels = c('MNND' = 'MNND',
                                 "MNRFD" = "MNRFD"),
                      name = ' ')+
    geom_point(data=(domin_data.inla.mnd_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MNRFD', 'MNND'))),
               mapping = aes(x=mean,y=rowname, color = significant),
               size = 4.4)+
    geom_errorbar(data=domin_data.inla.mnd_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MNRFD', 'MNND')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper, color = significant),
                  width = 0.12, linewidth = 1.2)+
    scale_color_manual(values = c('grey',
                                  'black'),
                       name = ' ')+
    ggnewscale::new_scale_color()+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNRFD', 'MNND'),
                     label = c(expression(RFD['nearest']),
                               expression(ND['nearest'])))+
    scale_x_continuous(limits=c(-0.15,0.5))+
    annotate(geom="text", x=0.35, y=2.45, family = 'Arial',
             label='Dominance', size = 6)+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

(domin_npd_nfd.all.varied.intercept.plot =
    ggplot()+
    geom_col(data=(domin_data.inla.mnd_conti_func_d.all_intercept %>% 
                     filter(rowname %in% c('MNFD_all', 'MNTD_all'))),
             mapping = aes(x=mean,y=rowname,fill=rowname),
             width = 0.4)+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MNTD_all" = "MNTD_all",
                                  "MNFD_all" = "MNFD_all"),
                       name = ' ')+
    scale_fill_manual(values = c(colors_4d[4],
                                 colors_4d[3]),
                      labels = c("MNTD_all" = "MNTD_all",
                                 "MNFD_all" = "MNFD_all"),
                      name = ' ')+
    geom_point(data=(domin_data.inla.mnd_conti_func_d.all_intercept %>% 
                       filter(rowname %in% c('MNFD_all', 'MNTD_all'))),
               mapping = aes(x=mean,y=rowname,color = significant),
               size = 4.4)+
    geom_errorbar(data=domin_data.inla.mnd_conti_func_d.all_intercept %>% 
                    filter(rowname %in% c('MNFD_all', 'MNTD_all')),
                  aes(x = mean,y = rowname,
                      xmin=lower,xmax=upper,color = significant),
                  width = 0.12, linewidth = 1.2)+
    scale_color_manual(values = c('black',
                                  'grey'),
                       name = ' ')+
    ggnewscale::new_scale_color()+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNFD_all','MNTD_all'),
                     label =c('NFD','NPD') )+
    scale_x_continuous(limits=c(-1,0.5), breaks = seq(-1, 0.5, 0.5))+
    annotate(geom="text", x=-0.69, y=2.45, family = 'Arial',
             label='Dominance', size = 6)+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

# R2%
(domin_mnd.all.varied.intercept.R2.plot = 
    ggplot(data=domin_data.inla.mnd_func_d.all_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MNFD_all','MNTD_all', 'MNRFD', 'MNND'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
dominance_all_mnd.all = ggarrange(domin_mnd.all.varied.intercept.plot,
                                  domin_mnd.all.varied.intercept.R2.plot,
                                  widths=c(2,1.2))



###### dominance predictive curves for md ####
### get dominance predict data for mnd
lincombs.data.domin.mnd = data.frame(mnd=seq(min(dat_dom_sp$mnd),max(dat_dom_sp$mnd),length=100),
                                     mlgfd=mean(dat_dom_sp$mlgfd),
                                     mpd_all = mean(dat_dom_sp$mpd_all),
                                     mfunc_d_all = mean(dat_dom_sp$mfunc_d_all))

lincombs.matrix.domin.mnd=model.matrix(~mnd+mlgfd+mpd_all+mfunc_d_all,
                                       data=lincombs.data.domin.mnd)
lincombs.matrix.domin.mnd=as.data.frame(lincombs.matrix.domin.mnd)
lincombs.domin.mnd=inla.make.lincombs(lincombs.matrix.domin.mnd)

inla.model_lincombs.domin.mnd = pglmm(domin ~ mnd+mlgfd+mpd_all+mfunc_d_all+(1|species) + 
                                        (1|f_p) + (1|field), data = dat_dom_sp,
                                      family = "binomial", cov_ranef = list(species = tree),
                                      bayes_options = list(control.compute = list(dic=T,
                                                                                  waic=T,
                                                                                  cpo=T,
                                                                                  config = TRUE),
                                                           quantiles=c(0.025,0.5,0.975),
                                                           lincomb=lincombs.domin.mnd,
                                                           control.predictor=list(compute=T)),
                                      bayes = T)

lincombs.posterior.domin.mnd = inla.model_lincombs.domin.mnd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.mnd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnd$predicted.value=unlist(lapply(lincombs.posterior.domin.mnd,
                                                      function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnd$lower=unlist(lapply(lincombs.posterior.domin.mnd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnd$upper=unlist(lapply(lincombs.posterior.domin.mnd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

save(lincombs.data.domin.mnd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance~mnd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnd.rdata")

(domin.mnd.partial.logistic=ggplot(data=lincombs.data.domin.mnd,
                                   aes(x=mnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1, linetype = 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd, y=domin),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x='Established-native niche difference',
         y="Dominance probability")+
    annotate(geom="text",x=c(0.33,0.33),y=c(0.80,0.70),
             label=c("italic(β)['E-N MND'] == -0.71",
                     "'95%CI' == '[-2.35, 0.94]'"),parse=T,size=3.5)
)


### get dominance predict data for mlgfd
lincombs.data.domin.mlgfd = data.frame(mlgfd=seq(min(dat_dom_sp$mlgfd),max(dat_dom_sp$mlgfd),length=100),
                                       mnd=mean(dat_dom_sp$mnd),
                                       mpd_all = mean(dat_dom_sp$mpd_all),
                                       mfunc_d_all = mean(dat_dom_sp$mfunc_d_all))

lincombs.matrix.domin.mlgfd=model.matrix(~mnd+mlgfd+mpd_all+mfunc_d_all,
                                         data=lincombs.data.domin.mlgfd)
lincombs.matrix.domin.mlgfd=as.data.frame(lincombs.matrix.domin.mlgfd)
lincombs.domin.mlgfd=inla.make.lincombs(lincombs.matrix.domin.mlgfd)

inla.model_lincombs.domin.mlgfd = pglmm(domin ~ mnd+mlgfd+mpd_all+mfunc_d_all+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_dom_sp,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.domin.mlgfd,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)

lincombs.posterior.domin.mlgfd = inla.model_lincombs.domin.mlgfd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mlgfd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mlgfd$predicted.value=unlist(lapply(lincombs.posterior.domin.mlgfd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mlgfd$lower=unlist(lapply(lincombs.posterior.domin.mlgfd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mlgfd$upper=unlist(lapply(lincombs.posterior.domin.mlgfd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mlgfd
save(lincombs.data.domin.mlgfd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mlgfd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd.rdata")
(domin.mlgfd.partial.logistic=ggplot(data=lincombs.data.domin.mlgfd,aes(x=mlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mlgfd, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native fitness difference', y='  ')+
    annotate(geom="text",x=c(0,0),y=c(0.80,0.70),
             label=c("italic(β)['E-N MRFD'] == 2.17",
                     "'95%CI' == '[0.60, 3.75]'"),
             parse=T,size=3.5)
)

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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mpd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mpd_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mpd_all.rdata")
(domin.mpd_all.partial.logistic=ggplot(data=lincombs.data.domin.mpd_all,aes(x=mpd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mpd_all, y=domin),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native phylogenetic difference', y='Dominance probability')+
    annotate(geom="text",x=c(260,260),y=c(0.80,0.70),
             label=c("italic(β)['E-N MPD'] == '-0.02'",
                     "'95%CI' == '[-0.04, 0.00]'"),
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mfunc_d_all.rdata")
(domin.mfunc_d_all.partial.logistic=ggplot(data=lincombs.data.domin.mfunc_d_all,
                                           aes(x=mfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1,
              linetype= 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfunc_d_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native functional difference', y='  ')+
    annotate(geom="text",x=c(0.3,0.3),y=c(0.80,0.70),
             label=c("italic(β)['E-N MFD'] == -0.47",
                     "'95%CI' == '[-9.97, 8.88]'"),
             parse=T,size=3.5)
)


###### Dominance predictive curves for mean nearest difference ####
### Get dominance predict data for mnnd
lincombs.data.domin.mnnd = data.frame(mnnd=seq(min(dat_dom_sp$mnnd),max(dat_dom_sp$mnnd),length=100),
                                      mnlgfd=mean(dat_dom_sp$mnlgfd),
                                      mntd_all = mean(dat_dom_sp$mntd_all),
                                      mnfunc_d_all = mean(dat_dom_sp$mnfunc_d_all))

lincombs.matrix.domin.mnnd=model.matrix(~mnnd+mnlgfd+mntd_all+mnfunc_d_all,
                                        data=lincombs.data.domin.mnnd)
lincombs.matrix.domin.mnnd=as.data.frame(lincombs.matrix.domin.mnnd)
lincombs.domin.mnnd=inla.make.lincombs(lincombs.matrix.domin.mnnd)

inla.model_lincombs.domin.mnnd = pglmm(domin ~ mnnd+mnlgfd+mntd_all+mnfunc_d_all+(1|species) + 
                                         (1|f_p) + (1|field), data = dat_dom_sp,
                                       family = "binomial", cov_ranef = list(species = tree),
                                       bayes_options = list(control.compute = list(dic=T,
                                                                                   waic=T,
                                                                                   cpo=T,
                                                                                   config = TRUE),
                                                            quantiles=c(0.025,0.5,0.975),
                                                            lincomb=lincombs.domin.mnnd,
                                                            control.predictor=list(compute=T)),
                                       bayes = T)

lincombs.posterior.domin.mnnd = inla.model_lincombs.domin.mnnd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnnd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.mnnd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnnd$predicted.value=unlist(lapply(lincombs.posterior.domin.mnnd,
                                                       function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnnd$lower=unlist(lapply(lincombs.posterior.domin.mnnd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnnd$upper=unlist(lapply(lincombs.posterior.domin.mnnd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

save(lincombs.data.domin.mnnd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance~mnnd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnnd.rdata")

(domin.mnnd.partial.logistic=ggplot(data=lincombs.data.domin.mnnd,
                                    aes(x=mnnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1, linetype = 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnnd, y=domin),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x='Established-native niche difference',
         y="Dominance probability")+
    annotate(geom="text",x=c(-0.25,-0.25),y=c(0.80,0.70),
             label=c("italic(β)['E-N MNND'] == -0.70",
                     "'95%CI' == '[-1.47, 0.08]'"),parse=T,size=3.5)
)


### get dominance predict data for mnlgfd
lincombs.data.domin.mnlgfd = data.frame(mnlgfd=seq(min(dat_dom_sp$mnlgfd),max(dat_dom_sp$mnlgfd),length=100),
                                        mnnd=mean(dat_dom_sp$mnnd),
                                        mntd_all = mean(dat_dom_sp$mntd_all),
                                        mnfunc_d_all = mean(dat_dom_sp$mnfunc_d_all))

lincombs.matrix.domin.mnlgfd=model.matrix(~mnnd+mnlgfd+mntd_all+mnfunc_d_all,
                                          data=lincombs.data.domin.mnlgfd)
lincombs.matrix.domin.mnlgfd=as.data.frame(lincombs.matrix.domin.mnlgfd)
lincombs.domin.mnlgfd=inla.make.lincombs(lincombs.matrix.domin.mnlgfd)

inla.model_lincombs.domin.mnlgfd = pglmm(domin ~ mnnd+mnlgfd+mntd_all+mnfunc_d_all+(1|species) + 
                                           (1|f_p) + (1|field), data = dat_dom_sp,
                                         family = "binomial", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = TRUE),
                                                              quantiles=c(0.025,0.5,0.975),
                                                              lincomb=lincombs.domin.mnlgfd,
                                                              control.predictor=list(compute=T)),
                                         bayes = T)

lincombs.posterior.domin.mnlgfd = inla.model_lincombs.domin.mnlgfd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnlgfd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnlgfd$predicted.value=unlist(lapply(lincombs.posterior.domin.mnlgfd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnlgfd$lower=unlist(lapply(lincombs.posterior.domin.mnlgfd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnlgfd$upper=unlist(lapply(lincombs.posterior.domin.mnlgfd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnlgfd
save(lincombs.data.domin.mnlgfd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnlgfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnlgfd
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnlgfd.rdata")
(domin.mnlgfd.partial.logistic=ggplot(data=lincombs.data.domin.mnlgfd,aes(x=mnlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnlgfd, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native fitness difference', y='  ')+
    annotate(geom="text",x=c(-0.25,-0.25),y=c(0.80,0.70),
             label=c("italic(β)['E-N MNRFD'] == 1.75",
                     "'95%CI' == '[0.63, 2.87]'"),
             parse=T,size=3.5)
)

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
inla.model_lincombs.domin.mntd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mntd_all$predicted.value=unlist(lapply(lincombs.posterior.domin.mntd_all,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mntd_all$lower=unlist(lapply(lincombs.posterior.domin.mntd_all,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mntd_all$upper=unlist(lapply(lincombs.posterior.domin.mntd_all,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mntd_all
save(lincombs.data.domin.mntd_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mntd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mntd_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mntd_all.rdata")
(domin.mntd_all.partial.logistic=ggplot(data=lincombs.data.domin.mntd_all,aes(x=mntd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mntd_all, y=domin),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native phylogenetic difference', y='Dominance probability')+
    annotate(geom="text",x=c(150,150),y=c(0.80,0.70),
             label=c("italic(β)['E-N MNTD'] == '-0.0054'",
                     "'95%CI' == '[-0.0103, -0.0005]'"),
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnfunc_d_all.rdata")
(domin.mnfunc_d_all.partial.logistic=ggplot(data=lincombs.data.domin.mnfunc_d_all,
                                            aes(x=mnfunc_d_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1,
              linetype= 3)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnfunc_d_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Established-native functional difference', y='  ')+
    annotate(geom="text",x=c(0.1,.1),y=c(0.80,0.70),
             label=c("italic(β)['E-N MNFD'] == -0.76",
                     "'95%CI' == '[-6.17, 4.66]'"),
             parse=T,size=3.5)
)



###### dominance predictive curves for md.ab ####
### get dominance predict data for mnd.ab
lincombs.data.domin.mnd.a = data.frame(mnd.a=seq(min(dat_dom_sp$mnd.a),max(dat_dom_sp$mnd.a),length=100),
                                       mlgfd.a=mean(dat_dom_sp$mlgfd.a),
                                       mpd.a_all = mean(dat_dom_sp$mpd.a_all),
                                       mfunc_d.a_all = mean(dat_dom_sp$mfunc_d.a_all))

lincombs.matrix.domin.mnd.a=model.matrix(~mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all,
                                         data=lincombs.data.domin.mnd.a)
lincombs.matrix.domin.mnd.a=as.data.frame(lincombs.matrix.domin.mnd.a)
lincombs.domin.mnd.a=inla.make.lincombs(lincombs.matrix.domin.mnd.a)

inla.model_lincombs.domin.mnd.a = pglmm(domin ~ mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_dom_sp,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,
                                                                                    waic=T,
                                                                                    cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.domin.mnd.a,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)

lincombs.posterior.domin.mnd.a = inla.model_lincombs.domin.mnd.a$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnd.a$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.mnd.a$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnd.a$predicted.value=unlist(lapply(lincombs.posterior.domin.mnd.a,
                                                        function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnd.a$lower=unlist(lapply(lincombs.posterior.domin.mnd.a,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnd.a$upper=unlist(lapply(lincombs.posterior.domin.mnd.a,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

save(lincombs.data.domin.mnd.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnd.a.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance~mnd.ab
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mnd.a.rdata")
(domin.mnd.a.partial.logistic=ggplot(data=lincombs.data.domin.mnd.a,
                                     aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1,
              linetype = 3)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd.a, y=domin),shape=1,
               color = colors_4d[1],
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x='Niche difference (ND)',
         y="Dominance probability")+
    annotate(geom="text",x=c(0.33,0.33),y=c(0.80,0.70),
             label=c("italic(β)['MND'[ab]] == -0.03",
                     "'95%CI' == '[-1.34, 1.29]'"),parse=T,size=3.5)+
    theme_regular()
)


### get dominance predict data for mlgfd.ab
lincombs.data.domin.mlgfd.a = data.frame(mlgfd.a=seq(min(dat_dom_sp$mlgfd.a),max(dat_dom_sp$mlgfd.a),length=100),
                                         mnd.a=mean(dat_dom_sp$mnd.a),
                                         mpd.a_all = mean(dat_dom_sp$mpd.a_all),
                                         mfunc_d.a_all = mean(dat_dom_sp$mfunc_d.a_all))

lincombs.matrix.domin.mlgfd.a=model.matrix(~mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all,
                                           data=lincombs.data.domin.mlgfd.a)
lincombs.matrix.domin.mlgfd.a=as.data.frame(lincombs.matrix.domin.mlgfd.a)
lincombs.domin.mlgfd.a=inla.make.lincombs(lincombs.matrix.domin.mlgfd.a)

inla.model_lincombs.domin.mlgfd.a = pglmm(domin ~ mnd.a+mlgfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_dom_sp,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.domin.mlgfd.a,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)

lincombs.posterior.domin.mlgfd.a = inla.model_lincombs.domin.mlgfd.a$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mlgfd.a$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mlgfd.a$predicted.value=unlist(lapply(lincombs.posterior.domin.mlgfd.a,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mlgfd.a$lower=unlist(lapply(lincombs.posterior.domin.mlgfd.a,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mlgfd.a$upper=unlist(lapply(lincombs.posterior.domin.mlgfd.a,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mlgfd.a
save(lincombs.data.domin.mlgfd.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd.a.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mlgfd.ab
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mlgfd.a.rdata")
(domin.mlgfd.a.partial.logistic=ggplot(data=lincombs.data.domin.mlgfd.a,aes(x=mlgfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mlgfd.a, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Relative fitness difference (RFD)', y='  ')+
    annotate(geom="text",x=c(0,0),y=c(0.80,0.70),
             label=c("italic(β)['MRFD'[ab]] == 2.49", "'95%CI' == '[1.17, 3.80]'"),
             parse=T,size=3.5)+
    theme_regular()
)

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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mpd.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mpd.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mpd.a_all.rdata")
(domin.mpd.a_all.partial.logistic=ggplot(data=lincombs.data.domin.mpd.a_all,aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mpd.a_all, y=domin),
               color = colors_4d[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Phylogenetic difference (PD)', y='Dominance probability')+
    annotate(geom="text",x=c(250,250),y=c(0.80,0.70),
             label=c("italic(β)['MPD'[ab]] == '-0.01'", "'95%CI' == '[-0.02, -0.01]'"),
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mfunc_d.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfunc_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.domin.mfunc_d.a_all.rdata")
(domin.mfunc_d.a_all.partial.logistic=ggplot(data=lincombs.data.domin.mfunc_d.a_all,
                                             aes(x=mfunc_d.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[4],alpha=0.2)+
    geom_line(color=colors_4d[4],size=1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfunc_d.a_all, y=domin),
               color = colors_4d[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Functional difference (FD)', y='  ')+
    annotate(geom="text",x=c(0.25,0.25),y=c(0.80,0.70),
             label=c("italic(β)['MFD'[ab]] == '11.99'", "'95%CI' == '[8.44, 15.54]'"),
             parse=T,size=3.5)+
    theme_regular()
)



#### Fig.1：Demonstration of how Weighted Mean ND/RFD affects successful colonization and dominance of species, including presentation of raw data and single & multi-regression results ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Adjust grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height  =  1 / num_rows

# Create a ggdraw object with all plots
Fig.1.all = ggdraw() +
  draw_plot(estab.mnd.a_single.logistic, x = 0, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mlgfd.a_single.logistic, x = plot_width, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab_nd_rfd.a.all.varied.intercept.plot, x = 2 * plot_width, y = 0.5,
            width = plot_width,height = plot_height) +
  draw_plot(domin.mnd.a_single.logistic, x = 0, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin_nd_rfd.a.all.varied.intercept.plot, x = 2 * plot_width, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mlgfd.a_single.logistic, x = plot_width, y = 0.05, 
            width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01,
          0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01),
    y = c(1, 1, 1, 0.55, 0.55, 0.55),
    hjust = 0, vjust = 1.1, size = 18
  )

Fig.1.all
emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.1.all.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.1.all
dev.off() #turn off device and finalize file



#### Fig.2：Demonstration of how Weighted Mean PD/FD affects successful colonization and dominance of species, including presentation of raw data and single & multi-regression results ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Adjust grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height  =  1 / num_rows

# Create a ggdraw object with all plots
Fig.2.all = ggdraw() +
  draw_plot(estab.mpd.a_single.logistic, x = 0, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mconti_func_d.a_single.logistic, x = plot_width, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab_pd_fd.a.all.varied.intercept.plot, x = 2 * plot_width, y = 0.5,
            width = plot_width,height = plot_height) +
  draw_plot(domin.mpd.a_single.logistic, x = 0, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mconti_func_d.a_single.logistic, x = plot_width, y = 0.05, 
            width = plot_width, height = plot_height) +
  draw_plot(domin_pd_fd.a.all.varied.intercept.plot, x = 2 * plot_width, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01,
          0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01),
    y = c(1, 1, 1, 0.55, 0.55, 0.55),
    hjust = 0, vjust = 1.1, size = 18
  )

Fig.2.all
emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.2.all.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.2.all
dev.off() #turn off device and finalize file




#### Fig.S5：Demonstration of how Mean ND/RFD affects successful colonization and dominance of species ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Adjust grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height  =  1 / num_rows

# Create a ggdraw object with all plots
Fig.S5.all = ggdraw() +
  draw_plot(estab.mnd_single.logistic, x = 0, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mlgfd_single.logistic, x = plot_width, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab_nd_rfd.all.varied.intercept.plot, x = (2 * plot_width)+0.03, y = 0.5,
            width = plot_width-0.03,height = plot_height) +
  draw_plot(domin.mnd_single.logistic, x = 0, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mlgfd_single.logistic, x = plot_width, y = 0.05, 
            width = plot_width, height = plot_height) +
  draw_plot(domin_nd_rfd.all.varied.intercept.plot, x = (2 * plot_width)+0.03, y = 0.05,
            width = plot_width-0.03, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01,
          0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01),
    y = c(1, 1, 1, 0.55, 0.55, 0.55),
    hjust = 0, vjust = 1.1, size = 18
  )


#Fig.S5.all
emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S5.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S5.all
dev.off() #turn off device and finalize file




#### Fig.S6：Demonstration of how Mean PD/FD affects successful colonization and dominance of species ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Adjust grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height  =  1 / num_rows

# Create a ggdraw object with all plots
Fig.S6.all = ggdraw() +
  draw_plot(estab.mpd_single.logistic, x = 0, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mconti_func_d_single.logistic, x = plot_width, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab_pd_fd.all.varied.intercept.plot, x = 2 * plot_width+0.03, y = 0.5,
            width = plot_width-0.03,height = plot_height) +
  draw_plot(domin.mpd_single.logistic, x = 0, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mconti_func_d_single.logistic, x = plot_width, y = 0.05, 
            width = plot_width, height = plot_height) +
  draw_plot(domin_pd_fd.all.varied.intercept.plot, x = 2 * plot_width+0.03, y = 0.05,
            width = plot_width-0.03, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01,
          0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01),
    y = c(1, 1, 1, 0.55, 0.55, 0.55),
    hjust = 0, vjust = 1.1, size = 18
  )


#Fig.S6.all
emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S6.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S6.all
dev.off() #turn off device and finalize file



#### Fig.：Demonstration of how Mean Nearest ND/RFD affects successful colonization and dominance of species ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Adjust grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height  =  1 / num_rows

# Create a ggdraw object with all plots
Fig.S5.all = ggdraw() +
  draw_plot(estab.mnnd_single.logistic, x = 0, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mnlgfd_single.logistic, x = plot_width-0.007, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab_nnd_nrfd.all.varied.intercept.plot,
            x = 2 * plot_width-0.003, y = 0.5,
            width = plot_width+0.01,height = plot_height) +
  draw_plot(domin.mnnd_single.logistic, x = 0, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mnlgfd_single.logistic, x = plot_width-0.007, y = 0.05, 
            width = plot_width, height = plot_height) +
  draw_plot(domin_nnd_nrfd.all.varied.intercept.plot,
            x = 2 * plot_width-0.003, y = 0.05,
            width = plot_width+0.01, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, (1 * plot_width)-0.001+0.01, (2 * plot_width)-0.008+0.01,
          0.01, (1 * plot_width)-0.001+0.01, (2 * plot_width)-0.008+0.01),
    y = c(1, 1, 1, 0.55, 0.55, 0.55),
    hjust = 0, vjust = 1.1, size = 18
  )


#Fig.S5.all
emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S5.all.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S5.all
dev.off() #turn off device and finalize file




#### Fig.：Demonstration of how Mean Nearest PD/FD affects successful colonization and dominance of species ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Adjust grid settings
num_rows  =  2
num_cols  =  3
plot_width  =  1 / num_cols
plot_height  =  1 / num_rows

# Create a ggdraw object with all plots
Fig.S6.all = ggdraw() +
  draw_plot(estab.mntd_single.logistic, x = 0, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mnconti_func_d_single.logistic, x = plot_width, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab_npd_nfd.all.varied.intercept.plot, x = 2 * plot_width+0.03, y = 0.5,
            width = plot_width-0.03,height = plot_height) +
  draw_plot(domin.mntd_single.logistic, x = 0, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mnconti_func_d_single.logistic, x = plot_width, y = 0.05, 
            width = plot_width, height = plot_height) +
  draw_plot(domin_npd_nfd.all.varied.intercept.plot, x = 2 * plot_width+0.03, y = 0.05,
            width = plot_width-0.03, height = plot_height) +
  draw_plot_label(
    label = c("a", "b", "c", "d", "e", "f"),
    x = c(0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01,
          0.01, 1 * plot_width+0.01, (2 * plot_width)+0.01),
    y = c(1, 1, 1, 0.55, 0.55, 0.55),
    hjust = 0, vjust = 1.1, size = 18
  )


#Fig.S6.all
emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig.S6.all.emf',
    width = 25.2*1.1, height = 16.2*1.1, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S6.all
dev.off() #turn off device and finalize file


#### SEM: Success ~ mND + mRFD ~ mfd + mpd for all species ####
### Original model ###
require(glmmTMB)
library(piecewiseSEM)
require(optimx)
require(dplyr)

### SEM for establishment
setwd("D:/R projects/BSS")
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/dat_suc_sp.rdata')

numcols = grep("^m.",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)


estab_sem1_md.a_all = psem(
  glmmTMB(estab~mnd.a+mlgfd.a
          +mpd.a_all
          +mconti_func_d.a_all
          #+(1|species)
          +(1|field)
          #+(1|f_p)
          , family=binomial, data=dat_suc_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a_all+
            #mconti_func_d.a_all+
            #(1|species)+
            #+(1|field)
            (1|f_p)
          , family=gaussian, data = dat_suc_sps),
  glmmTMB(mlgfd.a~mpd.a_all+
            mconti_func_d.a_all+
            #(1|species)+
            (1|field)
          +(1|f_p)
          , family=gaussian, data=dat_suc_sps),
  mlgfd.a %~~% mnd.a,
  #estab %~~% mpd,
  #estab %~~% mconti_func_d,
  data=dat_suc_sps)
summary(estab_sem1_md.a_all)

### SEM for Dominance ###
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c("establish", "dominant"))

domin_sem1_md.a_all = psem(
  glmmTMB(domin~mnd.a+
            mlgfd.a+
            mpd.a_all+
            mconti_func_d.a_all+
            #(1|species)+
            (1|field)
          +(1|f_p)
          , family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a_all+
            mconti_func_d.a_all+
            #(1|species)+
            (1|field)
          #+(1|f_p)
          , family=gaussian, data = dat_dom_sps),
  glmmTMB(mlgfd.a~mpd.a_all+
            #mconti_func_d.a_all+
            (1|species)+
            (1|field)
          #+(1|f_p)
          , family=gaussian, data=dat_dom_sps),
  mlgfd.a %~~% mnd.a,
  #domin %~~% mpd,
  #domin %~~% mconti_func_d,
  data=dat_dom_sps)
summary(domin_sem1_md.a_all)

domin_sem1_mnd_all = psem(
  glmmTMB(domin~mnnd+mnlgfd+
            mntd_all+
            mnconti_func_d_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnnd~mntd_all+
            #mnconti_func_d_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_dom_sps),
  glmmTMB(mnlgfd~#mntd_all+
            mnconti_func_d_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_dom_sps),
  mnlgfd %~~% mnnd,
  #domin %~~% mpd,
  #domin %~~% mconti_func_d,
  data=dat_dom_sps)
summary(domin_sem1_mnd_all)


#--------------Collect all the key parameters of SEM as the supplementary table
FisherC.all_mpd = rbind(fisherC(sem1_mpd),fisherC(sem2_mpd),fisherC(sem3_mpd),fisherC(sem4_mpd),fisherC(sem5_mpd))

AIC.K_mpd = rbind(infCrit(sem1_mpd),infCrit(sem2_mpd),infCrit(sem3_mpd),infCrit(sem4_mpd),infCrit(sem5_mpd))
Delta.AIC_mpd = dplyr::rename(AIC.K_mpd[1]-min(AIC.K_mpd[1]),Delta.AIC=AIC)

AIC.weight_mpd = c(AIC(sem1_mpd),AIC(sem2_mpd),AIC(sem3_mpd),AIC(sem4_mpd),AIC(sem5_mpd)) %>%
  MuMIn::Weights()

sem.table_mpd=cbind(FisherC.all_mpd[-2],
                    AIC.K_mpd[c(4,1)],
                    Delta.AIC_mpd,
                    AIC.weight_mpd)%>%round(2)
sem.table_mpd


##### Draw plot for md 
### Establishment corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_estab_sem1_md_all = sqrt(predict(estab_sem1_md_all[[1]],re.form=NA)%>%var+sum(as.numeric(unlist(VarCorr(estab_sem1_md_all[[1]]))))+pi^2/3)
sd.y_estab_sem1_md_all
std.estab_sem1_md_all=unlist(fixef(estab_sem1_md_all[[1]]))[-1]*
  c(sd(dat_suc_sps$mpd_all),
    sd(dat_suc_sps$mfunc_d_all),
    sd(dat_suc_sps$mnd),
    sd(dat_suc_sps$mlgfd))/sd.y_estab_sem1_md_all
std.estab_sem1_md_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_estab_sem1_md_all=coefs(estab_sem1_md_all)# original coefficients of piecewiseSEM
coefs_estab_sem1_md_all[1:4,8]=std.estab_sem1_md_all

### Dominance corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_domin_sem1_md_all=sqrt(predict(domin_sem1_md_all[[1]],re.form=NA)%>%var+sum(as.numeric(unlist(VarCorr(domin_sem1_md_all[[1]]))))+pi^2/3)
sd.y_domin_sem1_md_all
std.domin_sem1_md_all=unlist(fixef(domin_sem1_md_all[[1]]))[-1]*
  c(sd(dat_dom_sps$mpd_all),
    sd(dat_dom_sps$mfunc_d_all),
    sd(dat_dom_sps$mnd),
    sd(dat_dom_sps$mlgfd))/sd.y_domin_sem1_md_all
std.domin_sem1_md_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_domin_sem1_md_all=coefs(domin_sem1_md_all)# original coefficients of piecewiseSEM
coefs_domin_sem1_md_all[1:4,8]=std.domin_sem1_md_all
coefs_domin_sem1_md_all

### Calculate direct and indirect effects of SEM following Xu Meng (2022) GCB
### Estab
sem_estab_effects_md_all = data.frame(predictor = rep(unique(coefs_estab_sem1_md_all$Predictor)[1:4], 3),
                                      type = rep(c('direct', 'indirect', 'total'), each = 4),
                                      value = rep(0, 12),
                                      significant = rep('yes', 12))
coefs_estab_sem1_md_all = coefs_estab_sem1_md_all[-nrow(coefs_domin_sem1_md_all),]
coefs_estab_sem1_md_all_c = coefs_estab_sem1_md_all[coefs_estab_sem1_md_all$P.Value < 0.05,]
for (i in 1:length(unique(coefs_estab_sem1_md_all_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_estab_sem1_md_all_c$Predictor)[i]
  dat = coefs_estab_sem1_md_all_c[coefs_estab_sem1_md_all_c$Predictor == predictor_1,]
  dat_2 = coefs_estab_sem1_md_all_c[coefs_estab_sem1_md_all_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_estab_effects_md_all[sem_estab_effects_md_all$predictor == predictor_1&
                               sem_estab_effects_md_all$type == 'direct',]$value = dat_3$Std.Estimate
    sem_estab_effects_md_all[sem_estab_effects_md_all$predictor == predictor_1&
                               sem_estab_effects_md_all$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    sem_estab_effects_md_all[sem_estab_effects_md_all$predictor == predictor_1&
                               sem_estab_effects_md_all$type == 'direct',]$value = dat_3[dat_3$Response == 'estab' &
                                                                                           dat_3$Predictor == predictor_1,]$Std.Estimate
    sem_estab_effects_md_all[sem_estab_effects_md_all$predictor == predictor_1&
                               sem_estab_effects_md_all$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'estab',]$Std.Estimate*dat_3[dat_3$Response == 'estab' &
                                                                                                                                                dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_estab_effects_md_all[sem_estab_effects_md_all$predictor == predictor_1&
                               sem_estab_effects_md_all$type == 'total',]$value  = sem_estab_effects_md_all[sem_estab_effects_md_all$predictor == predictor_1&
                                                                                                              sem_estab_effects_md_all$type == 'direct',]$value + sem_estab_effects_md_all[sem_estab_effects_md_all$predictor == predictor_1&
                                                                                                                                                                                             sem_estab_effects_md_all$type == 'indirect',]$value
  }
}

### Domin
sem_domin_effects_md_all = data.frame(predictor = rep(unique(coefs_domin_sem1_md_all$Predictor)[1:4], 3),
                                      type = rep(c('direct', 'indirect', 'total'), each = 4),
                                      value = rep(0, 12))
coefs_domin_sem1_md_all = coefs_domin_sem1_md_all[-nrow(coefs_domin_sem1_md_all),]
coefs_domin_sem1_md_all_c = coefs_domin_sem1_md_all[coefs_domin_sem1_md_all$P.Value < 0.05,]
for (i in 1:length(unique(coefs_domin_sem1_md_all_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_domin_sem1_md_all_c$Predictor)[i]
  dat = coefs_domin_sem1_md_all_c[coefs_domin_sem1_md_all_c$Predictor == predictor_1,]
  dat_2 = coefs_domin_sem1_md_all_c[coefs_domin_sem1_md_all_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_domin_effects_md_all[sem_domin_effects_md_all$predictor == predictor_1&
                               sem_domin_effects_md_all$type == 'direct',]$value = dat_3$Std.Estimate
    sem_domin_effects_md_all[sem_domin_effects_md_all$predictor == predictor_1&
                               sem_domin_effects_md_all$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    sem_domin_effects_md_all[sem_domin_effects_md_all$predictor == predictor_1&
                               sem_domin_effects_md_all$type == 'direct',]$value = dat_3[dat_3$Response == 'domin' &
                                                                                           dat_3$Predictor == predictor_1,]$Std.Estimate
    sem_domin_effects_md_all[sem_domin_effects_md_all$predictor == predictor_1&
                               sem_domin_effects_md_all$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'domin',]$Std.Estimate*dat_3[dat_3$Response == 'domin' &
                                                                                                                                                dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_domin_effects_md_all[sem_domin_effects_md_all$predictor == predictor_1&
                               sem_domin_effects_md_all$type == 'total',]$value  = sem_domin_effects_md_all[sem_domin_effects_md_all$predictor == predictor_1&
                                                                                                              sem_domin_effects_md_all$type == 'direct',]$value + sem_domin_effects_md_all[sem_domin_effects_md_all$predictor == predictor_1&
                                                                                                                                                                                             sem_domin_effects_md_all$type == 'indirect',]$value
  }
}

#-----------Draw plots for SEMs' relative total effects, direct and indirect effects 
### Estab
require(ggthemes)
require(ggplot2)
require(ggpubr)

sem_estab_effects_md_all_total = sem_estab_effects_md_all %>% filter(type == 'total')
sem_estab_effects_md_all.percent = sem_estab_effects_md_all_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_estab_effects_md_all.longer = sem_estab_effects_md_all %>% filter(type != 'total' 
                                                                      #& value != 0
)

# Relative total effect 
sem_estab_effects_md_all.percent.plot = 
  ggplot(data=sem_estab_effects_md_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mfunc_d_all", "mpd_all", "mlgfd", "mnd"),
                   labels=c("I-N MFD_all", "I-N MPD_all", "I-N MRFD", "I-N MND"),
                   position="right")+
  scale_fill_stata()+# color palatte in ggthemes
  xlim(0,0.8)+
  xlab("")+
  #theme_test()+
  theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))+
  guides(fill="none")

# Direct and indirect effects
sem_estab_effects_md_all.dir_indir.mpd_all.plot = 
  ggplot(data=sem_estab_effects_md_all.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects","indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mfunc_d_all", "mpd_all", "mlgfd", "mnd"),
                   position="right")+
  scale_x_continuous(breaks = c(seq(-0.4, 0.1, 0.2)))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_estab_md_all_space = ggplot(NULL)+theme_void()

# Merge three plots
sem_estab_md_all_plot_all = ggarrange(sem_estab_md_all_space,
                                      sem_estab_effects_md_all.dir_indir.mpd_all.plot,
                                      sem_estab_effects_md_all.percent.plot,
                                      nrow = 1, ncol = 3,
                                      labels = c('(a)', '', ''),
                                      vjust = 1.8)


ggsave(plot = sem_estab_md_all_plot_all,
       "results/figures_ages1_35_top50_equal_interval_bh_partialb/all_species_d/sem_estab_md_all_plot_all.svg",
       width=9,height=3)

### domin
sem_domin_effects_md_all_total = sem_domin_effects_md_all %>% filter(type == 'total')
sem_domin_effects_md_all.percent = sem_domin_effects_md_all_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_domin_effects_md_all.longer = sem_domin_effects_md_all %>% filter(type != 'total' 
                                                                      #&value != 0
)

# Relative total effect
sem_domin_effects_md_all.percent.plot = 
  ggplot(data=sem_domin_effects_md_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mfunc_d_all", "mpd_all", "mlgfd", "mnd"),
                   labels=c("E-N MFD_all", "E-N MPD_all", "E-N MRFD", "E-N MND"),
                   position="right")+
  scale_fill_stata()+# color palatte in ggthemes
  xlim(0,0.8)+
  xlab("")+
  #theme_test()+
  theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))+
  guides(fill="none")

# Direct and indirect effects
sem_domin_effects_md_all.dir_indir.mpd_all.plot = 
  ggplot(data=sem_domin_effects_md_all.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects",
                                          "indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mfunc_d_all", "mpd_all", "mlgfd", "mnd"),
                   position="right")+
  scale_x_continuous(breaks = seq(-0.1, 0.1, 0.1),
                     limits = c(-0.11, 0.15))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_domin_md_all_space = ggplot()+
  facet_wrap(~"Dominance")+
  theme_test() + 
  theme(strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.015,1.2,0.2),units="lines"))

# Merge three plots
sem_domin_md_all_plot_all = ggarrange(sem_domin_md_all_space,
                                      sem_domin_effects_md_all.dir_indir.mpd_all.plot,
                                      sem_domin_effects_md_all.percent.plot,
                                      nrow = 1, ncol = 3,
                                      labels = c('(b)', '', ''),
                                      vjust = 1.8)


ggsave(plot = sem_domin_md_all_plot_all,
       "results/figures_ages1_35_top50_equal_interval_bh_partialb/all_species_d/sem_domin_md_all_plot_all.svg",
       width=9,height=3)


##### Draw plot for md.a ##### 
require(ggh4x)
### The followed code is adapted from Xu_et_al_2022_GCB
### Establishment corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_estab_sem1_md.a_all = sqrt(predict(estab_sem1_md.a_all[[1]],
                                        re.form=NA)%>%var+sum(unlist(VarCorr(estab_sem1_md.a_all[[1]])))+pi^2/3)
sd.y_estab_sem1_md.a_all
std.estab_sem1_md.a_all=unlist(fixef(estab_sem1_md.a_all[[1]]))[-1]*
  c(sd(dat_suc_sps$mpd.a_all),
    sd(dat_suc_sps$mconti_func_d.a_all),
    sd(dat_suc_sps$mnd),
    sd(dat_suc_sps$mlgfd))/sd.y_estab_sem1_md.a_all
std.estab_sem1_md.a_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_estab_sem1_md.a_all=coefs(estab_sem1_md.a_all)# original coefficients of piecewiseSEM
#coefs_estab_sem1_md.a_all[1:4,8]=std.estab_sem1_md.a_all

### Dominance corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_domin_sem1_md.a_all=sqrt(predict(domin_sem1_md.a_all[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(domin_sem1_md.a_all[[1]])))+pi^2/3)
sd.y_domin_sem1_md.a_all
std.domin_sem1_md.a_all=unlist(fixef(domin_sem1_md.a_all[[1]]))[-1]*
  c(sd(dat_dom_sps$mpd.a_all),
    sd(dat_dom_sps$mconti_func_d.a_all),
    sd(dat_dom_sps$mnd),
    sd(dat_dom_sps$mlgfd))/sd.y_domin_sem1_md.a_all
std.domin_sem1_md.a_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_domin_sem1_md.a_all=coefs(domin_sem1_md.a_all)# original coefficients of piecewiseSEM
#coefs_domin_sem1_md.a_all[1:4,8]=std.domin_sem1_md.a_all
coefs_domin_sem1_md.a_all

### Calculate direct and indirect effects of SEM following Xu Meng (2022) GCB
### Estab
coefs_estab_sem1_md.a_all = coefs_estab_sem1_md.a_all[-nrow(coefs_estab_sem1_md.a_all),]
sem_estab_effects_md.a_all = data.frame(predictor = rep(c('mnd.a', 'mlgfd.a',
                                                          'mpd.a_all', 'mconti_func_d.a_all'), 3),
                                        type = rep(c('direct', 'indirect', 'total'),
                                                   each = 4),
                                        value = 0,
                                        significant = 'yes')
coefs_estab_sem1_md.a_all_c = coefs_estab_sem1_md.a_all[coefs_estab_sem1_md.a_all$P.Value < 0.05,]
for (i in 1:length(unique(coefs_estab_sem1_md.a_all_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_estab_sem1_md.a_all_c$Predictor)[i]
  dat = coefs_estab_sem1_md.a_all_c[coefs_estab_sem1_md.a_all_c$Predictor == predictor_1,]
  dat_2 = coefs_estab_sem1_md.a_all_c[coefs_estab_sem1_md.a_all_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                 sem_estab_effects_md.a_all$type == 'direct',]$value = dat_3$Std.Estimate
    sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                 sem_estab_effects_md.a_all$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    if (nrow(dat_3[dat_3$Response == 'estab' &
                   dat_3$Predictor == predictor_1,]) > 0) {
      sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                   sem_estab_effects_md.a_all$type == 'direct',]$value = dat_3[dat_3$Response == 'estab' &
                                                                                                 dat_3$Predictor == predictor_1,]$Std.Estimate
    }
    
    sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                 sem_estab_effects_md.a_all$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'estab',]$Std.Estimate*dat_3[dat_3$Response == 'estab' &
                                                                                                                                                    dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                 sem_estab_effects_md.a_all$type == 'total',]$value  = sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                                                                                                    sem_estab_effects_md.a_all$type == 'direct',]$value + sem_estab_effects_md.a_all[sem_estab_effects_md.a_all$predictor == predictor_1&
                                                                                                                                                                                                       sem_estab_effects_md.a_all$type == 'indirect',]$value
  }
}

### Domin
coefs_domin_sem1_md.a_all = coefs_domin_sem1_md.a_all[-nrow(coefs_domin_sem1_md.a_all),]
sem_domin_effects_md.a_all = data.frame(predictor = rep(c('mnd.a', 'mlgfd.a',
                                                          'mpd.a_all', 'mconti_func_d.a_all'), 3),
                                        type = rep(c('direct', 'indirect', 'total'),
                                                   each = 4),
                                        value = 0,
                                        significant = 'yes')

coefs_domin_sem1_md.a_all_c = coefs_domin_sem1_md.a_all[coefs_domin_sem1_md.a_all$P.Value < 0.05,]
for (i in 1:length(unique(coefs_domin_sem1_md.a_all_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_domin_sem1_md.a_all_c$Predictor)[i]
  dat = coefs_domin_sem1_md.a_all_c[coefs_domin_sem1_md.a_all_c$Predictor == predictor_1,]
  dat_2 = coefs_domin_sem1_md.a_all_c[coefs_domin_sem1_md.a_all_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'direct',]$value = dat_3$Std.Estimate
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'direct',]$value = dat_3[dat_3$Response == 'domin' &
                                                                                               dat_3$Predictor == predictor_1,]$Std.Estimate
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'domin',]$Std.Estimate*dat_3[dat_3$Response == 'domin' &
                                                                                                                                                    dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                 sem_domin_effects_md.a_all$type == 'total',]$value  = sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                                                                                                    sem_domin_effects_md.a_all$type == 'direct',]$value + sem_domin_effects_md.a_all[sem_domin_effects_md.a_all$predictor == predictor_1&
                                                                                                                                                                                                       sem_domin_effects_md.a_all$type == 'indirect',]$value
  }
}

#-----------Draw plots for SEMs' relative total effects, direct and indirect effects 
### Estab
require(ggthemes)
require(ggplot2)
require(devEMF)
sem_estab_effects_md.a_all_total = sem_estab_effects_md.a_all %>% filter(type == 'total')
sem_estab_effects_md.a_all.percent = sem_estab_effects_md.a_all_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_estab_effects_md.a_all.longer = sem_estab_effects_md.a_all %>% filter(type != 'total' 
                                                                          #& value != 0
)

# Relative total effect 
sem_estab_effects_md.a_all.percent.plot = 
  ggplot(data=sem_estab_effects_md.a_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   labels=#c(bquote('I-N MFD'[ab]),bquote('I-N MPD'[ab]),bquote('I-N MRFD'[ab]),bquote('I-N MND'[ab]))
                     c(expression(MFD[ab]),
                       expression(MPD[ab]),
                       expression(RFD[ab]),
                       expression(ND[ab])) ,
                   position="right")+
  scale_fill_manual(values = c("#ccb974ff", "#dd8452ff", "#4c72b0ff", "#55a868ff"))+# color palatte in ggthemes
  xlim(0,0.5)+
  xlab("")+
  #theme_test()+
  theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))+
  guides(fill="none")
sem_estab_effects_md.a_all.percent.plot

# Direct and indirect effects
sem_estab_effects_md.a_all.dir_indir.mpd_all.plot = 
  ggplot(data=sem_estab_effects_md.a_all.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects","indirect"="Indirect effects")),
             scales = 'free')+
  ggh4x::facetted_pos_scales(
    x = list(
      `direct` = scale_x_continuous(limits = c(-0.44, 0.44), breaks = seq(-0.4, 0.4, 0.2)),
      `indirect` = scale_x_continuous(limits = c(-0.14, 0.14), breaks = seq(-0.1, 0.1, 0.1))  
    )
  )+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   position="right")+
  scale_fill_manual(values = c("#ccb974ff", "#dd8452ff", "#4c72b0ff", "#55a868ff"))+# color palatte in ggthemes
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))
sem_estab_effects_md.a_all.dir_indir.mpd_all.plot

### Space for SEM
sem_estab_md.a_all_space1 = ggplot(NULL)+theme_void()
sem_estab_md.a_all_space2 = ggplot(NULL)+theme_void()

### domin
sem_domin_effects_md.a_all_total = sem_domin_effects_md.a_all %>% filter(type == 'total')
sem_domin_effects_md.a_all.percent = sem_domin_effects_md.a_all_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_domin_effects_md.a_all.longer = sem_domin_effects_md.a_all %>% filter(type != 'total' 
                                                                          #&value != 0
)

# Relative total effect
sem_domin_effects_md.a_all.percent.plot = 
  ggplot(data=sem_domin_effects_md.a_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   labels=#c(bquote('E-N MFD'[ab]),bquote('E-N MPD'[ab]),bquote('E-N MRFD'[ab]),bquote('E-N MND'[ab]))
                     c(expression(MFD[ab]),
                       expression(MPD[ab]),
                       expression(RFD[ab]),
                       expression(ND[ab])),
                   position="right")+
  scale_fill_manual(values = c("#ccb974ff", "#dd8452ff", "#4c72b0ff", "#55a868ff"))+# color palatte in ggthemes
  xlim(0,0.72)+
  xlab("")+
  #theme_test()+
  theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))+
  guides(fill="none")

# Direct and indirect effects
sem_domin_effects_md.a_all.dir_indir.mpd_all.plot = 
  ggplot(data=sem_domin_effects_md.a_all.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects",
                                          "indirect"="Indirect effects")),
             scales = 'free')+
  ggh4x::facetted_pos_scales(
    x = list(
      `direct` = scale_x_continuous(limits = c(-0.44, 0.44), breaks = seq(-0.4, 0.4, 0.2)),
      `indirect` = scale_x_continuous(limits = c(-0.24, 0.24), breaks = seq(-0.2, 0.2, 0.1))  
    )
  )+
  scale_y_discrete(limits=c("mconti_func_d.a_all", "mpd.a_all", "mlgfd.a", "mnd.a"),
                   position="right")+
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2),
                     limits = c(-0.4, 0.4))+
  scale_fill_manual(values = c("#ccb974ff", "#dd8452ff", "#4c72b0ff", "#55a868ff"))+# color palatte in ggthemes
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_domin_md.a_all_space1 = ggplot(NULL)+theme_void()
sem_domin_md.a_all_space2 = ggplot(NULL)+theme_void()

# Merge and export eight plots
sem_md.a_all_plot_all = ggarrange(sem_estab_md.a_all_space1, sem_estab_md.a_all_space2,
                                  sem_estab_effects_md.a_all.dir_indir.mpd_all.plot,
                                  sem_estab_effects_md.a_all.percent.plot,
                                  sem_domin_md.a_all_space1,sem_domin_md.a_all_space2,
                                  sem_domin_effects_md.a_all.dir_indir.mpd_all.plot,
                                  sem_domin_effects_md.a_all.percent.plot,
                                  nrow = 4, ncol = 2,
                                  heights = c(0.35, 0.15, 0.35, 0.15),
                                  labels = c('(a)', '', '','',
                                             '(b)', '', '',''),
                                  vjust = 1.8)

emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/sem_md.a_all_plot_all.emf',
    width=16, height=24, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
sem_md.a_all_plot_all
dev.off() #turn off device and finalize file


#### Draw plot for mnd_all  
### Establishment corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_estab_sem1_mnd_all = sqrt(predict(estab_sem1_mnd_all[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(estab_sem1_mnd_all[[1]])))+pi^2/3)
sd.y_estab_sem1_mnd_all
std.estab_sem1_mnd_all=unlist(fixef(estab_sem1_mnd_all[[1]]))[-1]*
  c(sd(dat_suc_sps$mntd_all),
    sd(dat_suc_sps$mnfunc_d_all),
    sd(dat_suc_sps$mnnd),
    sd(dat_suc_sps$mnlgfd))/sd.y_estab_sem1_mnd_all
std.estab_sem1_mnd_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_estab_sem1_mnd_all=coefs(estab_sem1_mnd_all)# original coefficients of piecewiseSEM
coefs_estab_sem1_mnd_all[1:4,8]=std.estab_sem1_mnd_all

### Dominance corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_domin_sem1_mnd_all=sqrt(predict(domin_sem1_mnd_all[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(domin_sem1_mnd_all[[1]])))+pi^2/3)
sd.y_domin_sem1_mnd_all
std.domin_sem1_mnd_all=unlist(fixef(domin_sem1_mnd_all[[1]]))[-1]*
  c(sd(dat_dom_sps$mntd_all),
    sd(dat_dom_sps$mnfunc_d_all),
    sd(dat_dom_sps$mnnd),
    sd(dat_dom_sps$mnlgfd))/sd.y_domin_sem1_mnd_all
std.domin_sem1_mnd_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_domin_sem1_mnd_all=coefs(domin_sem1_mnd_all)# original coefficients of piecewiseSEM
coefs_domin_sem1_mnd_all[1:4,8]=std.domin_sem1_mnd_all
coefs_domin_sem1_mnd_all

### Calculate direct and indirect effects of SEM following Xu Meng (2022) GCB
### Estab
sem_estab_effects_mnd_all = data.frame(predictor = rep(unique(coefs_estab_sem1_mnd_all$Predictor)[1:4], 3),
                                       type = rep(c('direct', 'indirect', 'total'), each = 4),
                                       value = rep(0, 12),
                                       significant = rep('yes', 12))
coefs_estab_sem1_mnd_all = coefs_estab_sem1_mnd_all[-nrow(coefs_estab_sem1_mnd_all), ]
coefs_estab_sem1_mnd_all_c = coefs_estab_sem1_mnd_all[coefs_estab_sem1_mnd_all$P.Value < 0.05,]
for (i in 1:length(unique(coefs_estab_sem1_mnd_all_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_estab_sem1_mnd_all_c$Predictor)[i]
  dat = coefs_estab_sem1_mnd_all_c[coefs_estab_sem1_mnd_all_c$Predictor == predictor_1,]
  dat_2 = coefs_estab_sem1_mnd_all_c[coefs_estab_sem1_mnd_all_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_estab_effects_mnd_all[sem_estab_effects_mnd_all$predictor == predictor_1&
                                sem_estab_effects_mnd_all$type == 'direct',]$value = dat_3$Std.Estimate
    sem_estab_effects_mnd_all[sem_estab_effects_mnd_all$predictor == predictor_1&
                                sem_estab_effects_mnd_all$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    if (nrow(dat_3[dat_3$Response == 'estab' &
                   dat_3$Predictor == predictor_1,]) > 0){
      sem_estab_effects_mnd_all[sem_estab_effects_mnd_all$predictor == predictor_1&
                                  sem_estab_effects_mnd_all$type == 'direct',]$value = dat_3[dat_3$Response == 'estab' &
                                                                                               dat_3$Predictor == predictor_1,]$Std.Estimate
    }
    sem_estab_effects_mnd_all[sem_estab_effects_mnd_all$predictor == predictor_1&
                                sem_estab_effects_mnd_all$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'estab',]$Std.Estimate*dat_3[dat_3$Response == 'estab' &
                                                                                                                                                  dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_estab_effects_mnd_all[sem_estab_effects_mnd_all$predictor == predictor_1&
                                sem_estab_effects_mnd_all$type == 'total',]$value  = sem_estab_effects_mnd_all[sem_estab_effects_mnd_all$predictor == predictor_1&
                                                                                                                 sem_estab_effects_mnd_all$type == 'direct',]$value + sem_estab_effects_mnd_all[sem_estab_effects_mnd_all$predictor == predictor_1&
                                                                                                                                                                                                  sem_estab_effects_mnd_all$type == 'indirect',]$value
  }
}

### Domin
sem_domin_effects_mnd_all = data.frame(predictor = rep(unique(coefs_domin_sem1_mnd_all$Predictor)[1:4], 3),
                                       type = rep(c('direct', 'indirect', 'total'), each = 4),
                                       value = rep(0, 12))
coefs_domin_sem1_mnd_all = coefs_domin_sem1_mnd_all[-nrow(coefs_domin_sem1_mnd_all),]
coefs_domin_sem1_mnd_all_c = coefs_domin_sem1_mnd_all[coefs_domin_sem1_mnd_all$P.Value < 0.05,]
for (i in 1:length(unique(coefs_domin_sem1_mnd_all_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_domin_sem1_mnd_all_c$Predictor)[i]
  dat = coefs_domin_sem1_mnd_all_c[coefs_domin_sem1_mnd_all_c$Predictor == predictor_1,]
  dat_2 = coefs_domin_sem1_mnd_all_c[coefs_domin_sem1_mnd_all_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_domin_effects_mnd_all[sem_domin_effects_mnd_all$predictor == predictor_1&
                                sem_domin_effects_mnd_all$type == 'direct',]$value = dat_3$Std.Estimate
    sem_domin_effects_mnd_all[sem_domin_effects_mnd_all$predictor == predictor_1&
                                sem_domin_effects_mnd_all$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    if (nrow(dat_3[dat_3$Response == 'domin' &
                   dat_3$Predictor == predictor_1,]) > 0) {
      sem_domin_effects_mnd_all[sem_domin_effects_mnd_all$predictor == predictor_1&
                                  sem_domin_effects_mnd_all$type == 'direct',]$value = dat_3[dat_3$Response == 'domin' &
                                                                                               dat_3$Predictor == predictor_1,]$Std.Estimate
    }
    sem_domin_effects_mnd_all[sem_domin_effects_mnd_all$predictor == predictor_1&
                                sem_domin_effects_mnd_all$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'domin',]$Std.Estimate*dat_3[dat_3$Response == 'domin' &
                                                                                                                                                  dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_domin_effects_mnd_all[sem_domin_effects_mnd_all$predictor == predictor_1&
                                sem_domin_effects_mnd_all$type == 'total',]$value  = sem_domin_effects_mnd_all[sem_domin_effects_mnd_all$predictor == predictor_1&
                                                                                                                 sem_domin_effects_mnd_all$type == 'direct',]$value + sem_domin_effects_mnd_all[sem_domin_effects_mnd_all$predictor == predictor_1&
                                                                                                                                                                                                  sem_domin_effects_mnd_all$type == 'indirect',]$value
  }
}

#-----------Draw plots for SEMs' relative total effects, direct and indirect effects 
### Estab
require(ggthemes)
sem_estab_effects_mnd_all_total = sem_estab_effects_mnd_all %>% filter(type == 'total')
sem_estab_effects_mnd_all.percent = sem_estab_effects_mnd_all_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_estab_effects_mnd_all.longer = sem_estab_effects_mnd_all %>% filter(type != 'total' 
                                                                        #& value != 0
)

# Relative total effect 
sem_estab_effects_mnd_all.percent.plot = 
  ggplot(data=sem_estab_effects_mnd_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mnfunc_d_all", "mntd_all", "mnlgfd", "mnnd"),
                   labels=c("I-N MNFD_all", "I-N MNTD_all", "I-N MNRFD", "I-N MNND"),
                   position="right")+
  scale_fill_stata()+# color palatte in ggthemes
  xlim(0,0.8)+
  xlab("")+
  #theme_test()+
  theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))+
  guides(fill="none")

# Direct and indirect effects
sem_estab_effects_mnd_all.dir_indir.mpd_all.plot = 
  ggplot(data=sem_estab_effects_mnd_all.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects","indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mnfunc_d_all", "mntd_all", "mnlgfd", "mnnd"),
                   position="right")+
  scale_x_continuous(breaks = c(seq(-0.4, 0.1, 0.2)))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_estab_mnd_all_space = ggplot()+
  facet_wrap(~"Establishment")+
  theme_test() + 
  theme(strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.015,1.2,0.2),units="lines"))

# Merge three plots
sem_estab_mnd_all_plot_all = ggarrange(sem_estab_mnd_all_space,
                                       sem_estab_effects_mnd_all.dir_indir.mpd_all.plot,
                                       sem_estab_effects_mnd_all.percent.plot,
                                       nrow = 1, ncol = 3,
                                       labels = c('(a)', '', ''),
                                       vjust = 1.8)


ggsave(plot = sem_estab_mnd_all_plot_all,
       "results/figures_ages1_35_top50_equal_interval_bh_partialb/all_species_d/sem_estab_mnd_all_plot_all.svg",
       width=9,height=3)

### domin
sem_domin_effects_mnd_all_total = sem_domin_effects_mnd_all %>% filter(type == 'total')
sem_domin_effects_mnd_all.percent = sem_domin_effects_mnd_all_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_domin_effects_mnd_all.longer = sem_domin_effects_mnd_all %>% filter(type != 'total' 
                                                                        #&value != 0
)

# Relative total effect
sem_domin_effects_mnd_all.percent.plot = 
  ggplot(data=sem_domin_effects_mnd_all.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mnfunc_d_all", "mntd_all", "mnlgfd", "mnnd"),
                   labels=c("E-N MNFD_all", "E-N MNTD_all", "E-N MNRFD", "E-N MNND"),
                   position="right")+
  scale_fill_stata()+# color palatte in ggthemes
  xlim(0,0.8)+
  xlab("")+
  #theme_test()+
  theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))+
  guides(fill="none")

# Direct and indirect effects
sem_domin_effects_mnd_all.dir_indir.mpd_all.plot = 
  ggplot(data=sem_domin_effects_mnd_all.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects",
                                          "indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mnfunc_d_all", "mntd_all", "mnlgfd", "mnnd"),
                   position="right")+
  scale_x_continuous(breaks = seq(-0.1, 0.1, 0.1),
                     limits = c(-0.2, 0.2))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_domin_mnd_all_space = ggplot()+
  facet_wrap(~"Dominance")+
  theme_test() + 
  theme(strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.015,1.2,0.2),units="lines"))

# Merge three plots
sem_domin_mnd_all_plot_all = ggarrange(sem_domin_mnd_all_space,
                                       sem_domin_effects_mnd_all.dir_indir.mpd_all.plot,
                                       sem_domin_effects_mnd_all.percent.plot,
                                       nrow = 1, ncol = 3,
                                       labels = c('(b)', '', ''),
                                       vjust = 1.8)


ggsave(plot = sem_domin_mnd_all_plot_all,
       "results/figures_ages1_35_top50_equal_interval_bh_partialb/all_species_d/sem_domin_mnd_all_plot_all.svg",
       width=9,height=3)

