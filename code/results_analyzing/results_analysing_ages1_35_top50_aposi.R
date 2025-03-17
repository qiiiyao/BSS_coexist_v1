# Load packages
rm(list = ls())

require(stringr)
require(dplyr)
require(data.table)

load("D:/R projects/BSS/code/data preparation/transformed data/fit_fp_top50_ages1_35.RData")
load("D:/R projects/BSS/code/data preparation/transformed data/re_cover_ab_f1_ages1_35_fp.rdata")
dat_all_l = lapply(1:length(fit_fp_top50_ages1_35), function(x){
  #x = 1
  stage_sp = fit_fp_top50_ages1_35[[x]]$stage_sp_name
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
     file = 'results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_aposi/dat_all_alltime.rdata')
rm(list = ls()) ### free the memory

############### Fast Start ####################
rm(list = ls())
setwd("D:/R projects/BSS")
load('results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_aposi/inter_all_c_alltime.rdata')
load('results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50/dat_all_alltime.rdata')

library(dplyr)
library(betareg)
library(stringr)
library(data.table)
library(ape)
library(reticulate)
#inter_all_c = inter_all_c_alltime
#rm(inter_all_c_alltime)


gre = "#55a868ff"
ora = "#dd8452ff"
yel = "#ccb974ff"
blu = "#4c72b0ff"
colors_4d = c(blu,ora,gre,yel)
colors_2d = c(ora, blu)

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

# get the abs value of fd
inter_all_c$abfd = abs(inter_all_c$fd)

# calculate coexistence
do.coexistence = function(stabilizing, equalizing) {
  ni = stabilizing
  fi = equalizing
  if (ni > fi & ni > 0) {
    outcome = data.frame(coexists = 1, priority = 0, i_exclusion_j = 0,
                         j_exclusion_i = 0)
  } else if (ni < fi & ni < 0) {
    outcome = data.frame(coexists = 0, priority = 1, i_exclusion_j = 0,
                         j_exclusion_i = 0)
  } else if (ni > fi & ni < 0) {
    outcome = data.frame(coexists = 0, priority = 0, i_exclusion_j = 1,
                         j_exclusion_i = 0)
  } else {
    outcome = data.frame(coexists = 0, priority = 0, i_exclusion_j = 0,
                         j_exclusion_i = 1)
  }
  return(outcome)
}

#outcomes = rbindlist((apply(inter_all_c, 1, function(x){
 # outcomes = do.coexistence(stabilizing = as.numeric(x['nd']),
 #                          equalizing = as.numeric(x['fd']))     
#})))

#inter_all_c = cbind(inter_all_c, outcomes)

#save(inter_all_c,
  #  file = 'results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_aposi/inter_all_c_alltime.rdata')

#### Check the correlation between nd and fd
#cor.test(inter_all_c$nd, inter_all_c$fd) 
## Full data but corr coeff just 0.13
## correlation
#plot(inter_all_c$nd, inter_all_c$abfd)

## Add stage introduce
inter_all_c = as.data.frame(inter_all_c)
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

### Check the proportion of positive interactions
inter_all_c = inter_all_c %>% mutate(positive_inter = (aij < 0))
inter_all_c = inter_all_c %>% mutate(positive_intra = (aii < 0))
hist(as.numeric(inter_all_c$positive_inter))
nrow(inter_all_c %>% filter(positive_inter == T))/nrow(inter_all_c)
# positive_inter% = 47.75%
hist(as.numeric(inter_all_c$positive_intra))
nrow(inter_all_c %>% filter(positive_intra == T))/nrow(inter_all_c)
# positive_intra% = 12.61%


#### Draw coexsit_probs_pie plots ####
inter_all_c = as.data.frame(inter_all_c)
inter_all_forexotics = inter_all_c %>%
  filter(stage_i != 'native', stage_j == 'native')
length(unique(inter_all_forexotics$species_i))
length(unique(inter_all_forexotics$species_j))

inter_all_forestab = inter_all_c %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))

inter_all_forestab_msd = inter_all_forestab %>% 
  group_by(stage_ij_estab) %>% summarise_at(vars(nd, fd),
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

BoundaryValues = data.frame(x = seq(0,4.99,l=100), 
                            y1 = seq(0,4.99,l=100),
                            y2 = -5)      

BoundaryValues_P = data.frame(x = seq(-5,0,l=100), 
                              y1 = seq(-5,0,l=100),
                              y2 = 5)      
OneLine = data.frame(x= c(rev(BoundaryValues$x),BoundaryValues$x),
                     y= c(rev(BoundaryValues$y2),BoundaryValues$y1 ))


#"#4c72b0ff" "#dd8452ff" "#55a868ff" "#ccb974ff"
Fig.S1_compare_estab_original_pie = ggplot()+
  geom_rect(aes( xmin = 5, xmax = 0, ymax = 5, ymin =-5), fill = '#55a868ff',
            alpha = 0.4)+
  geom_rect(aes( xmin = 0, xmax =-5, ymax = 5, ymin =-5), fill = '#ccb974ff',
            alpha = 0.4)+
  geom_ribbon(data = BoundaryValues,
              aes(x = x, ymin = y2, ymax = y1), inherit.aes = FALSE, fill = 'grey75')+
  geom_ribbon(data = BoundaryValues_P,
              aes(x = x, ymin = y1, ymax = y2), inherit.aes = FALSE, fill = 'grey90')+
  geom_path(data = BoundaryValues,
            aes(x, y1), col = 'red', size = 1.5)+
  geom_line(data = BoundaryValues_P,
            aes(x, y1), col = 'red', size = 1.5)+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_vline(xintercept = 0)+
  geom_scatterpie(aes(x = nd, y = fd, group = stage_ij_estab), 
                  data = inter_all_forestab_outcomes,
                  pie_scale = 1500,
                  cols= c('j_exclusion_i_probs', 'i_exclusion_j_probs',
                          'coexists_probs' ,'priority_probs'))+
  coord_fixed(xlim =c(-2, 3), ylim= c(-3.5, 1.5))+
  geom_point(aes(x = nd, y = fd, col = stage_ij_estab), size = 6,
             inter_all_forestab_outcomes)+ 
  geom_path(aes(nd, fd), arrow = arrow(length = unit(0.3, "cm"), ends = 'first'), size =1.5,
            data = inter_all_forestab_outcomes)+
  labs(x = 'Niche difference (ND)', 
       y = 'Relative fitness difference\n(RFD)') +
  annotate(geom="text", x=-1.2, y=1.4, family = 'Arial',
           label='Establishment', fontface = 'bold')+
  scale_colour_manual(values = c(colors_4d[2],colors_4d[1]),
                      labels = c("intro.estab_native" = "Successful",
                                 "intro.noestab_native" = "Failed"),
                      name = ' ')+ 
  scale_fill_manual(values = c('#55a868ff', '#ccb974ff', 'grey55','grey85'),
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
  theme(legend.position = c(0.2, 0.5),
        legend.direction = "vertical", legend.box = "vertical",
        legend.box.background = element_rect(colour = 'black',
                                             fill = "white",
                                          linewidth = 0.5),
        legend.title = element_text(margin = margin(t = 0, r = 0, b = 5, l = 0)),  # create more white space above the legend titles
        legend.box.margin = margin(t = 5, r = -10, b = 5, l = 5)
  ) 
Fig.S1_compare_estab_original_pie

# Dominance
inter_all_fordomin = inter_all_c %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))
inter_all_fordomin_msd = inter_all_fordomin %>% 
  group_by(stage_ij_domin) %>% summarise_at(vars(nd, fd),
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

Fig.S1_compare_domin_original_pie =  ggplot()+
  geom_rect(aes( xmin = 5, xmax = 0, ymax = 5, ymin =-5), fill = '#55a868ff',
            alpha = 0.4)+
  geom_rect(aes( xmin = 0, xmax =-5, ymax = 5, ymin =-5), fill = '#ccb974ff',
            alpha = 0.4)+
  geom_ribbon(data = BoundaryValues,
              aes(x = x, ymin = y2, ymax = y1), inherit.aes = FALSE, fill = 'grey75')+
  geom_ribbon(data = BoundaryValues_P,
              aes(x = x, ymin = y1, ymax = y2), inherit.aes = FALSE, fill = 'grey90')+
  geom_path(data = BoundaryValues,
            aes(x, y1), col = 'red', size = 1.5)+
  geom_line(data = BoundaryValues_P,
            aes(x, y1), col = 'red', size = 1.5)+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_vline(xintercept = 0)+
  geom_scatterpie(aes(x = nd, y = fd, group = stage_ij_domin), 
                  data = inter_all_fordomin_outcomes,
                  pie_scale = 600,
                  cols= c('j_exclusion_i_probs', 'i_exclusion_j_probs',
                          'coexists_probs' ,'priority_probs'))+
  coord_fixed(xlim =c(-2, 3), ylim= c(-4.5, 0.5))+
  geom_point(aes(x = nd, y = fd, col = stage_ij_domin), size = 6,
             inter_all_fordomin_outcomes)+ 
  geom_path(aes(nd, fd), arrow = arrow(length = unit(0.3, "cm"), ends = 'first'), size =1.5,
            data = inter_all_fordomin_outcomes)+
  labs(x = 'Niche difference (ND)', 
       y = '  ') +
  annotate(geom="text", x=-1.25, y=0.4, family = 'Arial',
           label='Dominance', fontface = 'bold')+
  scale_colour_manual(values = c(colors_4d[2],colors_4d[1]),
                      labels = c("estab.domin_native" = "Successful",
                                 "estab.nodomin_native" = "Failed"),
                      name = ' ')+ 
  scale_fill_manual(values = c( '#55a868ff', '#ccb974ff', 'grey55','grey85'),
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
  theme(#legend.position = c(0.4, 0.7),
        legend.position = 'none',
        legend.direction = "vertical", legend.box = "horizontal",
        legend.box.background = element_rect(colour = 'black',
                                             fill = "white",
                                             linewidth = 0.1)#,
        #legend.title = element_text(margin = margin(t = 140)),  # create more white space above the legend titles
        #legend.box.margin = margin(b = -80)
  ) 
Fig.S1_compare_domin_original_pie


#### Stage split two parts: estab noestab domin nodomin ####
##### just mean and sd of raw data #####
inter_all_forexotics = inter_all_c %>%
  filter(stage_i != 'native', stage_j == 'native')
length(unique(inter_all_forexotics$species_i))
length(unique(inter_all_forexotics$species_j))

inter_all_forestab = inter_all_c %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))

spss.f = function(x){y=x
                     return(y)}
x = seq(-1, 1, 0.001)

length(unique(inter_all_forestab$sp_pair))
length(unique(inter_all_forestab$species_i))
length(unique(inter_all_forestab$species_j))

inter_all_forestab_msd = inter_all_forestab %>% 
  group_by(stage_ij_estab) %>% summarise_at(vars(nd, fd),
                                            c(mean, sd), na.rm = T)

colnames(inter_all_forestab_msd) = c('stage', 'nd_mean', 'fd_mean',
                                     'nd_sd',
                                     'fd_sd')

Fig.1_compare_estab_original = ggplot(NULL) +
  geom_point(data = inter_all_forestab, aes(x = nd, y = fd,
                                            color = stage_ij_estab,
                                            shape = stage_ij_estab,
                                            alpha = stage_ij_estab),
             size = 1)+
  geom_pointrange(data = inter_all_forestab_msd,
                  mapping = aes(x = nd_mean,y = fd_mean,
                                ymin = fd_mean-fd_sd,
                                ymax = fd_mean+fd_sd,
                                color = stage,
                                shape = stage),
                  alpha = 1.2
                  ,fatten = 8,linewidth = 0.8
  )+
  geom_pointrange(data = inter_all_forestab_msd,
                  mapping = aes(x = nd_mean,y = fd_mean,
                                xmin = nd_mean-nd_sd,
                                xmax = nd_mean+nd_sd,
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
  #scale_color_brewer(palette="Set2")+
  scale_x_continuous(limits=c(-10, 10)) +
  scale_y_continuous(limits=c(-12, 12)) +
  #geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.3) +
  annotate(geom="text", x=-5, y=10, family = 'Arial',
           label='Establishment')+
  labs(x = 'Niche difference (ND)', 
       y = 'Relative fitness difference\n(RFD)') +
  stat_function(fun=spss.f, colour="black", linewidth = 1.5, linetype = 'dashed')+
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(linewidth = 0.3),
                             family = 'Arial'),
         shape=guide_legend(title=NULL,
                            family = 'Arial')) +
  theme_for_coe_plot() + 
  theme(legend.position = c(0.18, 0.12))
Fig.1_compare_estab_original


inter_all_fordomin = inter_all_c %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))

inter_all_fordomin_msd = inter_all_fordomin %>% 
  group_by(stage_ij_domin) %>% summarise_at(vars(nd, fd),
                                            c(mean, sd), na.rm = T)

colnames(inter_all_fordomin_msd) = c('stage', 'nd_mean',
                                     'fd_mean', 'nd_sd',
                                     'fd_sd')

Fig.1_compare_domin_original = ggplot(NULL) +
  geom_point(data = inter_all_fordomin, aes(x = nd, y = fd,
                                            color = stage_ij_domin,
                                          shape = stage_ij_domin),
            alpha = 0.05, size = 1)+
  geom_pointrange(data = inter_all_fordomin_msd,
                  mapping = aes(x = nd_mean,y = fd_mean,
                                ymin = fd_mean-fd_sd,
                                ymax = fd_mean+fd_sd,
                                color = stage,
                                shape = stage)
                  ,fatten = 8,linewidth = 0.8
  )+
  geom_pointrange(data = inter_all_fordomin_msd,
                  mapping = aes(x = nd_mean, y = fd_mean,
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
  scale_x_continuous(limits=c(-10, 10)) +
  scale_y_continuous(limits=c(-12.5, 12.5)) +
  #geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.3) +
  annotate(geom="text", x=-5, y=10, family = 'Arial',
           label='Dominance')+
  labs(x = 'Niche difference (ND)', 
       y = '  ') +
  stat_function(fun=spss.f, colour="black", linewidth = 1.5, linetype = 'dashed')+
  guides(colour=guide_legend(title=NULL,
                             override.aes = list(linewidth = 0.3)),
         shape=guide_legend(title=NULL,
                            override.aes = list(linewidth = 0.3))) +
  theme_for_coe_plot() + 
  theme(legend.position = c(0.18, 0.12)) 
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

mod_fd_compare_estab_1 = lmer(fd ~ stage_ij_estab+(1|field/plot),
                              data = inter_all_forestab,
                              REML = T)
mod_fd_compare_estab = lmer(fd ~ stage_ij_estab+(1|field/plot)+
                              (1|species_i) + (1|species_j) +
                              (1|sp_pair),
                            data = inter_all_forestab,
                            REML = T)
mod_fd_compare_estab_2 = lmer(fd ~ stage_ij_estab+(1|field/plot)+
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
                                    fd = pre_mod_fd_compare_estab$predicted,
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
  geom_point(data = inter_all_forestab, aes(x = nd, y = fd,
                                            color = stage_ij_estab),
             alpha = 0.1, size = 30)+
  geom_point(data = pre_ndfd_compare_estab,
             mapping = aes(x = nd,y = fd),
             color = 'black',
             size = 50)+
  geom_pointrange(data = pre_ndfd_compare_estab,
                  mapping = aes(x = nd,y = fd, ymin = fd_low,
                                ymax = fd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  geom_pointrange(data = pre_ndfd_compare_estab,
                  mapping = aes(y = fd,x = nd, xmin = nd_low,
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

mod_fd_compare_domin_1 = lmer(fd ~ stage_ij_domin+(1|field/plot),
                              data = inter_all_fordomin,
                              REML = T)
mod_fd_compare_domin_2 = lmer(fd ~ stage_ij_domin+(1|field/plot)+(1|sp_pair),
                              data = inter_all_fordomin,
                              REML = T)
mod_fd_compare_domin = lmer(fd ~ stage_ij_domin+(1|field/plot)+
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
                                    fd = pre_mod_fd_compare_domin$predicted,
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
  geom_point(data = inter_all_fordomin, aes(x = nd, y = fd,
                                            color = stage_ij_domin
  ), alpha = 0.1,
  size = 30)+
  geom_point(data = pre_ndfd_compare_domin,
             mapping = aes(x = nd,y = fd),
             color = 'black',
             size = 50)+
  geom_pointrange(data = pre_ndfd_compare_domin,
                  mapping = aes(x = nd,y = fd, ymin = fd_low,
                                ymax = fd_high,
                                color = stage),
                  fatten = 70,linewidth = 3)+
  geom_pointrange(data = pre_ndfd_compare_domin,
                  mapping = aes(y = fd,x = nd, xmin = nd_low,
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

#ggsave('results/figures_ages1_35_top50_aposi/Fig.1_compare_domin.svg',
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
       'results/figures_ages1_35_top50_aposi/Fig.1_split_stage.svg',
       width = 320,height = 100, dpi = 300, units = 'cm',
       limitsize = F)

Fig.1_split_stage_2 = ggarrange(gap,Fig.1_compare_estab,gap,
                                Fig.1_compare_domin, nrow = 4, ncol = 1,
                                labels = c('','a)','','b)'), hjust = -0.5,
                                vjust = 0, heights = c(0.1, 1, 0.1, 1),
                                font.label = list(size = 150))

ggsave(plot = Fig.1_split_stage_2,
       'results/figures_ages1_35_top50_aposi/Fig.1_split_stage_2.svg',
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
    #j = 1
    inv_suc_sp = inv_suc %>% filter(species_i == inv_sp[j])
    inv_suc_sp_all = inv_suc_all %>% filter(species_i == inv_sp[j])
    mnd = mean(inv_suc_sp$nd)
    mfd = mean(inv_suc_sp$fd)
    mabfd = mean(inv_suc_sp$abfd)
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
    
    mnd.a = sum(inv_suc_sp$nd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mfd.a = sum(inv_suc_sp$fd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mabfd.a = sum(inv_suc_sp$abfd*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mpd.a = sum(inv_suc_sp$Phylo_dis*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mfunc_d.a = sum(inv_suc_sp$Multi_traits*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mconti_func_d.a = sum(inv_suc_sp$Multi_conti_traits*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mgrowth.a = sum(inv_suc_sp$growth*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mspan.a = sum(inv_suc_sp$span*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mpollination.a = sum(inv_suc_sp$pollination*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mdispersal.a = sum(inv_suc_sp$dispersal*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mclonality.a = sum(inv_suc_sp$clonality*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mheight.a = sum(inv_suc_sp$height*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mldmc.a = sum(inv_suc_sp$ldmc*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    msla.a = sum(inv_suc_sp$sla*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    mseedmass.a = sum(inv_suc_sp$seedmass*inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)/sum(inv_suc_sp$ra_m_real_t_i*inv_suc_sp$ra_m_real_t_j)
    
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
    
    mnnd = min(inv_suc_sp$nd)
    mnfd = min(inv_suc_sp$fd)
    mnabfd = min(inv_suc_sp$abfd)
    mntd = min(inv_suc_sp$Phylo_dis)
    mnfunc_d = min(inv_suc_sp$Multi_traits)
    mnconti_func_d = min(inv_suc_sp$Multi_conti_traits)
    mngrowth = min(inv_suc_sp$growth)
    mnspan = min(inv_suc_sp$span)
    mnpollination = min(inv_suc_sp$pollination)
    mndispersal = min(inv_suc_sp$dispersal)
    mnclonality = min(inv_suc_sp$clonality)
    mnheight = min(inv_suc_sp$height)
    mnldmc = min(inv_suc_sp$ldmc)
    mnsla = min(inv_suc_sp$sla)
    mnseedmass = min(inv_suc_sp$seedmass)
    
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
    
    
    dat_suc_sp_1 = data.frame(f_p = unique(trans_plot$f_p), plot = unique(trans_plot$f_p), field = unique(trans_plot$field),
                              species = inv_sp[j], 
                              stage = unique(inv_suc_sp$stage_i),
                              estab = unique(inv_suc_sp$estab_i),
                              domin = unique(inv_suc_sp$domin_i),
                              mnd = mnd, mfd = mfd, mabfd = mabfd,
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
                              
                              mnd.a = mnd.a, mfd.a = mfd.a,
                              mabfd.a = mabfd.a, mpd.a = mpd.a,
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
                              
                              mnnd = mnnd, mnfd = mnfd,
                              mnabfd = mnabfd, mntd = mntd,
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/dat_suc_sp.rdata')



##### Community level #####
######## top25 exotics from ages 3-51 ########
tree = read.tree('data/original data/phylo_tree332.txt')
tree_dis = cophenetic(tree)
trait = read.csv('data/original data/traits332.csv',
                 header = T)
source('code/invasion impact/measure_impact_func.R')
source('code/data preparation/code/function/data collation_function.R')
source('others_example/Li_et_al_2015_EL/function/Li_suc_prepare_ages_4_distance.R')
source('others_example/Li_et_al_2015_EL/function/mpd_calculate.R')
load("code/data preparation/transformed data/sp_cover_ages1_35_top25_exotic_fp.rdata")
load("code/Phylo_Func/fd_gower_mat.rdata")

function_dis = fd_gower_mat
inter_all_c_alltime_l = split(inter_all_c_alltime, inter_all_c_alltime$f_p)
nd_rfd_matrix = lapply(inter_all_c_alltime_l, function(x){
  niche_dat = tidyr::spread(x[,c('species_i', 'species_j', 'nd')], species_i, nd, fill = 0)
  rownames(niche_dat) = niche_dat[,1]
  niche_dis = as.matrix(niche_dat[,-1])
  fitness_dat = tidyr::spread(x[,c('species_i', 'species_j', 'fd')], species_i, fd, fill = 0)
  rownames(fitness_dat) = fitness_dat[,1]
  fitness_dis = as.matrix(fitness_dat[,-1])
  fit_natives = unique((x %>% filter(stage_i == 'native'))$species_i)
  fit_exotics = unique((x %>% filter(stage_i != 'native'))$species_i)
  results = list(niche_dis = niche_dis, fitness_dis = fitness_dis, 
                 fit_natives = fit_natives, fit_exotics = fit_exotics)
  return(results)
})

sp_cover_ages1_35_top25_exotic_fp_all = sapply(c(1:length(sp_cover_ages1_35_top25_exotic_fp)),
                                               function(x){f_p = names(sp_cover_ages1_35_top25_exotic_fp)[x]
                                               list(data_pd = sp_cover_ages1_35_top25_exotic_fp[[x]],
                                                    data_coexist = nd_rfd_matrix[[which(names(nd_rfd_matrix) == f_p)]])},
                                               simplify=FALSE)
names(sp_cover_ages1_35_top25_exotic_fp_all) = names(sp_cover_ages1_35_top25_exotic_fp)

dat_suc_com = lapply(sp_cover_ages1_35_top25_exotic_fp_all,
                     Li_suc_prepare_ages_4_distance)
save(dat_suc_com,
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/dat_suc_com.rdata')
require(lme4)
require(dplyr)


##### Single regression #####
###### estab ######
load('code/results_analyzing/analysing_ages1_35_top50_aposi_data/dat_suc_com.rdata')
dat_estab_com = lapply(dat_suc_com, function(x){
  estab = x$estab
  return(estab)
})

dat_estab_com = rbindlist(dat_estab_com)
dat_estab_com = as.data.frame(dat_estab_com %>% filter(!is.infinite(mnnd)))
dat_estab_com$plot = sapply(stringr::str_split(dat_estab_com$f_p, '_'),
                            function(x){x[2]})

dat_estab_com_ages15_35 = dat_estab_com %>% filter(fake_age %in% c(15:35))
####### Phylogenetic differences
estab.mpd_single = glmer(estab~mpd+(1|field/f_p)+(1|fake_age) ,
                         family=binomial(link = "logit"),
                         data=dat_estab_com)
summary(estab.mpd_single)

estab.mpd.a_single = glmer(estab~mpd.a+(1|field/f_p)+(1|fake_age) ,
                           family=binomial(link = "logit"),
                           data=dat_estab_com)
summary(estab.mpd.a_single)

estab.mntd_single = glmer(estab~mntd+(1|field/f_p)+(1|fake_age) ,
                          family=binomial(link = "logit"),
                          data=dat_estab_com)
summary(estab.mntd_single)

####### Functional differences
estab.mfd_single = glmer(estab~mfd+(1|field/f_p)+(1|fake_age) ,
                         family=binomial(link = "logit"),
                         data=dat_estab_com)
summary(estab.mfd_single)

estab.mfd.a_single = glmer(estab~mfd.a+(1|field/f_p)+(1|fake_age) ,
                           family=binomial(link = "logit"),
                           data=dat_estab_com)
summary(estab.mfd.a_single)

estab.mnfd_single = glmer(estab~mnfd+(1|field/f_p)+(1|fake_age) ,
                          family=binomial(link = "logit"),
                          data=dat_estab_com)
summary(estab.mnfd_single)

####### Niche differences
estab.mnd_single = glmer(estab~mnd+(1|field/f_p)+(1|fake_age) ,
                         family=binomial(link = "logit"),
                         data=dat_estab_com)
summary(estab.mnd_single)

estab.mnd.a_single = glmer(estab~mnd.a+(1|field/f_p)+(1|fake_age) ,
                           family=binomial(link = "logit"),
                           data=dat_estab_com)
summary(estab.mnd.a_single)

estab.mnnd_single = glmer(estab~mnnd+(1|field/f_p)+(1|fake_age) ,
                          family=binomial(link = "logit"),
                          data=dat_estab_com)
summary(estab.mnnd_single)

####### Relative fitness differences
estab.mrfd_single = glmer(estab~mrfd+(1|field/f_p)+(1|fake_age) ,
                          family=binomial(link = "logit"),
                          data=dat_estab_com)
summary(estab.mrfd_single)

estab.mrfd.a_single = glmer(estab~mrfd.a+(1|field/f_p)+(1|fake_age) ,
                            family=binomial(link = "logit"),
                            data=dat_estab_com)
summary(estab.mrfd.a_single)

estab.mnrfd_single = glmer(estab~mnrfd+(1|field/f_p)+(1|fake_age) ,
                           family=binomial(link = "logit"),
                           data=dat_estab_com)
summary(estab.mnrfd_single)

####### domin #######
dat_domin_com = lapply(dat_suc_ages1_35_top25_exotic_fp_all, function(x){
  domin = x$domin
  return(domin)
})

dat_domin_com = rbindlist(dat_domin_com)
dat_domin_com = as.data.frame(dat_domin_com %>% filter(!is.infinite(mnnd)))
dat_domin_com$plot = as.integer(dat_domin_com$plot)
dat_domin_com$field = as.integer(dat_domin_com$field)
dat_domin_com$fake_age = as.integer(dat_domin_com$fake_age)

dat_domin_com$plot = sapply(stringr::str_split(dat_domin_com$f_p, '_'),
                            function(x){x[2]})

####### Phylogenetic differences
domin.mpd_single = glmer(domin~mpd+(1|field/f_p)+(1|fake_age) ,
                         family=binomial(link = "logit"),
                         data=dat_domin_com)
summary(domin.mpd_single)

domin.mpd.a_single = glmer(domin~mpd.a+(1|field/f_p)+(1|fake_age) ,
                           family=binomial(link = "logit"),
                           data=dat_domin_com)
summary(domin.mpd.a_single)

domin.mntd_single = glmer(domin~mntd+(1|field/f_p)+(1|fake_age) ,
                          family=binomial(link = "logit"),
                          data=dat_domin_com)
summary(domin.mntd_single)

####### Functional differences
domin.mfd_single = glmer(domin~mfd+(1|field/f_p)+(1|fake_age) ,
                         family=binomial(link = "logit"),
                         data=dat_domin_com)
summary(domin.mfd_single)

domin.mfd.a_single = glmer(domin~mfd.a+(1|field/f_p)+(1|fake_age) ,
                           family=binomial(link = "logit"),
                           data=dat_domin_com)
summary(domin.mfd.a_single)

domin.mnfd_single = glmer(domin~mnfd+(1|field/f_p)+(1|fake_age) ,
                          family=binomial(link = "logit"),
                          data=dat_domin_com)
summary(domin.mnfd_single)

####### Niche differences
domin.mnd_single = glmer(domin~mnd+(1|field/f_p)+(1|fake_age) ,
                         family=binomial(link = "logit"),
                         data=dat_domin_com)
summary(domin.mnd_single)

domin.mnd.a_single = glmer(domin~mnd.a+(1|field/f_p)+(1|fake_age) ,
                           family=binomial(link = "logit"),
                           data=dat_domin_com)
summary(domin.mnd.a_single)

domin.mnnd_single = glmer(domin~mnnd+(1|field/f_p)+(1|fake_age) ,
                          family=binomial(link = "logit"),
                          data=dat_domin_com)
summary(domin.mnnd_single)

####### Relative fitness differences
domin.mrfd_single = glmer(domin~mrfd+(1|field/f_p)+(1|fake_age) ,
                          family=binomial(link = "logit"),
                          data=dat_domin_com)
summary(domin.mrfd_single)

domin.mrfd.a_single = glmer(domin~mrfd.a+(1|field/f_p)+(1|fake_age) ,
                            family=binomial(link = "logit"),
                            data=dat_domin_com)
summary(domin.mrfd.a_single)

domin.mnrfd_single = glmer(domin~mnrfd+(1|field/f_p)+(1|fake_age) ,
                           family=binomial(link = "logit"),
                           data=dat_domin_com)
summary(domin.mnrfd_single)



#### All stages
# Normal method
cor.test(dat_suc_sps$mpd, dat_suc_sps$mfd)
mod_sp_all_ab = clm(stage_level ~ mnd.a + mfd.a,
                    data = dat_suc_sps)
mod_sp_all_ab_1 = clmm(stage_level ~ mnd.a + mfd.a+ (1|species)+
                         (1|ffield/fplot),
                       data = dat_suc_sps)
summary(mod_sp_all_ab)
summary(mod_sp_all_ab_1)

# Bayesian method given the phylogenetic independence 
options(mc.cores = parallel::detectCores(logical = F)) 
mod_sp_all_md_brms = brm(stage_level ~ mnd + mfd + (1|species)+
                           (1|field/plot), data=dat_suc_sps,
                         cov_ranef = list(species = estab_vcv_tree), # phylo
                         family=cumulative("logit", link_disc = "log",
                                           threshold = "flexible"),
                         chains = 4,warmup = 1000, iter = 3000, refresh = 100,   
                         control = list(max_treedepth = 20,
                                        adapt_delta = 0.99), cores = 10)
save(mod_sp_all_md_brms,
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/mod_sp_all_md_brms.rdata')

mod_sp_all_md.a_brms = brm(stage_level ~ mnd.a + mfd.a + (1|species)+
                             (1|field/plot), data=dat_suc_sps,
                           cov_ranef = list(species = estab_vcv_tree), # phylo
                           family=cumulative("logit", link_disc = "log",
                                             threshold = "flexible"),
                           chains = 4,warmup = 1000, iter = 3000, refresh = 100,   
                           control = list(max_treedepth = 20,
                                          adapt_delta = 0.99), cores = 10)
save(mod_sp_all_md.a_brms,
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/mod_sp_all_md.a_brms.rdata')

mod_sp_all_mnd_brms = brm(stage_level ~ mnnd + mnfd + (1|species)+
                            (1|field/plot), data=dat_suc_sps,
                          cov_ranef = list(species = estab_vcv_tree), # phylo
                          family=cumulative("logit", link_disc = "log",
                                            threshold = "flexible"),
                          chains = 4,warmup = 1000, iter = 3000, refresh = 100,   
                          control = list(max_treedepth = 20,
                                         adapt_delta = 0.99), cores = 10)
save(mod_sp_all_mnd_brms,
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/mod_sp_all_mnd_brms.rdata')

mod_sp_all_md.a_fixed = fixef(mod_sp_all_md.a_brms,
                              probs = c(0.025, 0.975),
                              summary = T,
                              robust = T)
mod_sp_all_md_fixed = fixef(mod_sp_all_md_brms,
                            probs = c(0.025, 0.975),
                            summary = T,
                            robust = T)
mod_sp_all_mnd_fixed = fixef(mod_sp_all_mnd_brms,
                             probs = c(0.025, 0.975),
                             summary = T,
                             robust = T)


############### Fast start for analyzing mnd ~ mpd, mnd ~ mfunc_d, mfd ~ mpd, mfd ~ mfunc_d ###########
require(INLA)
require(phyr)
require(inlabru)
require(tibble)
require(lme4)
require(nlme)
require(lmerTest)
require(minpack.lm)

load('code/results_analyzing/analysing_ages1_35_top50_aposi_data/dat_suc_sp.rdata')

numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))
dat_suc_sps$species_1 = dat_suc_sps$species

max_e = function(a, b, dist){ 
  f = a*(1-exp(-(b*dist)))
  return(f)}
max_e_2 = function(b, dist){ 
  f = (1-exp(-(b*dist)))
  return(f)}
power = function(a, b, dist) {
  f = a*(dist^b)
  return(f)
}
hyperbola = function(a, b, dist) {
  f = dist/(a+b*dist)
  return(f)
}
logistic = function(a, b, c, dist) {
  f = a*dist/(1+exp(1)*((b-dist)/c))
}


####### mnd.a ~ mpd.a ######
mpd.a = dat_suc_sp$mpd.a
#plot(dist, max_e(1, 0.01, dist))
#plot(dist, max_e_2(0.002740836, dist))
#plot(dist, power(1, 0.01, dist))
#plot(dist, hyperbola(1, 0.01, dist))
#plot(dist, logistic(a = 2, b = 1, c = 5,
#                    dist))

model_e1 = nlsLM(mnd.a ~ max_e(a, b, dist = mpd.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnd.a ~ max_e_2(b, dist = mpd.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnd.a ~ power(a, b, dist = mpd.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnd.a ~ hyperbola(a, b, dist = mpd.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnd.a ~ logistic(a, b, c,
                                        dist = mpd.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

plot(mpd.a, max_e(nls_coff_e1[1,1],nls_coff_e1[2,1],mpd.a))
plot(mpd.a, max_e_2(nls_coff_e2[1,1],mpd.a))
#options(error=recover)

## lm
(mod_mnd.a_mpd.a_lmer = lmer(mnd.a ~ mpd.a + (1|field/f_p) + (1|species),
                             data = dat_suc_sp, REML = TRUE))
summary(mod_mnd.a_mpd.a_lmer)
(mod_mnd.a_mpd.a_lmer = lmer(mnd.a ~ mpd.a + (1|field/f_p) + (1|species),
                             data = dat_suc_sp, REML = TRUE))
summary(mod_mnd.a_mpd.a_lmer)
require(ggeffects)
ggpredict(mod_mnd.a_mpd.a_lmer, terms = 'mpd.a')
ggpredict(mod_mnd.a_mpd.a_lmer, terms = 'mpd.a')


(mod_mnd.a_mpd.a_pglmm = pglmm(mnd.a ~ mpd.a + (1|species) + (1|f_p) + (1|field),
                               data = dat_suc_sp,
                               family = "gaussian", cov_ranef = list(species = tree),
                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                           config = TRUE),
                                                    quantiles=c(0.025,0.5,0.975)),
                               bayes = T))
(mod_mnd.a_mpd.a_pglmm_nobayes = pglmm(mnd.a ~ mpd.a + (1|species) + (1|f_p) + (1|field),
                                       data = dat_suc_sp,
                                       family = "gaussian", cov_ranef = list(species = tree),
                                       bayes = F))
### mnd.a correlate with mpd.a positively 

#### predictive curve for mnd.a ~ mpd.a
lincombs.data.nd.a.pd.a = data.frame(mpd.a=seq(0.0001,max(dat_suc_sp$mpd.a),length=100))

lincombs.matrix.nd.a.pd.a=model.matrix(~mpd.a,
                                       data=lincombs.data.nd.a.pd.a)
lincombs.matrix.nd.a.pd.a=as.data.frame(lincombs.matrix.nd.a.pd.a)
lincombs.nd.a.pd.a=inla.make.lincombs(lincombs.matrix.nd.a.pd.a)

inla.model_lincombs.nd.a.pd.a = pglmm(mnd.a ~ mpd.a+(1|species) + 
                                        (1|f_p) + (1|field), data = dat_suc_sp,
                                      family = "gaussian", cov_ranef = list(species = tree),
                                      bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                  config = TRUE),
                                                           quantiles=c(0.025,0.5,0.975),
                                                           lincomb=lincombs.nd.a.pd.a,
                                                           control.predictor=list(compute=T)),
                                      bayes = T)


inla.model_lincombs.nd.a.pd.a$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nd.a.pd.a$predicted.value=inla.model_lincombs.nd.a.pd.a$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nd.a.pd.a$lower=inla.model_lincombs.nd.a.pd.a$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nd.a.pd.a$upper=inla.model_lincombs.nd.a.pd.a$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nd.a.pd.a

save(lincombs.data.nd.a.pd.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nd.a.pd.a.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnd.a_mpd.a_e = nlmer(mnd.a ~ max_e_fun(dist = mpd.a, a, b) ~ (a|field/f_p) + (a|species),
                           dat_suc_sp,
                           start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnd.a_mpd.a_e) 

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnd.a_mpd.a_e2 = nlmer(mnd.a ~ max_e2_fun(dist = mpd.a, b) ~ (b|field/f_p) + (b|species),
                            dat_suc_sp,
                            start = c(b = nls_coff_e2[1,1])))
summary(mod_mnd.a_mpd.a_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnd.a_mpd.a_power = nlmer(mnd.a ~ power_fun(dist = mpd.a, a, b) ~ (a|field/f_p) + (a|species),
                               dat_suc_sp,
                               start = c(a = nls_coff_power[1,1],
                                         b = nls_coff_power[2,1])))
summary(mod_mnd.a_mpd.a_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnd.a_mpd.a_hyperbola = nlmer(mnd.a ~ hyperbola_fun(dist = mpd.a, a, b) ~ a|field/f_p,
                                   dat_suc_sp,
                                   start = c(a = nls_coff_hyperbola[1,1],
                                             b = nls_coff_hyperbola[2,1])))
summary(mod_mnd.a_mpd.a_hyperbola)

anova(mod_mnd.a_mpd.a_e, mod_mnd.a_mpd.a_e2,
      mod_mnd.a_mpd.a_power, mod_mnd.a_mpd.a_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_mnd.a_mpd.a_e, mod_mnd.a_mpd.a_e2,
    mod_mnd.a_mpd.a_power, mod_mnd.a_mpd.a_hyperbola)



####### mfd.a ~ mpd.a ######
mpd.a = dat_suc_sp$mpd.a
#plot(dist, max_e(1, 0.01, dist))
#plot(dist, max_e_2(0.002740836, dist))
#plot(dist, power(1, 0.01, dist))
#plot(dist, hyperbola(1, 0.01, dist))
#plot(dist, logistic(a = 2, b = 1, c = 5,
#                    dist))

model_e1 = nlsLM(mfd.a ~ max_e(a, b, dist = mpd.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mfd.a ~ max_e_2(b, dist = mpd.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mfd.a ~ power(a, b, dist = mpd.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mfd.a ~ hyperbola(a, b, dist = mpd.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mfd.a ~ logistic(a, b, c,
                                          dist = mpd.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mfd.a_mpd.a_lmer = lmer(mfd.a ~ mpd.a + (1|field/f_p),
                               data = dat_suc_sp, REML = TRUE))
summary(mod_mfd.a_mpd.a_lmer)
ggpredict(mod_mfd.a_mpd.a_lmer, terms = 'mpd.a')

(mod_mfd.a_mpd.a_pglmm = pglmm(mfd.a ~ mpd.a + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                 family = "gaussian", cov_ranef = list(species = tree),
                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                             config = TRUE),
                                                      quantiles=c(0.025,0.5,0.975)),
                                 bayes = T))
mod_mfd.a_mpd.a_pglmm$inla.model$summary.fixed ### mfd.a correlate with mpd.a positively 

#### predictive curve for mfd.a ~ mpd.a
lincombs.data.fd.a.pd.a = data.frame(mpd.a=seq(0.0001,max(dat_suc_sp$mpd.a),length=100))

lincombs.matrix.fd.a.pd.a=model.matrix(~mpd.a,
                                         data=lincombs.data.fd.a.pd.a)
lincombs.matrix.fd.a.pd.a=as.data.frame(lincombs.matrix.fd.a.pd.a)
lincombs.fd.a.pd.a=inla.make.lincombs(lincombs.matrix.fd.a.pd.a)

inla.model_lincombs.fd.a.pd.a = pglmm(mfd.a ~  mpd.a+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                        family = "gaussian", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.fd.a.pd.a,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)


inla.model_lincombs.fd.a.pd.a$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.fd.a.pd.a$predicted.value=inla.model_lincombs.fd.a.pd.a$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.fd.a.pd.a$lower=inla.model_lincombs.fd.a.pd.a$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.fd.a.pd.a$upper=inla.model_lincombs.fd.a.pd.a$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.fd.a.pd.a

save(lincombs.data.fd.a.pd.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.fd.a.pd.a.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mfd.a_mpd.a_e = nlmer(mfd.a ~ max_e_fun(dist = mpd.a, a, b) ~ (a|field/f_p) + (a|species),
                             dat_suc_sp,
                             start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mfd.a_mpd.a_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mfd.a_mpd.a_e2 = nlmer(mfd.a ~ max_e2_fun(dist = mpd.a, b) ~ (b|field/f_p) + (b|species),
                              dat_suc_sp,
                              start = c(b = nls_coff_e2[1,1])))
summary(mod_mfd.a_mpd.a_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mfd.a_mpd.a_power = nlmer(mfd.a ~ power_fun(dist = mpd.a, a, b) ~ (a|field/f_p) + (a|species),
                                 dat_suc_sp,
                                 start = c(a = nls_coff_power[1,1],
                                           b = nls_coff_power[2,1])))
summary(mod_mfd.a_mpd.a_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mfd.a_mpd.a_hyperbola = nlmer(mfd.a ~ hyperbola_fun(dist = mpd.a, a, b) ~ a|field/f_p,
                                     dat_suc_sp,
                                     start = c(a = nls_coff_hyperbola[1,1],
                                               b = nls_coff_hyperbola[2,1])))
summary(mod_mfd.a_mpd.a_hyperbola)

anova(mod_lmer, mod_mfd.a_mpd.a_e, mod_mfd.a_mpd.a_e2,
      mod_mfd.a_mpd.a_power, mod_mfd.a_mpd.a_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)



####### mabfd.a ~ mpd.a ######
mpd.a = dat_suc_sp$mpd.a
#plot(dist, max_e(1, 0.01, dist))
#plot(dist, max_e_2(0.002740836, dist))
#plot(dist, power(1, 0.01, dist))
#plot(dist, hyperbola(1, 0.01, dist))
#plot(dist, logistic(a = 2, b = 1, c = 5,
#                    dist))

model_e1 = nlsLM(mabfd.a ~ max_e(a, b, dist = mpd.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mabfd.a ~ max_e_2(b, dist = mpd.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mabfd.a ~ power(a, b, dist = mpd.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mabfd.a ~ hyperbola(a, b, dist = mpd.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mabfd.a ~ logistic(a, b, c,
                                            dist = mpd.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mabfd.a_mpd.a_lmer = lmer(mabfd.a ~ mpd.a + (1|field/f_p),
                                 data = dat_suc_sp, REML = TRUE))
summary(mod_mabfd.a_mpd.a_lmer)
ggpredict(mod_mabfd.a_mpd.a_lmer, terms = 'mpd.a')

(mod_mabfd.a_mpd.a_pglmm = pglmm(mabfd.a ~ mpd.a + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                   family = "gaussian", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T))
mod_mabfd.a_mpd.a_pglmm$inla.model$summary.fixed ### mnd.a correlate with mpd.a positively 

#### predictive curve for mabfd.a ~ mpd.a
lincombs.data.abfd.a.pd.a = data.frame(mpd.a=seq(0.0001,max(dat_suc_sp$mpd.a),length=100))

lincombs.matrix.abfd.a.pd.a=model.matrix(~mpd.a,
                                           data=lincombs.data.abfd.a.pd.a)
lincombs.matrix.abfd.a.pd.a=as.data.frame(lincombs.matrix.abfd.a.pd.a)
lincombs.abfd.a.pd.a=inla.make.lincombs(lincombs.matrix.abfd.a.pd.a)

inla.model_lincombs.abfd.a.pd.a = pglmm(mabfd.a ~  mpd.a+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_suc_sp,
                                          family = "gaussian", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.abfd.a.pd.a,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)


inla.model_lincombs.abfd.a.pd.a$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.abfd.a.pd.a$predicted.value=inla.model_lincombs.abfd.a.pd.a$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.abfd.a.pd.a$lower=inla.model_lincombs.abfd.a.pd.a$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.abfd.a.pd.a$upper=inla.model_lincombs.abfd.a.pd.a$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.abfd.a.pd.a

save(lincombs.data.abfd.a.pd.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.abfd.a.pd.a.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mabfd.a_mpd.a_e = nlmer(mabfd.a ~ max_e_fun(dist = mpd.a, a, b) ~ (a|field/f_p) + (a|species),
                               dat_suc_sp,
                               start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mabfd.a_mpd.a_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mabfd.a_mpd.a_e2 = nlmer(mabfd.a ~ max_e2_fun(dist = mpd.a, b) ~ (b|field/f_p) + (b|species),
                                dat_suc_sp,
                                start = c(b = nls_coff_e2[1,1])))
summary(mod_mabfd.a_mpd.a_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mabfd.a_mpd.a_power = nlmer(mabfd.a ~ power_fun(dist = mpd.a, a, b) ~ (a|field/f_p) + (a|species),
                                   dat_suc_sp,
                                   start = c(a = nls_coff_power[1,1],
                                             b = nls_coff_power[2,1])))
summary(mod_mabfd.a_mpd.a_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mabfd.a_mpd.a_hyperbola = nlmer(mabfd.a ~ hyperbola_fun(dist = mpd.a, a, b) ~ a|field/f_p,
                                       dat_suc_sp,
                                       start = c(a = nls_coff_hyperbola[1,1],
                                                 b = nls_coff_hyperbola[2,1])))
summary(mod_mabfd.a_mpd.a_hyperbola)

anova(mod_lmer, mod_mabfd.a_mpd.a_e, mod_mabfd.a_mpd.a_e2,
      mod_mabfd.a_mpd.a_power, mod_mabfd.a_mpd.a_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)



####### mnd.a ~ mconti_func_d.a ######
mconti_func_d.a = dat_suc_sp$mconti_func_d.a
#plot(dist, max_e(1, 0.01, dist))
#plot(dist, max_e_2(0.002740836, dist))
#plot(dist, power(1, 0.01, dist))
#plot(dist, hyperbola(1, 0.01, dist))
#plot(dist, logistic(a = 2, b = 1, c = 5,
#                    dist))

model_e1 = nlsLM(mnd.a ~ max_e(a, b, dist = mconti_func_d.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnd.a ~ max_e_2(b, dist = mconti_func_d.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnd.a ~ power(a, b, dist = mconti_func_d.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnd.a ~ hyperbola(a, b, dist = mconti_func_d.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnd.a ~ logistic(a, b, c,
                                        dist = mconti_func_d.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnd.a_mconti_func_d.almer = lmer(mnd.a ~ mconti_func_d.a + (1|field/f_p),
                                      data = dat_suc_sp, REML = TRUE))
summary(mod_mnd.a_mconti_func_d.almer)

(mod_mnd.a_mconti_func_d.a_pglmm = pglmm(mnd.a ~ mconti_func_d.a + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                         family = "gaussian", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = TRUE),
                                                              quantiles=c(0.025,0.5,0.975)),
                                         bayes = T))
### mnd.a correlate with mconti_func_d.a positively 

#### predictive curve for mnd.a ~ mconti_func_d.a
lincombs.data.nd.a.conti_func_d.a = data.frame(mconti_func_d.a=seq(0.0001,
                                                                   max(dat_suc_sp$mconti_func_d.a),length=100))

lincombs.matrix.nd.a.conti_func_d.a=model.matrix(~mconti_func_d.a,
                                                 data=lincombs.data.nd.a.conti_func_d.a)
lincombs.matrix.nd.a.conti_func_d.a=as.data.frame(lincombs.matrix.nd.a.conti_func_d.a)
lincombs.nd.a.conti_func_d.a=inla.make.lincombs(lincombs.matrix.nd.a.conti_func_d.a)

inla.model_lincombs.nd.a.conti_func_d.a = pglmm(mnd.a ~  mconti_func_d.a+(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_suc_sp,
                                                family = "gaussian", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.nd.a.conti_func_d.a,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)


inla.model_lincombs.nd.a.conti_func_d.a$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nd.a.conti_func_d.a$predicted.value=inla.model_lincombs.nd.a.conti_func_d.a$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nd.a.conti_func_d.a$lower=inla.model_lincombs.nd.a.conti_func_d.a$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nd.a.conti_func_d.a$upper=inla.model_lincombs.nd.a.conti_func_d.a$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nd.a.conti_func_d.a

save(lincombs.data.nd.a.conti_func_d.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nd.a.conti_func_d.a.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnd.a_mconti_func_d.a_e = nlmer(mnd.a ~ max_e_fun(dist = mconti_func_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                     dat_suc_sp,
                                     start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnd.a_mconti_func_d.a_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnd.a_mconti_func_d.a_e2 = nlmer(mnd.a ~ max_e2_fun(dist = mconti_func_d.a, b) ~ (b|field/f_p) + (b|species),
                                      dat_suc_sp,
                                      start = c(b = nls_coff_e2[1,1])))
summary(mod_mnd.a_mconti_func_d.a_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnd.a_mconti_func_d.a_power = nlmer(mnd.a ~ power_fun(dist = mconti_func_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                         dat_suc_sp,
                                         start = c(a = nls_coff_power[1,1],
                                                   b = nls_coff_power[2,1])))
summary(mod_mnd.a_mconti_func_d.a_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnd.a_mconti_func_d.a_hyperbola = nlmer(mnd.a ~ hyperbola_fun(dist = mconti_func_d.a, a, b) ~ a|field/f_p,
                                             dat_suc_sp,
                                             start = c(a = nls_coff_hyperbola[1,1],
                                                       b = nls_coff_hyperbola[2,1])))
summary(mod_mnd.a_mconti_func_d.a_hyperbola)

anova(mod_lmer, mod_mnd.a_mconti_func_d.a_e, mod_mnd.a_mconti_func_d.a_e2,
      mod_mnd.a_mconti_func_d.a_power, mod_mnd.a_mconti_func_d.a_hyperbola, test="Chisq") ## mod_lmer the best!



####### mfd.a ~ mconti_func_d.a ######
mconti_func_d.a = dat_suc_sp$mconti_func_d.a
#plot(dist, max_e(1, 0.01, dist))
#plot(dist, max_e_2(0.002740836, dist))
#plot(dist, power(1, 0.01, dist))
#plot(dist, hyperbola(1, 0.01, dist))
#plot(dist, logistic(a = 2, b = 1, c = 5,
#                    dist))

model_e1 = nlsLM(mfd.a ~ max_e(a, b, dist = mconti_func_d.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mfd.a ~ max_e_2(b, dist = mconti_func_d.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mfd.a ~ power(a, b, dist = mconti_func_d.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mfd.a ~ hyperbola(a, b, dist = mconti_func_d.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mfd.a ~ logistic(a, b, c,
                                          dist = mconti_func_d.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mfd.a_mconti_func_d.almer = lmer(mfd.a ~ mconti_func_d.a + (1|field/f_p),
                                        data = dat_suc_sp, REML = TRUE))
summary(mod_mfd.a_mconti_func_d.almer)
ggpredict(mod_mfd.a_mconti_func_d.almer, terms = 'mconti_func_d.a')

(mod_mfd.a_mconti_func_d.a_pglmm = pglmm(mfd.a ~ mconti_func_d.a + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                           family = "gaussian", cov_ranef = list(species = tree),
                                           bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                       config = TRUE),
                                                                quantiles=c(0.025,0.5,0.975)),
                                           bayes = T))
### mfd.a correlate with mconti_func_d.a positively 

#### predictive curve for mfd.a ~ mconti_func_d.a
lincombs.data.fd.a.conti_func_d.a = data.frame(mconti_func_d.a=seq(0.0001,
                                                                     max(dat_suc_sp$mconti_func_d.a),length=100))

lincombs.matrix.fd.a.conti_func_d.a=model.matrix(~mconti_func_d.a,
                                                   data=lincombs.data.fd.a.conti_func_d.a)
lincombs.matrix.fd.a.conti_func_d.a=as.data.frame(lincombs.matrix.fd.a.conti_func_d.a)
lincombs.fd.a.conti_func_d.a=inla.make.lincombs(lincombs.matrix.fd.a.conti_func_d.a)

inla.model_lincombs.fd.a.conti_func_d.a = pglmm(mfd.a ~  mconti_func_d.a+(1|species) + 
                                                    (1|f_p) + (1|field), data = dat_suc_sp,
                                                  family = "gaussian", cov_ranef = list(species = tree),
                                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                              config = TRUE),
                                                                       quantiles=c(0.025,0.5,0.975),
                                                                       lincomb=lincombs.fd.a.conti_func_d.a,
                                                                       control.predictor=list(compute=T)),
                                                  bayes = T)

inla.model_lincombs.fd.a.conti_func_d.a$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.fd.a.conti_func_d.a$predicted.value=inla.model_lincombs.fd.a.conti_func_d.a$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.fd.a.conti_func_d.a$lower=inla.model_lincombs.fd.a.conti_func_d.a$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.fd.a.conti_func_d.a$upper=inla.model_lincombs.fd.a.conti_func_d.a$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.fd.a.conti_func_d.a

save(lincombs.data.fd.a.conti_func_d.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.fd.a.conti_func_d.a.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mfd.a_mconti_func_d.a_e = nlmer(mfd.a ~ max_e_fun(dist = mconti_func_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                       dat_suc_sp,
                                       start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mfd.a_mconti_func_d.a_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mfd.a_mconti_func_d.a_e2 = nlmer(mfd.a ~ max_e2_fun(dist = mconti_func_d.a, b) ~ (b|field/f_p) + (b|species),
                                        dat_suc_sp,
                                        start = c(b = nls_coff_e2[1,1])))
summary(mod_mfd.a_mconti_func_d.a_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mfd.a_mconti_func_d.a_power = nlmer(mfd.a ~ power_fun(dist = mconti_func_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                           dat_suc_sp,
                                           start = c(a = nls_coff_power[1,1],
                                                     b = nls_coff_power[2,1])))
summary(mod_mfd.a_mconti_func_d.a_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mfd.a_mconti_func_d.a_hyperbola = nlmer(mfd.a ~ hyperbola_fun(dist = mconti_func_d.a, a, b) ~ a|field/f_p,
                                               dat_suc_sp,
                                               start = c(a = nls_coff_hyperbola[1,1],
                                                         b = nls_coff_hyperbola[2,1])))
summary(mod_mfd.a_mconti_func_d.a_hyperbola)

anova(mod_lmer, mod_mfd.a_mconti_func_d.a_e, mod_mfd.a_mconti_func_d.a_e2,
      mod_mfd.a_mconti_func_d.a_power, mod_mfd.a_mconti_func_d.a_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)



####### mabfd.a ~ mconti_func_d.a ######
mconti_func_d.a = dat_suc_sp$mconti_func_d.a
#plot(dist, max_e(1, 0.01, dist))
#plot(dist, max_e_2(0.002740836, dist))
#plot(dist, power(1, 0.01, dist))
#plot(dist, hyperbola(1, 0.01, dist))
#plot(dist, logistic(a = 2, b = 1, c = 5,
#                    dist))

model_e1 = nlsLM(mabfd.a ~ max_e(a, b, dist = mconti_func_d.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mabfd.a ~ max_e_2(b, dist = mconti_func_d.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mabfd.a ~ power(a, b, dist = mconti_func_d.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mabfd.a ~ hyperbola(a, b, dist = mconti_func_d.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mabfd.a ~ logistic(a, b, c,
                                            dist = mconti_func_d.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mabfd.a_mconti_func_d.almer = lmer(mabfd.a ~ mconti_func_d.a + (1|field/f_p),
                                          data = dat_suc_sp, REML = TRUE))
summary(mod_mabfd.a_mconti_func_d.almer)
ggpredict(mod_mabfd.a_mconti_func_d.almer, terms = 'mconti_func_d.a')

(mod_mabfd.a_mconti_func_d.a_pglmm = pglmm(mabfd.a ~ mconti_func_d.a + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                             family = "gaussian", cov_ranef = list(species = tree),
                                             bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                         config = TRUE),
                                                                  quantiles=c(0.025,0.5,0.975)),
                                             bayes = T))
### mabfd.a correlate with mconti_func_d.a positively 

#### predictive curve for mabfd.a ~ mconti_func_d.a
lincombs.data.abfd.a.conti_func_d.a = data.frame(mconti_func_d.a=seq(0.0001,
                                                                       max(dat_suc_sp$mconti_func_d.a),length=100))

lincombs.matrix.abfd.a.conti_func_d.a=model.matrix(~mconti_func_d.a,
                                                     data=lincombs.data.abfd.a.conti_func_d.a)
lincombs.matrix.abfd.a.conti_func_d.a=as.data.frame(lincombs.matrix.abfd.a.conti_func_d.a)
lincombs.abfd.a.conti_func_d.a=inla.make.lincombs(lincombs.matrix.abfd.a.conti_func_d.a)

inla.model_lincombs.abfd.a.conti_func_d.a = pglmm(mabfd.a ~  mconti_func_d.a+(1|species) + 
                                                      (1|f_p) + (1|field), data = dat_suc_sp,
                                                    family = "gaussian", cov_ranef = list(species = tree),
                                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                config = TRUE),
                                                                         quantiles=c(0.025,0.5,0.975),
                                                                         lincomb=lincombs.abfd.a.conti_func_d.a,
                                                                         control.predictor=list(compute=T)),
                                                    bayes = T)

inla.model_lincombs.abfd.a.conti_func_d.a$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.abfd.a.conti_func_d.a$predicted.value=inla.model_lincombs.abfd.a.conti_func_d.a$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.abfd.a.conti_func_d.a$lower=inla.model_lincombs.abfd.a.conti_func_d.a$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.abfd.a.conti_func_d.a$upper=inla.model_lincombs.abfd.a.conti_func_d.a$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.abfd.a.conti_func_d.a

save(lincombs.data.abfd.a.conti_func_d.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.abfd.a.conti_func_d.a.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mabfd.a_mconti_func_d.a_e = nlmer(mabfd.a ~ max_e_fun(dist = mconti_func_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                         dat_suc_sp,
                                         start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mabfd.a_mconti_func_d.a_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mabfd.a_mconti_func_d.a_e2 = nlmer(mabfd.a ~ max_e2_fun(dist = mconti_func_d.a, b) ~ (b|field/f_p) + (b|species),
                                          dat_suc_sp,
                                          start = c(b = nls_coff_e2[1,1])))
summary(mod_mabfd.a_mconti_func_d.a_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mabfd.a_mconti_func_d.a_power = nlmer(mabfd.a ~ power_fun(dist = mconti_func_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                             dat_suc_sp,
                                             start = c(a = nls_coff_power[1,1],
                                                       b = nls_coff_power[2,1])))
summary(mod_mabfd.a_mconti_func_d.a_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mabfd.a_mconti_func_d.a_hyperbola = nlmer(mabfd.a ~ hyperbola_fun(dist = mconti_func_d.a, a, b) ~ a|field/f_p,
                                                 dat_suc_sp,
                                                 start = c(a = nls_coff_hyperbola[1,1],
                                                           b = nls_coff_hyperbola[2,1])))
summary(mod_mabfd.a_mconti_func_d.a_hyperbola)

anova(mod_lmer, mod_mabfd.a_mconti_func_d.a_e, mod_mabfd.a_mconti_func_d.a_e2,
      mod_mabfd.a_mconti_func_d.a_power, mod_mabfd.a_mconti_func_d.a_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)


####### mnd.a ~ mfunc_d.a ######
model_e1 = nlsLM(mnd.a ~ max_e(a, b, dist = mfunc_d.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnd.a ~ max_e_2(b, dist = mfunc_d.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnd.a ~ power(a, b, dist = mfunc_d.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnd.a ~ hyperbola(a, b, dist = mfunc_d.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnd.a ~ logistic(a, b, c,
                                        dist = mfunc_d.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnd.a_mfunc_d.almer = lmer(mnd.a ~ mfunc_d.a + (1|field/f_p),
                                data = dat_suc_sp, REML = TRUE))
summary(mod_mnd.a_mfunc_d.almer)

(mod_mnd.a_mfunc_d.a_pglmm = pglmm(mnd.a ~ mfunc_d.a + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                   family = "gaussian", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T))
### mnd.a correlate with mfunc_d.a positively 

#### predictive curve for mnd.a ~ mfunc_d.a
lincombs.data.nd.a.func_d.a = data.frame(mfunc_d.a=seq(0.0001,
                                                       max(dat_suc_sp$mfunc_d.a),length=100))

lincombs.matrix.nd.a.func_d.a=model.matrix(~mfunc_d.a,
                                           data=lincombs.data.nd.a.func_d.a)
lincombs.matrix.nd.a.func_d.a=as.data.frame(lincombs.matrix.nd.a.func_d.a)
lincombs.nd.a.func_d.a=inla.make.lincombs(lincombs.matrix.nd.a.func_d.a)

inla.model_lincombs.nd.a.func_d.a = pglmm(mnd.a ~  mfunc_d.a+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_suc_sp,
                                          family = "gaussian", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.nd.a.func_d.a,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)


inla.model_lincombs.nd.a.func_d.a$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nd.a.func_d.a$predicted.value=inla.model_lincombs.nd.a.func_d.a$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nd.a.func_d.a$lower=inla.model_lincombs.nd.a.func_d.a$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nd.a.func_d.a$upper=inla.model_lincombs.nd.a.func_d.a$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nd.a.func_d.a

save(lincombs.data.nd.a.func_d.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nd.a.func_d.a.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnd.a_mfunc_d.a_e = nlmer(mnd.a ~ max_e_fun(dist = mfunc_d.a, a, b) ~ (a|field/f_p) + (a|species),
                               dat_suc_sp,
                               start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnd.a_mfunc_d.a_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnd.a_mfunc_d.a_e2 = nlmer(mnd.a ~ max_e2_fun(dist = mfunc_d.a, b) ~ (b|field/f_p) + (b|species),
                                dat_suc_sp,
                                start = c(b = nls_coff_e2[1,1])))
summary(mod_mnd.a_mfunc_d.a_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnd.a_mfunc_d.a_power = nlmer(mnd.a ~ power_fun(dist = mfunc_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                   dat_suc_sp,
                                   start = c(a = nls_coff_power[1,1],
                                             b = nls_coff_power[2,1])))
summary(mod_mnd.a_mfunc_d.a_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnd.a_mfunc_d.a_hyperbola = nlmer(mnd.a ~ hyperbola_fun(dist = mfunc_d.a, a, b) ~ a|field/f_p,
                                       dat_suc_sp,
                                       start = c(a = nls_coff_hyperbola[1,1],
                                                 b = nls_coff_hyperbola[2,1])))
summary(mod_mnd.a_mfunc_d.a_hyperbola)

anova(mod_lmer, mod_mnd.a_mfunc_d.a_e, mod_mnd.a_mfunc_d.a_e2,
      mod_mnd.a_mfunc_d.a_power, mod_mnd.a_mfunc_d.a_hyperbola, test="Chisq") ## mod_lmer the best!



####### mfd.a ~ mfunc_d.a ######
model_e1 = nlsLM(mfd.a ~ max_e(a, b, dist = mfunc_d.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mfd.a ~ max_e_2(b, dist = mfunc_d.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mfd.a ~ power(a, b, dist = mfunc_d.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mfd.a ~ hyperbola(a, b, dist = mfunc_d.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mfd.a ~ logistic(a, b, c,
                                          dist = mfunc_d.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mfd.a_mfunc_d.almer = lmer(mfd.a ~ mfunc_d.a + (1|field/f_p),
                                  data = dat_suc_sp, REML = TRUE))
summary(mod_mfd.a_mfunc_d.almer)
ggpredict(mod_mfd.a_mfunc_d.almer, terms = 'mfunc_d.a')

(mod_mfd.a_mfunc_d.a_pglmm = pglmm(mfd.a ~ mfunc_d.a + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                     family = "gaussian", cov_ranef = list(species = tree),
                                     bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                 config = TRUE),
                                                          quantiles=c(0.025,0.5,0.975)),
                                     bayes = T))
### mfd.a correlate with mfunc_d.a positively 

#### predictive curve for mfd.a ~ mfunc_d.a
lincombs.data.fd.a.func_d.a = data.frame(mfunc_d.a=seq(0.0001,
                                                         max(dat_suc_sp$mfunc_d.a),
                                                         length=100))

lincombs.matrix.fd.a.func_d.a=model.matrix(~mfunc_d.a,
                                             data=lincombs.data.fd.a.func_d.a)
lincombs.matrix.fd.a.func_d.a=as.data.frame(lincombs.matrix.fd.a.func_d.a)
lincombs.fd.a.func_d.a=inla.make.lincombs(lincombs.matrix.fd.a.func_d.a)

inla.model_lincombs.fd.a.func_d.a = pglmm(mfd.a ~  mfunc_d.a+(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.fd.a.func_d.a,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.fd.a.func_d.a$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.fd.a.func_d.a$predicted.value=inla.model_lincombs.fd.a.func_d.a$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.fd.a.func_d.a$lower=inla.model_lincombs.fd.a.func_d.a$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.fd.a.func_d.a$upper=inla.model_lincombs.fd.a.func_d.a$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.fd.a.func_d.a

save(lincombs.data.fd.a.func_d.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.fd.a.func_d.a.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mfd.a_mfunc_d.a_e = nlmer(mfd.a ~ max_e_fun(dist = mfunc_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                 dat_suc_sp,
                                 start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mfd.a_mfunc_d.a_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mfd.a_mfunc_d.a_e2 = nlmer(mfd.a ~ max_e2_fun(dist = mfunc_d.a, b) ~ (b|field/f_p) + (b|species),
                                  dat_suc_sp,
                                  start = c(b = nls_coff_e2[1,1])))
summary(mod_mfd.a_mfunc_d.a_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mfd.a_mfunc_d.a_power = nlmer(mfd.a ~ power_fun(dist = mfunc_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                     dat_suc_sp,
                                     start = c(a = nls_coff_power[1,1],
                                               b = nls_coff_power[2,1])))
summary(mod_mfd.a_mfunc_d.a_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mfd.a_mfunc_d.a_hyperbola = nlmer(mfd.a ~ hyperbola_fun(dist = mfunc_d.a, a, b) ~ a|field/f_p,
                                         dat_suc_sp,
                                         start = c(a = nls_coff_hyperbola[1,1],
                                                   b = nls_coff_hyperbola[2,1])))
summary(mod_mfd.a_mfunc_d.a_hyperbola)

anova(mod_lmer, mod_mfd.a_mfunc_d.a_e, mod_mfd.a_mfunc_d.a_e2,
      mod_mfd.a_mfunc_d.a_power, mod_mfd.a_mfunc_d.a_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)

####### mabfd.a ~ mfunc_d.a ######
model_e1 = nlsLM(mabfd.a ~ max_e(a, b, dist = mfunc_d.a),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mabfd.a ~ max_e_2(b, dist = mfunc_d.a),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mabfd.a ~ power(a, b, dist = mfunc_d.a),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mabfd.a ~ hyperbola(a, b, dist = mfunc_d.a),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mabfd.a ~ logistic(a, b, c,
                                            dist = mfunc_d.a),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mabfd.a_mfunc_d.almer = lmer(mabfd.a ~ mfunc_d.a + (1|field/f_p),
                                    data = dat_suc_sp, REML = TRUE))
summary(mod_mabfd.a_mfunc_d.almer)
ggpredict(mod_mabfd.a_mfunc_d.almer, terms = 'mfunc_d.a')

(mod_mabfd.a_mfunc_d.a_pglmm = pglmm(mabfd.a ~ mfunc_d.a + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                       family = "gaussian", cov_ranef = list(species = tree),
                                       bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                   config = TRUE),
                                                            quantiles=c(0.025,0.5,0.975)),
                                       bayes = T))
### mabfd.a correlate with mfunc_d.a positively 

#### predictive curve for mabfd.a ~ mfunc_d.a 
lincombs.data.abfd.a.func_d.a = data.frame(mfunc_d.a=seq(0.0001,
                                                           max(dat_suc_sp$mfunc_d.a),length=100))

lincombs.matrix.abfd.a.func_d.a=model.matrix(~mfunc_d.a,
                                               data=lincombs.data.abfd.a.func_d.a)
lincombs.matrix.abfd.a.func_d.a=as.data.frame(lincombs.matrix.abfd.a.func_d.a)
lincombs.abfd.a.func_d.a=inla.make.lincombs(lincombs.matrix.abfd.a.func_d.a)

inla.model_lincombs.abfd.a.func_d.a = pglmm(mabfd.a ~  mfunc_d.a+(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "gaussian", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.abfd.a.func_d.a,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)

inla.model_lincombs.abfd.a.func_d.a$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.abfd.a.func_d.a$predicted.value=inla.model_lincombs.abfd.a.func_d.a$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.abfd.a.func_d.a$lower=inla.model_lincombs.abfd.a.func_d.a$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.abfd.a.func_d.a$upper=inla.model_lincombs.abfd.a.func_d.a$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.abfd.a.func_d.a

save(lincombs.data.abfd.a.func_d.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.abfd.a.func_d.a.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mabfd.a_mfunc_d.a_e = nlmer(mabfd.a ~ max_e_fun(dist = mfunc_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                   dat_suc_sp,
                                   start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mabfd.a_mfunc_d.a_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mabfd.a_mfunc_d.a_e2 = nlmer(mabfd.a ~ max_e2_fun(dist = mfunc_d.a, b) ~ (b|field/f_p) + (b|species),
                                    dat_suc_sp,
                                    start = c(b = nls_coff_e2[1,1])))
summary(mod_mabfd.a_mfunc_d.a_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mabfd.a_mfunc_d.a_power = nlmer(mabfd.a ~ power_fun(dist = mfunc_d.a, a, b) ~ (a|field/f_p) + (a|species),
                                       dat_suc_sp,
                                       start = c(a = nls_coff_power[1,1],
                                                 b = nls_coff_power[2,1])))
summary(mod_mabfd.a_mfunc_d.a_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mabfd.a_mfunc_d.a_hyperbola = nlmer(mabfd.a ~ hyperbola_fun(dist = mfunc_d.a, a, b) ~ a|field/f_p,
                                           dat_suc_sp,
                                           start = c(a = nls_coff_hyperbola[1,1],
                                                     b = nls_coff_hyperbola[2,1])))
summary(mod_mabfd.a_mfunc_d.a_hyperbola)

anova(mod_lmer, mod_mabfd.a_mfunc_d.a_e, mod_mabfd.a_mfunc_d.a_e2,
      mod_mabfd.a_mfunc_d.a_power, mod_mabfd.a_mfunc_d.a_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)


####### mnd ~ mpd ######
model_e1 = nlsLM(mnd ~ max_e(a, b, dist = mpd),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnd ~ max_e_2(b, dist = mpd),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnd ~ power(a, b, dist = mpd),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnd ~ hyperbola(a, b, dist = mpd),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnd ~ logistic(a, b, c,
                                      dist = mpd),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

plot(mpd, max_e(nls_coff_e1[1,1],nls_coff_e1[2,1],mpd))
plot(mpd, max_e_2(nls_coff_e2[1,1],mpd))
#options(error=recover)

## lm
(mod_mnd_mpd_lmer = lmer(mnd ~ mpd + (1|field/f_p) + (1|species),
                         data = dat_suc_sp, REML = TRUE))
summary(mod_mnd_mpd_lmer)
(mod_mnd_mpd_lmer = lmer(mnd ~ mpd + (1|field/f_p) + (1|species),
                         data = dat_suc_sp, REML = TRUE))
summary(mod_mnd_mpd_lmer)
require(ggeffects)
ggpredict(mod_mnd_mpd_lmer, terms = 'mpd')
ggpredict(mod_mnd_mpd_lmer, terms = 'mpd')


(mod_mnd_mpd_pglmm = pglmm(mnd ~ mpd + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                           family = "gaussian", cov_ranef = list(species = tree),
                           bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                       config = TRUE),
                                                quantiles=c(0.025,0.5,0.975)),
                           bayes = T))
mod_mnd_mpd_pglmm$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

### mnd correlate with mpd positively 

#### predictive curve for mnd ~ mpd 
lincombs.data.nd.pd = data.frame(mpd=seq(0.0001,max(dat_suc_sp$mpd),length=100))

lincombs.matrix.nd.pd=model.matrix(~mpd,
                                   data=lincombs.data.nd.pd)
lincombs.matrix.nd.pd=as.data.frame(lincombs.matrix.nd.pd)
lincombs.nd.pd=inla.make.lincombs(lincombs.matrix.nd.pd)

inla.model_lincombs.nd.pd = pglmm(mnd ~  mpd+(1|species) + 
                                    (1|f_p) + (1|field), data = dat_suc_sp,
                                  family = "gaussian", cov_ranef = list(species = tree),
                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                              config = TRUE),
                                                       quantiles=c(0.025,0.5,0.975),
                                                       lincomb=lincombs.nd.pd,
                                                       control.predictor=list(compute=T)),
                                  bayes = T)


inla.model_lincombs.nd.pd$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nd.pd$predicted.value=inla.model_lincombs.nd.pd$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nd.pd$lower=inla.model_lincombs.nd.pd$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nd.pd$upper=inla.model_lincombs.nd.pd$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nd.pd

save(lincombs.data.nd.pd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nd.pd.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnd_mpd_e = nlmer(mnd ~ max_e_fun(dist = mpd, a, b) ~ (a|field/f_p) + (a|species),
                       dat_suc_sp,
                       start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnd_mpd_e) 

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnd_mpd_e2 = nlmer(mnd ~ max_e2_fun(dist = mpd, b) ~ (b|field/f_p) + (b|species),
                        dat_suc_sp,
                        start = c(b = nls_coff_e2[1,1])))
summary(mod_mnd_mpd_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnd_mpd_power = nlmer(mnd ~ power_fun(dist = mpd, a, b) ~ (a|field/f_p) + (a|species),
                           dat_suc_sp,
                           start = c(a = nls_coff_power[1,1],
                                     b = nls_coff_power[2,1])))
summary(mod_mnd_mpd_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnd_mpd_hyperbola = nlmer(mnd ~ hyperbola_fun(dist = mpd, a, b) ~ a|field/f_p,
                               dat_suc_sp,
                               start = c(a = nls_coff_hyperbola[1,1],
                                         b = nls_coff_hyperbola[2,1])))
summary(mod_mnd_mpd_hyperbola)

anova(mod_mnd_mpd_e, mod_mnd_mpd_e2,
      mod_mnd_mpd_power, mod_mnd_mpd_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_mnd_mpd_e, mod_mnd_mpd_e2,
    mod_mnd_mpd_power, mod_mnd_mpd_hyperbola)


####### mfd ~ mpd ######
model_e1 = nlsLM(mfd ~ max_e(a, b, dist = mpd),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mfd ~ max_e_2(b, dist = mpd),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mfd ~ power(a, b, dist = mpd),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mfd ~ hyperbola(a, b, dist = mpd),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mfd ~ logistic(a, b, c,
                                        dist = mpd),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mfd_mpd_lmer = lmer(mfd ~ mpd + (1|field/f_p),
                           data = dat_suc_sp, REML = TRUE))
summary(mod_mfd_mpd_lmer)
ggpredict(mod_mfd_mpd_lmer, terms = 'mpd')

(mod_mfd_mpd_pglmm = pglmm(mfd ~ mpd + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                             family = "gaussian", cov_ranef = list(species = tree),
                             bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                         config = TRUE),
                                                  quantiles=c(0.025,0.5,0.975)),
                             bayes = T))
mod_mfd_mpd_pglmm$inla.model$summary.fixed ### mfd correlate with mpd 

#### predictive curve for mfd ~ mpd 
lincombs.data.fd.pd = data.frame(mpd=seq(0.0001,max(dat_suc_sp$mpd),length=100))

lincombs.matrix.fd.pd=model.matrix(~mpd,
                                     data=lincombs.data.fd.pd)
lincombs.matrix.fd.pd=as.data.frame(lincombs.matrix.fd.pd)
lincombs.fd.pd=inla.make.lincombs(lincombs.matrix.fd.pd)

inla.model_lincombs.fd.pd = pglmm(mfd ~  mpd+(1|species) + 
                                      (1|f_p) + (1|field), data = dat_suc_sp,
                                    family = "gaussian", cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975),
                                                         lincomb=lincombs.fd.pd,
                                                         control.predictor=list(compute=T)),
                                    bayes = T)


inla.model_lincombs.fd.pd$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.fd.pd$predicted.value=inla.model_lincombs.fd.pd$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.fd.pd$lower=inla.model_lincombs.fd.pd$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.fd.pd$upper=inla.model_lincombs.fd.pd$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.fd.pd

save(lincombs.data.fd.pd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.fd.pd.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mfd_mpd_e = nlmer(mfd ~ max_e_fun(dist = mpd, a, b) ~ (a|field/f_p) + (a|species),
                         dat_suc_sp,
                         start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mfd_mpd_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mfd_mpd_e2 = nlmer(mfd ~ max_e2_fun(dist = mpd, b) ~ (b|field/f_p) + (b|species),
                          dat_suc_sp,
                          start = c(b = nls_coff_e2[1,1])))
summary(mod_mfd_mpd_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mfd_mpd_power = nlmer(mfd ~ power_fun(dist = mpd, a, b) ~ (a|field/f_p) + (a|species),
                             dat_suc_sp,
                             start = c(a = nls_coff_power[1,1],
                                       b = nls_coff_power[2,1])))
summary(mod_mfd_mpd_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mfd_mpd_hyperbola = nlmer(mfd ~ hyperbola_fun(dist = mpd, a, b) ~ a|field/f_p,
                                 dat_suc_sp,
                                 start = c(a = nls_coff_hyperbola[1,1],
                                           b = nls_coff_hyperbola[2,1])))
summary(mod_mfd_mpd_hyperbola)

anova(mod_lmer, mod_mfd_mpd_e, mod_mfd_mpd_e2,
      mod_mfd_mpd_power, mod_mfd_mpd_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)


####### mabfd ~ mpd ######
model_e1 = nlsLM(mabfd ~ max_e(a, b, dist = mpd),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mabfd ~ max_e_2(b, dist = mpd),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mabfd ~ power(a, b, dist = mpd),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mabfd ~ hyperbola(a, b, dist = mpd),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mabfd ~ logistic(a, b, c,
                                          dist = mpd),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mabfd_mpd_lmer = lmer(mabfd ~ mpd + (1|field/f_p),
                             data = dat_suc_sp, REML = TRUE))
summary(mod_mabfd_mpd_lmer)
ggpredict(mod_mabfd_mpd_lmer, terms = 'mpd')

(mod_mabfd_mpd_pglmm = pglmm(mabfd ~ mpd + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                               family = "gaussian", cov_ranef = list(species = tree),
                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                           config = TRUE),
                                                    quantiles=c(0.025,0.5,0.975)),
                               bayes = T))
mod_mabfd_mpd_pglmm$inla.model$summary.fixed 

### mabfd correlate with mpd 

#### predictive curve for mabfd ~ mpd 
lincombs.data.abfd.pd = data.frame(mpd=seq(0.0001,max(dat_suc_sp$mpd),length=100))

lincombs.matrix.abfd.pd=model.matrix(~mpd,
                                       data=lincombs.data.abfd.pd)
lincombs.matrix.abfd.pd=as.data.frame(lincombs.matrix.abfd.pd)
lincombs.abfd.pd=inla.make.lincombs(lincombs.matrix.abfd.pd)

inla.model_lincombs.abfd.pd = pglmm(mabfd ~  mpd+(1|species) + 
                                        (1|f_p) + (1|field), data = dat_suc_sp,
                                      family = "gaussian", cov_ranef = list(species = tree),
                                      bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                  config = TRUE),
                                                           quantiles=c(0.025,0.5,0.975),
                                                           lincomb=lincombs.abfd.pd,
                                                           control.predictor=list(compute=T)),
                                      bayes = T)


inla.model_lincombs.abfd.pd$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.abfd.pd$predicted.value=inla.model_lincombs.abfd.pd$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.abfd.pd$lower=inla.model_lincombs.abfd.pd$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.abfd.pd$upper=inla.model_lincombs.abfd.pd$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.abfd.pd

save(lincombs.data.abfd.pd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.abfd.pd.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mabfd_mpd_e = nlmer(mabfd ~ max_e_fun(dist = mpd, a, b) ~ (a|field/f_p) + (a|species),
                           dat_suc_sp,
                           start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mabfd_mpd_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mabfd_mpd_e2 = nlmer(mabfd ~ max_e2_fun(dist = mpd, b) ~ (b|field/f_p) + (b|species),
                            dat_suc_sp,
                            start = c(b = nls_coff_e2[1,1])))
summary(mod_mabfd_mpd_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mabfd_mpd_power = nlmer(mabfd ~ power_fun(dist = mpd, a, b) ~ (a|field/f_p) + (a|species),
                               dat_suc_sp,
                               start = c(a = nls_coff_power[1,1],
                                         b = nls_coff_power[2,1])))
summary(mod_mabfd_mpd_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mabfd_mpd_hyperbola = nlmer(mabfd ~ hyperbola_fun(dist = mpd, a, b) ~ a|field/f_p,
                                   dat_suc_sp,
                                   start = c(a = nls_coff_hyperbola[1,1],
                                             b = nls_coff_hyperbola[2,1])))
summary(mod_mabfd_mpd_hyperbola)

anova(mod_lmer, mod_mabfd_mpd_e, mod_mabfd_mpd_e2,
      mod_mabfd_mpd_power, mod_mabfd_mpd_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)



####### mnd ~ mconti_func_d ######
model_e1 = nlsLM(mnd ~ max_e(a, b, dist = mconti_func_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnd ~ max_e_2(b, dist = mconti_func_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnd ~ power(a, b, dist = mconti_func_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnd ~ hyperbola(a, b, dist = mconti_func_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnd ~ logistic(a, b, c,
                                      dist = mconti_func_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnd_mconti_func_dlmer = lmer(mnd ~ mconti_func_d + (1|field/f_p),
                                  data = dat_suc_sp, REML = TRUE))
summary(mod_mnd_mconti_func_dlmer)

(mod_mnd_mconti_func_d_pglmm = pglmm(mnd ~ mconti_func_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                     family = "gaussian", cov_ranef = list(species = tree),
                                     bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                 config = TRUE),
                                                          quantiles=c(0.025,0.5,0.975)),
                                     bayes = T))
mod_mnd_mconti_func_d_pglmm$inla.model$summary.fixed[c(1,3,5)]%>%round(5) ## Extracting effects and confidence intervals for prediction curves from raw data

### mnd correlate with mconti_func_d 

#### predictive curve for mnd ~ mconti_func_d 
lincombs.data.nd.conti_func_d = data.frame(mconti_func_d=seq(0.0001,
                                                             max(dat_suc_sp$mconti_func_d),length=100))

lincombs.matrix.nd.conti_func_d=model.matrix(~mconti_func_d,
                                             data=lincombs.data.nd.conti_func_d)
lincombs.matrix.nd.conti_func_d=as.data.frame(lincombs.matrix.nd.conti_func_d)
lincombs.nd.conti_func_d=inla.make.lincombs(lincombs.matrix.nd.conti_func_d)

inla.model_lincombs.nd.conti_func_d = pglmm(mnd ~  mconti_func_d+(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.nd.conti_func_d,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)


inla.model_lincombs.nd.conti_func_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nd.conti_func_d$predicted.value=inla.model_lincombs.nd.conti_func_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nd.conti_func_d$lower=inla.model_lincombs.nd.conti_func_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nd.conti_func_d$upper=inla.model_lincombs.nd.conti_func_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nd.conti_func_d

save(lincombs.data.nd.conti_func_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nd.conti_func_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnd_mconti_func_d_e = nlmer(mnd ~ max_e_fun(dist = mconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                 dat_suc_sp,
                                 start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnd_mconti_func_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnd_mconti_func_d_e2 = nlmer(mnd ~ max_e2_fun(dist = mconti_func_d, b) ~ (b|field/f_p) + (b|species),
                                  dat_suc_sp,
                                  start = c(b = nls_coff_e2[1,1])))
summary(mod_mnd_mconti_func_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnd_mconti_func_d_power = nlmer(mnd ~ power_fun(dist = mconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                     dat_suc_sp,
                                     start = c(a = nls_coff_power[1,1],
                                               b = nls_coff_power[2,1])))
summary(mod_mnd_mconti_func_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnd_mconti_func_d_hyperbola = nlmer(mnd ~ hyperbola_fun(dist = mconti_func_d, a, b) ~ a|field/f_p,
                                         dat_suc_sp,
                                         start = c(a = nls_coff_hyperbola[1,1],
                                                   b = nls_coff_hyperbola[2,1])))
summary(mod_mnd_mconti_func_d_hyperbola)

anova(mod_lmer, mod_mnd_mconti_func_d_e, mod_mnd_mconti_func_d_e2,
      mod_mnd_mconti_func_d_power, mod_mnd_mconti_func_d_hyperbola, test="Chisq") ## mod_lmer the best!



####### mfd ~ mconti_func_d ######
model_e1 = nlsLM(mfd ~ max_e(a, b, dist = mconti_func_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mfd ~ max_e_2(b, dist = mconti_func_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mfd ~ power(a, b, dist = mconti_func_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mfd ~ hyperbola(a, b, dist = mconti_func_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mfd ~ logistic(a, b, c,
                                        dist = mconti_func_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mfd_mconti_func_dlmer = lmer(mfd ~ mconti_func_d + (1|field/f_p),
                                    data = dat_suc_sp, REML = TRUE))
summary(mod_mfd_mconti_func_dlmer)
ggpredict(mod_mfd_mconti_func_dlmer, terms = 'mconti_func_d')

(mod_mfd_mconti_func_d_pglmm = pglmm(mfd ~ mconti_func_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                       family = "gaussian", cov_ranef = list(species = tree),
                                       bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                   config = TRUE),
                                                            quantiles=c(0.025,0.5,0.975)),
                                       bayes = T))
### mfd correlate with mconti_func_d 

#### predictive curve for mfd ~ mconti_func_d 
lincombs.data.fd.conti_func_d = data.frame(mconti_func_d=seq(0.0001,
                                                               max(dat_suc_sp$mconti_func_d),length=100))

lincombs.matrix.fd.conti_func_d=model.matrix(~mconti_func_d,
                                               data=lincombs.data.fd.conti_func_d)
lincombs.matrix.fd.conti_func_d=as.data.frame(lincombs.matrix.fd.conti_func_d)
lincombs.fd.conti_func_d=inla.make.lincombs(lincombs.matrix.fd.conti_func_d)

inla.model_lincombs.fd.conti_func_d = pglmm(mfd ~  mconti_func_d+(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "gaussian", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.fd.conti_func_d,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)

inla.model_lincombs.fd.conti_func_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.fd.conti_func_d$predicted.value=inla.model_lincombs.fd.conti_func_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.fd.conti_func_d$lower=inla.model_lincombs.fd.conti_func_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.fd.conti_func_d$upper=inla.model_lincombs.fd.conti_func_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.fd.conti_func_d

save(lincombs.data.fd.conti_func_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.fd.conti_func_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mfd_mconti_func_d_e = nlmer(mfd ~ max_e_fun(dist = mconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                   dat_suc_sp,
                                   start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mfd_mconti_func_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mfd_mconti_func_d_e2 = nlmer(mfd ~ max_e2_fun(dist = mconti_func_d, b) ~ (b|field/f_p) + (b|species),
                                    dat_suc_sp,
                                    start = c(b = nls_coff_e2[1,1])))
summary(mod_mfd_mconti_func_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mfd_mconti_func_d_power = nlmer(mfd ~ power_fun(dist = mconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                       dat_suc_sp,
                                       start = c(a = nls_coff_power[1,1],
                                                 b = nls_coff_power[2,1])))
summary(mod_mfd_mconti_func_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mfd_mconti_func_d_hyperbola = nlmer(mfd ~ hyperbola_fun(dist = mconti_func_d, a, b) ~ a|field/f_p,
                                           dat_suc_sp,
                                           start = c(a = nls_coff_hyperbola[1,1],
                                                     b = nls_coff_hyperbola[2,1])))
summary(mod_mfd_mconti_func_d_hyperbola)

anova(mod_lmer, mod_mfd_mconti_func_d_e, mod_mfd_mconti_func_d_e2,
      mod_mfd_mconti_func_d_power, mod_mfd_mconti_func_d_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)


####### mabfd ~ mconti_func_d ######
model_e1 = nlsLM(mabfd ~ max_e(a, b, dist = mconti_func_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mabfd ~ max_e_2(b, dist = mconti_func_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mabfd ~ power(a, b, dist = mconti_func_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mabfd ~ hyperbola(a, b, dist = mconti_func_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mabfd ~ logistic(a, b, c,
                                          dist = mconti_func_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mabfd_mconti_func_dlmer = lmer(mabfd ~ mconti_func_d + (1|field/f_p),
                                      data = dat_suc_sp, REML = TRUE))
summary(mod_mabfd_mconti_func_dlmer)
ggpredict(mod_mabfd_mconti_func_dlmer, terms = 'mconti_func_d')

(mod_mabfd_mconti_func_d_pglmm = pglmm(mabfd ~ mconti_func_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                         family = "gaussian", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = TRUE),
                                                              quantiles=c(0.025,0.5,0.975)),
                                         bayes = T))
mod_mabfd_mconti_func_d_pglmm$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

### mabfd correlate with mconti_func_d 

#### predictive curve for mnd ~ mconti_func_d 
lincombs.data.abfd.conti_func_d = data.frame(mconti_func_d=seq(0.0001,
                                                                 max(dat_suc_sp$mconti_func_d),length=100))

lincombs.matrix.abfd.conti_func_d=model.matrix(~mconti_func_d,
                                                 data=lincombs.data.abfd.conti_func_d)
lincombs.matrix.abfd.conti_func_d=as.data.frame(lincombs.matrix.abfd.conti_func_d)
lincombs.abfd.conti_func_d=inla.make.lincombs(lincombs.matrix.abfd.conti_func_d)

inla.model_lincombs.abfd.conti_func_d = pglmm(mabfd ~  mconti_func_d+(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_suc_sp,
                                                family = "gaussian", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.abfd.conti_func_d,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

inla.model_lincombs.abfd.conti_func_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.abfd.conti_func_d$predicted.value=inla.model_lincombs.abfd.conti_func_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.abfd.conti_func_d$lower=inla.model_lincombs.abfd.conti_func_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.abfd.conti_func_d$upper=inla.model_lincombs.abfd.conti_func_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.abfd.conti_func_d

save(lincombs.data.abfd.conti_func_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.abfd.conti_func_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mabfd_mconti_func_d_e = nlmer(mabfd ~ max_e_fun(dist = mconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                     dat_suc_sp,
                                     start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mabfd_mconti_func_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mabfd_mconti_func_d_e2 = nlmer(mabfd ~ max_e2_fun(dist = mconti_func_d, b) ~ (b|field/f_p) + (b|species),
                                      dat_suc_sp,
                                      start = c(b = nls_coff_e2[1,1])))
summary(mod_mabfd_mconti_func_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mabfd_mconti_func_d_power = nlmer(mabfd ~ power_fun(dist = mconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                         dat_suc_sp,
                                         start = c(a = nls_coff_power[1,1],
                                                   b = nls_coff_power[2,1])))
summary(mod_mabfd_mconti_func_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mabfd_mconti_func_d_hyperbola = nlmer(mabfd ~ hyperbola_fun(dist = mconti_func_d, a, b) ~ a|field/f_p,
                                             dat_suc_sp,
                                             start = c(a = nls_coff_hyperbola[1,1],
                                                       b = nls_coff_hyperbola[2,1])))
summary(mod_mabfd_mconti_func_d_hyperbola)

anova(mod_lmer, mod_mabfd_mconti_func_d_e, mod_mabfd_mconti_func_d_e2,
      mod_mabfd_mconti_func_d_power, mod_mabfd_mconti_func_d_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)


####### mnd ~ mfunc_d ######
model_e1 = nlsLM(mnd ~ max_e(a, b, dist = mfunc_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnd ~ max_e_2(b, dist = mfunc_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnd ~ power(a, b, dist = mfunc_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnd ~ hyperbola(a, b, dist = mfunc_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnd ~ logistic(a, b, c,
                                      dist = mfunc_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnd_mfunc_dlmer = lmer(mnd ~ mfunc_d + (1|field/f_p),
                            data = dat_suc_sp, REML = TRUE))
summary(mod_mnd_mfunc_dlmer)

(mod_mnd_mfunc_d_pglmm = pglmm(mnd ~ mfunc_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                               family = "gaussian", cov_ranef = list(species = tree),
                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                           config = TRUE),
                                                    quantiles=c(0.025,0.5,0.975)),
                               bayes = T))
### mnd correlate with mfunc_d 

#### predictive curve for mnd ~ mfunc_d 
lincombs.data.nd.func_d = data.frame(mfunc_d=seq(0.0001,
                                                 max(dat_suc_sp$mfunc_d),length=100))

lincombs.matrix.nd.func_d=model.matrix(~mfunc_d,
                                       data=lincombs.data.nd.func_d)
lincombs.matrix.nd.func_d=as.data.frame(lincombs.matrix.nd.func_d)
lincombs.nd.func_d=inla.make.lincombs(lincombs.matrix.nd.func_d)

inla.model_lincombs.nd.func_d = pglmm(mnd ~  mfunc_d+(1|species) + 
                                        (1|f_p) + (1|field), data = dat_suc_sp,
                                      family = "gaussian", cov_ranef = list(species = tree),
                                      bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                  config = TRUE),
                                                           quantiles=c(0.025,0.5,0.975),
                                                           lincomb=lincombs.nd.func_d,
                                                           control.predictor=list(compute=T)),
                                      bayes = T)


inla.model_lincombs.nd.func_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nd.func_d$predicted.value=inla.model_lincombs.nd.func_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nd.func_d$lower=inla.model_lincombs.nd.func_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nd.func_d$upper=inla.model_lincombs.nd.func_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nd.func_d

save(lincombs.data.nd.func_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nd.func_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnd_mfunc_d_e = nlmer(mnd ~ max_e_fun(dist = mfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                           dat_suc_sp,
                           start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnd_mfunc_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnd_mfunc_d_e2 = nlmer(mnd ~ max_e2_fun(dist = mfunc_d, b) ~ (b|field/f_p) + (b|species),
                            dat_suc_sp,
                            start = c(b = nls_coff_e2[1,1])))
summary(mod_mnd_mfunc_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnd_mfunc_d_power = nlmer(mnd ~ power_fun(dist = mfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                               dat_suc_sp,
                               start = c(a = nls_coff_power[1,1],
                                         b = nls_coff_power[2,1])))
summary(mod_mnd_mfunc_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnd_mfunc_d_hyperbola = nlmer(mnd ~ hyperbola_fun(dist = mfunc_d, a, b) ~ a|field/f_p,
                                   dat_suc_sp,
                                   start = c(a = nls_coff_hyperbola[1,1],
                                             b = nls_coff_hyperbola[2,1])))
summary(mod_mnd_mfunc_d_hyperbola)

anova(mod_lmer, mod_mnd_mfunc_d_e, mod_mnd_mfunc_d_e2,
      mod_mnd_mfunc_d_power, mod_mnd_mfunc_d_hyperbola, test="Chisq") ## mod_lmer the best!



####### mfd ~ mfunc_d ######
model_e1 = nlsLM(mfd ~ max_e(a, b, dist = mfunc_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mfd ~ max_e_2(b, dist = mfunc_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mfd ~ power(a, b, dist = mfunc_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mfd ~ hyperbola(a, b, dist = mfunc_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mfd ~ logistic(a, b, c,
                                        dist = mfunc_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mfd_mfunc_dlmer = lmer(mfd ~ mfunc_d + (1|field/f_p),
                              data = dat_suc_sp, REML = TRUE))
summary(mod_mfd_mfunc_dlmer)
ggpredict(mod_mfd_mfunc_dlmer, terms = 'mfunc_d')

(mod_mfd_mfunc_d_pglmm = pglmm(mfd ~ mfunc_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                 family = "gaussian", cov_ranef = list(species = tree),
                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                             config = TRUE),
                                                      quantiles=c(0.025,0.5,0.975)),
                                 bayes = T))
### mabfd correlate with mfunc_d 

#### predictive curve for mabfd ~ mfunc_d 
lincombs.data.fd.func_d = data.frame(mfunc_d=seq(0.0001,
                                                   max(dat_suc_sp$mfunc_d),length=100))

lincombs.matrix.fd.func_d=model.matrix(~mfunc_d,
                                         data=lincombs.data.fd.func_d)
lincombs.matrix.fd.func_d=as.data.frame(lincombs.matrix.fd.func_d)
lincombs.fd.func_d=inla.make.lincombs(lincombs.matrix.fd.func_d)

inla.model_lincombs.fd.func_d = pglmm(mfd ~  mfunc_d+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                        family = "gaussian", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.fd.func_d,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)

inla.model_lincombs.fd.func_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.fd.func_d$predicted.value=inla.model_lincombs.fd.func_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.fd.func_d$lower=inla.model_lincombs.fd.func_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.fd.func_d$upper=inla.model_lincombs.fd.func_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.fd.func_d

save(lincombs.data.fd.func_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.fd.func_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mfd_mfunc_d_e = nlmer(mfd ~ max_e_fun(dist = mfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                             dat_suc_sp,
                             start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mfd_mfunc_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mfd_mfunc_d_e2 = nlmer(mfd ~ max_e2_fun(dist = mfunc_d, b) ~ (b|field/f_p) + (b|species),
                              dat_suc_sp,
                              start = c(b = nls_coff_e2[1,1])))
summary(mod_mfd_mfunc_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mfd_mfunc_d_power = nlmer(mfd ~ power_fun(dist = mfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                                 dat_suc_sp,
                                 start = c(a = nls_coff_power[1,1],
                                           b = nls_coff_power[2,1])))
summary(mod_mfd_mfunc_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mfd_mfunc_d_hyperbola = nlmer(mfd ~ hyperbola_fun(dist = mfunc_d, a, b) ~ a|field/f_p,
                                     dat_suc_sp,
                                     start = c(a = nls_coff_hyperbola[1,1],
                                               b = nls_coff_hyperbola[2,1])))
summary(mod_mfd_mfunc_d_hyperbola)

anova(mod_lmer, mod_mfd_mfunc_d_e, mod_mfd_mfunc_d_e2,
      mod_mfd_mfunc_d_power, mod_mfd_mfunc_d_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)


####### mabfd ~ mfunc_d ######
model_e1 = nlsLM(mabfd ~ max_e(a, b, dist = mfunc_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mabfd ~ max_e_2(b, dist = mfunc_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mabfd ~ power(a, b, dist = mfunc_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mabfd ~ hyperbola(a, b, dist = mfunc_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mabfd ~ logistic(a, b, c,
                                          dist = mfunc_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mabfd_mfunc_dlmer = lmer(mabfd ~ mfunc_d + (1|field/f_p),
                                data = dat_suc_sp, REML = TRUE))
summary(mod_mabfd_mfunc_dlmer)
ggpredict(mod_mabfd_mfunc_dlmer, terms = 'mfunc_d')

(mod_mabfd_mfunc_d_pglmm = pglmm(mabfd ~ mfunc_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                   family = "gaussian", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T))
### mabfd correlate with mfunc_d 

#### predictive curve for mabfd ~ mfunc_d 
lincombs.data.abfd.func_d = data.frame(mfunc_d=seq(0.0001,
                                                     max(dat_suc_sp$mfunc_d),length=100))

lincombs.matrix.abfd.func_d=model.matrix(~mfunc_d,
                                           data=lincombs.data.abfd.func_d)
lincombs.matrix.abfd.func_d=as.data.frame(lincombs.matrix.abfd.func_d)
lincombs.abfd.func_d=inla.make.lincombs(lincombs.matrix.abfd.func_d)

inla.model_lincombs.abfd.func_d = pglmm(mabfd ~  mfunc_d+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_suc_sp,
                                          family = "gaussian", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.abfd.func_d,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)

inla.model_lincombs.abfd.func_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.abfd.func_d$predicted.value=inla.model_lincombs.abfd.func_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.abfd.func_d$lower=inla.model_lincombs.abfd.func_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.abfd.func_d$upper=inla.model_lincombs.abfd.func_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.abfd.func_d

save(lincombs.data.abfd.func_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.abfd.func_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mabfd_mfunc_d_e = nlmer(mabfd ~ max_e_fun(dist = mfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                               dat_suc_sp,
                               start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mabfd_mfunc_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mabfd_mfunc_d_e2 = nlmer(mabfd ~ max_e2_fun(dist = mfunc_d, b) ~ (b|field/f_p) + (b|species),
                                dat_suc_sp,
                                start = c(b = nls_coff_e2[1,1])))
summary(mod_mabfd_mfunc_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mabfd_mfunc_d_power = nlmer(mabfd ~ power_fun(dist = mfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                                   dat_suc_sp,
                                   start = c(a = nls_coff_power[1,1],
                                             b = nls_coff_power[2,1])))
summary(mod_mabfd_mfunc_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mabfd_mfunc_d_hyperbola = nlmer(mabfd ~ hyperbola_fun(dist = mfunc_d, a, b) ~ a|field/f_p,
                                       dat_suc_sp,
                                       start = c(a = nls_coff_hyperbola[1,1],
                                                 b = nls_coff_hyperbola[2,1])))
summary(mod_mabfd_mfunc_d_hyperbola)

anova(mod_lmer, mod_mabfd_mfunc_d_e, mod_mabfd_mfunc_d_e2,
      mod_mabfd_mfunc_d_power, mod_mabfd_mfunc_d_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)


####### mnnd ~ mntd ######
model_e1 = nlsLM(mnnd ~ max_e(a, b, dist = mntd),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnnd ~ max_e_2(b, dist = mntd),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnnd ~ power(a, b, dist = mntd),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnnd ~ hyperbola(a, b, dist = mntd),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnnd ~ logistic(a, b, c,
                                       dist = mntd),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

plot(mntd, max_e(nls_coff_e1[1,1],nls_coff_e1[2,1],mntd))
plot(mntd, max_e_2(nls_coff_e2[1,1],mntd))
#options(error=recover)

## lm
(mod_mnnd_mntd_lmer = lmer(mnnd ~ mntd + (1|field/f_p) + (1|species),
                           data = dat_suc_sp, REML = TRUE))
summary(mod_mnnd_mntd_lmer)
(mod_mnnd_mntd_lmer = lmer(mnnd ~ mntd + (1|field/f_p) + (1|species),
                           data = dat_suc_sp, REML = TRUE))
summary(mod_mnnd_mntd_lmer)
require(ggeffects)
ggpredict(mod_mnnd_mntd_lmer, terms = 'mntd')
ggpredict(mod_mnnd_mntd_lmer, terms = 'mntd')


(mod_mnnd_mntd_pglmm = pglmm(mnnd ~ mntd + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                             family = "gaussian", cov_ranef = list(species = tree),
                             bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                         config = TRUE),
                                                  quantiles=c(0.025,0.5,0.975)),
                             bayes = T))
### mnnd correlate with mntd 

#### predictive curve for mnnd ~ mntd 
lincombs.data.nnd.ntd = data.frame(mntd=seq((0.0001),max(dat_suc_sp$mntd),length=100))

lincombs.matrix.nnd.ntd=model.matrix(~mntd,
                                     data=lincombs.data.nnd.ntd)
lincombs.matrix.nnd.ntd=as.data.frame(lincombs.matrix.nnd.ntd)
lincombs.nnd.ntd=inla.make.lincombs(lincombs.matrix.nnd.ntd)

inla.model_lincombs.nnd.ntd = pglmm(mnnd ~  mntd+(1|species) + 
                                      (1|f_p) + (1|field), data = dat_suc_sp,
                                    family = "gaussian", cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975),
                                                         lincomb=lincombs.nnd.ntd,
                                                         control.predictor=list(compute=T)),
                                    bayes = T)


inla.model_lincombs.nnd.ntd$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nnd.ntd$predicted.value=inla.model_lincombs.nnd.ntd$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nnd.ntd$lower=inla.model_lincombs.nnd.ntd$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nnd.ntd$upper=inla.model_lincombs.nnd.ntd$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nnd.ntd

save(lincombs.data.nnd.ntd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nnd.ntd.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnnd_mntd_e = nlmer(mnnd ~ max_e_fun(dist = mntd, a, b) ~ (a|field/f_p) + (a|species),
                         dat_suc_sp,
                         start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnnd_mntd_e) 

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnnd_mntd_e2 = nlmer(mnnd ~ max_e2_fun(dist = mntd, b) ~ (b|field/f_p) + (b|species),
                          dat_suc_sp,
                          start = c(b = nls_coff_e2[1,1])))
summary(mod_mnnd_mntd_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnnd_mntd_power = nlmer(mnnd ~ power_fun(dist = mntd, a, b) ~ (a|field/f_p) + (a|species),
                             dat_suc_sp,
                             start = c(a = nls_coff_power[1,1],
                                       b = nls_coff_power[2,1])))
summary(mod_mnnd_mntd_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnnd_mntd_hyperbola = nlmer(mnnd ~ hyperbola_fun(dist = mntd, a, b) ~ a|field/f_p,
                                 dat_suc_sp,
                                 start = c(a = nls_coff_hyperbola[1,1],
                                           b = nls_coff_hyperbola[2,1])))
summary(mod_mnnd_mntd_hyperbola)

anova(mod_mnnd_mntd_e, mod_mnnd_mntd_e2,
      mod_mnnd_mntd_power, mod_mnnd_mntd_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_mnnd_mntd_e, mod_mnnd_mntd_e2,
    mod_mnnd_mntd_power, mod_mnnd_mntd_hyperbola)


####### mnfd ~ mntd ######
model_e1 = nlsLM(mnfd ~ max_e(a, b, dist = mntd),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnfd ~ max_e_2(b, dist = mntd),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnfd ~ power(a, b, dist = mntd),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnfd ~ hyperbola(a, b, dist = mntd),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnfd ~ logistic(a, b, c,
                                         dist = mntd),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnfd_mntd_lmer = lmer(mnfd ~ mntd + (1|field/f_p),
                             data = dat_suc_sp, REML = TRUE))
summary(mod_mnfd_mntd_lmer)
ggpredict(mod_mnfd_mntd_lmer, terms = 'mntd')

(mod_mnfd_mntd_pglmm = pglmm(mnfd ~ mntd + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                               family = "gaussian", cov_ranef = list(species = tree),
                               bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                           config = TRUE),
                                                    quantiles=c(0.025,0.5,0.975)),
                               bayes = T))
mod_mnfd_mntd_pglmm$inla.model$summary.fixed

### mnfd correlate with mntd 

#### predictive curve for mnfd ~ mntd 
lincombs.data.nfd.ntd = data.frame(mntd=seq(0.0001,max(dat_suc_sp$mntd),length=100))

lincombs.matrix.nfd.ntd=model.matrix(~mntd,
                                       data=lincombs.data.nfd.ntd)
lincombs.matrix.nfd.ntd=as.data.frame(lincombs.matrix.nfd.ntd)
lincombs.nfd.ntd=inla.make.lincombs(lincombs.matrix.nfd.ntd)

inla.model_lincombs.nfd.ntd = pglmm(mnfd ~  mntd+(1|species) + 
                                        (1|f_p) + (1|field), data = dat_suc_sp,
                                      family = "gaussian", cov_ranef = list(species = tree),
                                      bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                  config = TRUE),
                                                           quantiles=c(0.025,0.5,0.975),
                                                           lincomb=lincombs.nfd.ntd,
                                                           control.predictor=list(compute=T)),
                                      bayes = T)


inla.model_lincombs.nfd.ntd$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nfd.ntd$predicted.value=inla.model_lincombs.nfd.ntd$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nfd.ntd$lower=inla.model_lincombs.nfd.ntd$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nfd.ntd$upper=inla.model_lincombs.nfd.ntd$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nfd.ntd

save(lincombs.data.nfd.ntd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nfd.ntd.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnfd_mntd_e = nlmer(mnfd ~ max_e_fun(dist = mntd, a, b) ~ (a|field/f_p) + (a|species),
                           dat_suc_sp,
                           start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnfd_mntd_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnfd_mntd_e2 = nlmer(mnfd ~ max_e2_fun(dist = mntd, b) ~ (b|field/f_p) + (b|species),
                            dat_suc_sp,
                            start = c(b = nls_coff_e2[1,1])))
summary(mod_mnfd_mntd_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnfd_mntd_power = nlmer(mnfd ~ power_fun(dist = mntd, a, b) ~ (a|field/f_p) + (a|species),
                               dat_suc_sp,
                               start = c(a = nls_coff_power[1,1],
                                         b = nls_coff_power[2,1])))
summary(mod_mnfd_mntd_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnfd_mntd_hyperbola = nlmer(mnfd ~ hyperbola_fun(dist = mntd, a, b) ~ a|field/f_p,
                                   dat_suc_sp,
                                   start = c(a = nls_coff_hyperbola[1,1],
                                             b = nls_coff_hyperbola[2,1])))
summary(mod_mnfd_mntd_hyperbola)

anova(mod_lmer, mod_mnfd_mntd_e, mod_mnfd_mntd_e2,
      mod_mnfd_mntd_power, mod_mnfd_mntd_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)




####### mnabfd ~ mntd ######
model_e1 = nlsLM(mnabfd ~ max_e(a, b, dist = mntd),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnabfd ~ max_e_2(b, dist = mntd),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnabfd ~ power(a, b, dist = mntd),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnabfd ~ hyperbola(a, b, dist = mntd),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnabfd ~ logistic(a, b, c,
                                           dist = mntd),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnabfd_mntd_lmer = lmer(mnabfd ~ mntd + (1|field/f_p),
                               data = dat_suc_sp, REML = TRUE))
summary(mod_mnabfd_mntd_lmer)
ggpredict(mod_mnabfd_mntd_lmer, terms = 'mntd')

(mod_mnabfd_mntd_pglmm = pglmm(mnabfd ~ mntd + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                 family = "gaussian", cov_ranef = list(species = tree),
                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                             config = TRUE),
                                                      quantiles=c(0.025,0.5,0.975)),
                                 bayes = T))
mod_mnabfd_mntd_pglmm$inla.model$summary.fixed

### nabfd correlate with mntd 

#### predictive curve for nabfd ~ mntd 
lincombs.data.nabfd.ntd = data.frame(mntd=seq(0.0001,max(dat_suc_sp$mntd),length=100))

lincombs.matrix.nabfd.ntd=model.matrix(~mntd,
                                         data=lincombs.data.nabfd.ntd)
lincombs.matrix.nabfd.ntd=as.data.frame(lincombs.matrix.nabfd.ntd)
lincombs.nabfd.ntd=inla.make.lincombs(lincombs.matrix.nabfd.ntd)

inla.model_lincombs.nabfd.ntd = pglmm(mnabfd ~  mntd+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                        family = "gaussian", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.nabfd.ntd,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)


inla.model_lincombs.nabfd.ntd$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nabfd.ntd$predicted.value=inla.model_lincombs.nabfd.ntd$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nabfd.ntd$lower=inla.model_lincombs.nabfd.ntd$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nabfd.ntd$upper=inla.model_lincombs.nabfd.ntd$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nabfd.ntd

save(lincombs.data.nabfd.ntd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nabfd.ntd.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnabfd_mntd_e = nlmer(mnabfd ~ max_e_fun(dist = mntd, a, b) ~ (a|field/f_p) + (a|species),
                             dat_suc_sp,
                             start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnabfd_mntd_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnabfd_mntd_e2 = nlmer(mnabfd ~ max_e2_fun(dist = mntd, b) ~ (b|field/f_p) + (b|species),
                              dat_suc_sp,
                              start = c(b = nls_coff_e2[1,1])))
summary(mod_mnabfd_mntd_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnabfd_mntd_power = nlmer(mnabfd ~ power_fun(dist = mntd, a, b) ~ (a|field/f_p) + (a|species),
                                 dat_suc_sp,
                                 start = c(a = nls_coff_power[1,1],
                                           b = nls_coff_power[2,1])))
summary(mod_mnabfd_mntd_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnabfd_mntd_hyperbola = nlmer(mnabfd ~ hyperbola_fun(dist = mntd, a, b) ~ a|field/f_p,
                                     dat_suc_sp,
                                     start = c(a = nls_coff_hyperbola[1,1],
                                               b = nls_coff_hyperbola[2,1])))
summary(mod_mnabfd_mntd_hyperbola)

anova(mod_lmer, mod_mnabfd_mntd_e, mod_mnabfd_mntd_e2,
      mod_mnabfd_mntd_power, mod_mnabfd_mntd_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)



####### mnnd ~ mnconti_func_d ######
model_e1 = nlsLM(mnnd ~ max_e(a, b, dist = mnconti_func_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnnd ~ max_e_2(b, dist = mnconti_func_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnnd ~ power(a, b, dist = mnconti_func_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnnd ~ hyperbola(a, b, dist = mnconti_func_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnnd ~ logistic(a, b, c,
                                       dist = mnconti_func_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnnd_mnconti_func_dlmer = lmer(mnnd ~ mnconti_func_d + (1|field/f_p),
                                    data = dat_suc_sp, REML = TRUE))
summary(mod_mnnd_mnconti_func_dlmer)

(mod_mnnd_mnconti_func_d_pglmm = pglmm(mnnd ~ mnconti_func_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                       family = "gaussian", cov_ranef = list(species = tree),
                                       bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                   config = TRUE),
                                                            quantiles=c(0.025,0.5,0.975)),
                                       bayes = T))
### nnd correlate with nconti_func_d

#### predictive curve for nnd ~ nconti_func_d
lincombs.data.nnd.nconti_func_d = data.frame(mnconti_func_d=seq(0.0001,
                                                                max(dat_suc_sp$mnconti_func_d),length=100))

lincombs.matrix.nnd.nconti_func_d=model.matrix(~mnconti_func_d,
                                               data=lincombs.data.nnd.nconti_func_d)
lincombs.matrix.nnd.nconti_func_d=as.data.frame(lincombs.matrix.nnd.nconti_func_d)
lincombs.nnd.nconti_func_d=inla.make.lincombs(lincombs.matrix.nnd.nconti_func_d)

inla.model_lincombs.nnd.nconti_func_d = pglmm(mnnd ~  mnconti_func_d+(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "gaussian", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.nnd.nconti_func_d,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)


inla.model_lincombs.nnd.nconti_func_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nnd.nconti_func_d$predicted.value=inla.model_lincombs.nnd.nconti_func_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nnd.nconti_func_d$lower=inla.model_lincombs.nnd.nconti_func_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nnd.nconti_func_d$upper=inla.model_lincombs.nnd.nconti_func_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nnd.nconti_func_d

save(lincombs.data.nnd.nconti_func_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nnd.nconti_func_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnnd_mnconti_func_d_e = nlmer(mnnd ~ max_e_fun(dist = mnconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                   dat_suc_sp,
                                   start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnnd_mnconti_func_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnnd_mnconti_func_d_e2 = nlmer(mnnd ~ max_e2_fun(dist = mnconti_func_d, b) ~ (b|field/f_p) + (b|species),
                                    dat_suc_sp,
                                    start = c(b = nls_coff_e2[1,1])))
summary(mod_mnnd_mnconti_func_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnnd_mnconti_func_d_power = nlmer(mnnd ~ power_fun(dist = mnconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                       dat_suc_sp,
                                       start = c(a = nls_coff_power[1,1],
                                                 b = nls_coff_power[2,1])))
summary(mod_mnnd_mnconti_func_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnnd_mnconti_func_d_hyperbola = nlmer(mnnd ~ hyperbola_fun(dist = mnconti_func_d, a, b) ~ a|field/f_p,
                                           dat_suc_sp,
                                           start = c(a = nls_coff_hyperbola[1,1],
                                                     b = nls_coff_hyperbola[2,1])))
summary(mod_mnnd_mnconti_func_d_hyperbola)

anova(mod_lmer, mod_mnnd_mnconti_func_d_e, mod_mnnd_mnconti_func_d_e2,
      mod_mnnd_mnconti_func_d_power, mod_mnnd_mnconti_func_d_hyperbola, test="Chisq") ## mod_lmer the best!



####### mnfd ~ mnconti_func_d ######
model_e1 = nlsLM(mnfd ~ max_e(a, b, dist = mnconti_func_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnfd ~ max_e_2(b, dist = mnconti_func_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnfd ~ power(a, b, dist = mnconti_func_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnfd ~ hyperbola(a, b, dist = mnconti_func_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnfd ~ logistic(a, b, c,
                                         dist = mnconti_func_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnfd_mnconti_func_dlmer = lmer(mnfd ~ mnconti_func_d + (1|field/f_p),
                                      data = dat_suc_sp, REML = TRUE))
summary(mod_mnfd_mnconti_func_dlmer)
ggpredict(mod_mnfd_mnconti_func_dlmer, terms = 'mnconti_func_d')

(mod_mnfd_mnconti_func_d_pglmm = pglmm(mnfd ~ mnconti_func_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                         family = "gaussian", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = TRUE),
                                                              quantiles=c(0.025,0.5,0.975)),
                                         bayes = T))
### nfd correlate with nconti_func_d

#### predictive curve for nfd ~ nconti_func_d
lincombs.data.nfd.nconti_func_d = data.frame(mnconti_func_d=seq(0.0001,
                                                                  max(dat_suc_sp$mnconti_func_d),
                                                                  length=100))

lincombs.matrix.nfd.nconti_func_d=model.matrix(~mnconti_func_d,
                                                 data=lincombs.data.nfd.nconti_func_d)
lincombs.matrix.nfd.nconti_func_d=as.data.frame(lincombs.matrix.nfd.nconti_func_d)
lincombs.nfd.nconti_func_d=inla.make.lincombs(lincombs.matrix.nfd.nconti_func_d)

inla.model_lincombs.nfd.nconti_func_d = pglmm(mnfd ~  mnconti_func_d+(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_suc_sp,
                                                family = "gaussian", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.nfd.nconti_func_d,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

inla.model_lincombs.nfd.nconti_func_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nfd.nconti_func_d$predicted.value=inla.model_lincombs.nfd.nconti_func_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nfd.nconti_func_d$lower=inla.model_lincombs.nfd.nconti_func_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nfd.nconti_func_d$upper=inla.model_lincombs.nfd.nconti_func_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nfd.nconti_func_d

save(lincombs.data.nfd.nconti_func_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nfd.nconti_func_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnfd_mnconti_func_d_e = nlmer(mnfd ~ max_e_fun(dist = mnconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                     dat_suc_sp,
                                     start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnfd_mnconti_func_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnfd_mnconti_func_d_e2 = nlmer(mnfd ~ max_e2_fun(dist = mnconti_func_d, b) ~ (b|field/f_p) + (b|species),
                                      dat_suc_sp,
                                      start = c(b = nls_coff_e2[1,1])))
summary(mod_mnfd_mnconti_func_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnfd_mnconti_func_d_power = nlmer(mnfd ~ power_fun(dist = mnconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                         dat_suc_sp,
                                         start = c(a = nls_coff_power[1,1],
                                                   b = nls_coff_power[2,1])))
summary(mod_mnfd_mnconti_func_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnfd_mnconti_func_d_hyperbola = nlmer(mnfd ~ hyperbola_fun(dist = mnconti_func_d, a, b) ~ a|field/f_p,
                                             dat_suc_sp,
                                             start = c(a = nls_coff_hyperbola[1,1],
                                                       b = nls_coff_hyperbola[2,1])))
summary(mod_mnfd_mnconti_func_d_hyperbola)

anova(mod_lmer, mod_mnfd_mnconti_func_d_e, mod_mnfd_mnconti_func_d_e2,
      mod_mnfd_mnconti_func_d_power, mod_mnfd_mnconti_func_d_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)



####### mnabfd ~ mnconti_func_d ######
model_e1 = nlsLM(mnabfd ~ max_e(a, b, dist = mnconti_func_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnabfd ~ max_e_2(b, dist = mnconti_func_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnabfd ~ power(a, b, dist = mnconti_func_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnabfd ~ hyperbola(a, b, dist = mnconti_func_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnabfd ~ logistic(a, b, c,
                                           dist = mnconti_func_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnabfd_mnconti_func_dlmer = lmer(mnabfd ~ mnconti_func_d + (1|field/f_p),
                                        data = dat_suc_sp, REML = TRUE))
summary(mod_mnabfd_mnconti_func_dlmer)
ggpredict(mod_mnabfd_mnconti_func_dlmer, terms = 'mnconti_func_d')

(mod_mnabfd_mnconti_func_d_pglmm = pglmm(mnabfd ~ mnconti_func_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                           family = "gaussian", cov_ranef = list(species = tree),
                                           bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                       config = TRUE),
                                                                quantiles=c(0.025,0.5,0.975)),
                                           bayes = T))
### nabfd correlate with nconti_func_d

#### predictive curve for nabfd ~ nconti_func_d
lincombs.data.nabfd.nconti_func_d = data.frame(mnconti_func_d=seq(0.0001,
                                                                    max(dat_suc_sp$mnconti_func_d),length=100))

lincombs.matrix.nabfd.nconti_func_d=model.matrix(~mnconti_func_d,
                                                   data=lincombs.data.nabfd.nconti_func_d)
lincombs.matrix.nabfd.nconti_func_d=as.data.frame(lincombs.matrix.nabfd.nconti_func_d)
lincombs.nabfd.nconti_func_d=inla.make.lincombs(lincombs.matrix.nabfd.nconti_func_d)

inla.model_lincombs.nabfd.nconti_func_d = pglmm(mnabfd ~  mnconti_func_d+(1|species) + 
                                                    (1|f_p) + (1|field), data = dat_suc_sp,
                                                  family = "gaussian", cov_ranef = list(species = tree),
                                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                              config = TRUE),
                                                                       quantiles=c(0.025,0.5,0.975),
                                                                       lincomb=lincombs.nabfd.nconti_func_d,
                                                                       control.predictor=list(compute=T)),
                                                  bayes = T)

inla.model_lincombs.nabfd.nconti_func_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nabfd.nconti_func_d$predicted.value=inla.model_lincombs.nabfd.nconti_func_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nabfd.nconti_func_d$lower=inla.model_lincombs.nabfd.nconti_func_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nabfd.nconti_func_d$upper=inla.model_lincombs.nabfd.nconti_func_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nabfd.nconti_func_d

save(lincombs.data.nabfd.nconti_func_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nabfd.nconti_func_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnabfd_mnconti_func_d_e = nlmer(mnabfd ~ max_e_fun(dist = mnconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                       dat_suc_sp,
                                       start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnabfd_mnconti_func_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnabfd_mnconti_func_d_e2 = nlmer(mnabfd ~ max_e2_fun(dist = mnconti_func_d, b) ~ (b|field/f_p) + (b|species),
                                        dat_suc_sp,
                                        start = c(b = nls_coff_e2[1,1])))
summary(mod_mnabfd_mnconti_func_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnabfd_mnconti_func_d_power = nlmer(mnabfd ~ power_fun(dist = mnconti_func_d, a, b) ~ (a|field/f_p) + (a|species),
                                           dat_suc_sp,
                                           start = c(a = nls_coff_power[1,1],
                                                     b = nls_coff_power[2,1])))
summary(mod_mnabfd_mnconti_func_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnabfd_mnconti_func_d_hyperbola = nlmer(mnabfd ~ hyperbola_fun(dist = mnconti_func_d, a, b) ~ a|field/f_p,
                                               dat_suc_sp,
                                               start = c(a = nls_coff_hyperbola[1,1],
                                                         b = nls_coff_hyperbola[2,1])))
summary(mod_mnabfd_mnconti_func_d_hyperbola)

anova(mod_lmer, mod_mnabfd_mnconti_func_d_e, mod_mnabfd_mnconti_func_d_e2,
      mod_mnabfd_mnconti_func_d_power, mod_mnabfd_mnconti_func_d_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)



####### mnnd ~ mnfunc_d ######
model_e1 = nlsLM(mnnd ~ max_e(a, b, dist = mnfunc_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnnd ~ max_e_2(b, dist = mnfunc_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnnd ~ power(a, b, dist = mnfunc_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnnd ~ hyperbola(a, b, dist = mnfunc_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnnd ~ logistic(a, b, c,
                                       dist = mnfunc_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnnd_mnfunc_dlmer = lmer(mnnd ~ mnfunc_d + (1|field/f_p),
                              data = dat_suc_sp, REML = TRUE))
summary(mod_mnnd_mnfunc_dlmer)

(mod_mnnd_mnfunc_d_pglmm = pglmm(mnnd ~ mnfunc_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                 family = "gaussian", cov_ranef = list(species = tree),
                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                             config = TRUE),
                                                      quantiles=c(0.025,0.5,0.975)),
                                 bayes = T))
### nnd correlate with nfunc_d

#### predictive curve for nnd ~ nfunc_d
lincombs.data.nnd.nfunc_d = data.frame(mnfunc_d=seq(0.0001,
                                                    max(dat_suc_sp$mnfunc_d),length=100))

lincombs.matrix.nnd.nfunc_d=model.matrix(~mnfunc_d,
                                         data=lincombs.data.nnd.nfunc_d)
lincombs.matrix.nnd.nfunc_d=as.data.frame(lincombs.matrix.nnd.nfunc_d)
lincombs.nnd.nfunc_d=inla.make.lincombs(lincombs.matrix.nnd.nfunc_d)

inla.model_lincombs.nnd.nfunc_d = pglmm(mnnd ~  mnfunc_d+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                        family = "gaussian", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.nnd.nfunc_d,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)


inla.model_lincombs.nnd.nfunc_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nnd.nfunc_d$predicted.value=inla.model_lincombs.nnd.nfunc_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nnd.nfunc_d$lower=inla.model_lincombs.nnd.nfunc_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nnd.nfunc_d$upper=inla.model_lincombs.nnd.nfunc_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nnd.nfunc_d

save(lincombs.data.nnd.nfunc_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nnd.nfunc_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnnd_mnfunc_d_e = nlmer(mnnd ~ max_e_fun(dist = mnfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                             dat_suc_sp,
                             start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnnd_mnfunc_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnnd_mnfunc_d_e2 = nlmer(mnnd ~ max_e2_fun(dist = mnfunc_d, b) ~ (b|field/f_p) + (b|species),
                              dat_suc_sp,
                              start = c(b = nls_coff_e2[1,1])))
summary(mod_mnnd_mnfunc_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnnd_mnfunc_d_power = nlmer(mnnd ~ power_fun(dist = mnfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                                 dat_suc_sp,
                                 start = c(a = nls_coff_power[1,1],
                                           b = nls_coff_power[2,1])))
summary(mod_mnnd_mnfunc_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnnd_mnfunc_d_hyperbola = nlmer(mnnd ~ hyperbola_fun(dist = mnfunc_d, a, b) ~ a|field/f_p,
                                     dat_suc_sp,
                                     start = c(a = nls_coff_hyperbola[1,1],
                                               b = nls_coff_hyperbola[2,1])))
summary(mod_mnnd_mnfunc_d_hyperbola)

anova(mod_lmer, mod_mnnd_mnfunc_d_e, mod_mnnd_mnfunc_d_e2,
      mod_mnnd_mnfunc_d_power, mod_mnnd_mnfunc_d_hyperbola, test="Chisq") ## mod_lmer the best!



####### mnfd ~ mnfunc_d ######
model_e1 = nlsLM(mnfd ~ max_e(a, b, dist = mnfunc_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnfd ~ max_e_2(b, dist = mnfunc_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnfd ~ power(a, b, dist = mnfunc_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnfd ~ hyperbola(a, b, dist = mnfunc_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnfd ~ logistic(a, b, c,
                                         dist = mnfunc_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnfd_mnfunc_dlmer = lmer(mnfd ~ mnfunc_d + (1|field/f_p),
                                data = dat_suc_sp, REML = TRUE))
summary(mod_mnfd_mnfunc_dlmer)
ggpredict(mod_mnfd_mnfunc_dlmer, terms = 'mnfunc_d')

(mod_mnfd_mnfunc_d_pglmm = pglmm(mnfd ~ mnfunc_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                   family = "gaussian", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T))
### nabfd correlate with nfunc_d negatively ###

#### predictive curve for nfd ~ nfunc_d
lincombs.data.nfd.nfunc_d = data.frame(mnfunc_d=seq(0.0001,
                                                      max(dat_suc_sp$mnfunc_d),length=100))

lincombs.matrix.nfd.nfunc_d=model.matrix(~mnfunc_d,
                                           data=lincombs.data.nfd.nfunc_d)
lincombs.matrix.nfd.nfunc_d=as.data.frame(lincombs.matrix.nfd.nfunc_d)
lincombs.nfd.nfunc_d=inla.make.lincombs(lincombs.matrix.nfd.nfunc_d)

inla.model_lincombs.nfd.nfunc_d = pglmm(mnfd ~ mnfunc_d+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_suc_sp,
                                          family = "gaussian", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.nfd.nfunc_d,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)

inla.model_lincombs.nfd.nfunc_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nfd.nfunc_d$predicted.value=inla.model_lincombs.nfd.nfunc_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nfd.nfunc_d$lower=inla.model_lincombs.nfd.nfunc_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nfd.nfunc_d$upper=inla.model_lincombs.nfd.nfunc_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nfd.nfunc_d

save(lincombs.data.nfd.nfunc_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nfd.nfunc_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnfd_mnfunc_d_e = nlmer(mnfd ~ max_e_fun(dist = mnfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                               dat_suc_sp,
                               start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnfd_mnfunc_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnfd_mnfunc_d_e2 = nlmer(mnfd ~ max_e2_fun(dist = mnfunc_d, b) ~ (b|field/f_p) + (b|species),
                                dat_suc_sp,
                                start = c(b = nls_coff_e2[1,1])))
summary(mod_mnfd_mnfunc_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnfd_mnfunc_d_power = nlmer(mnfd ~ power_fun(dist = mnfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                                   dat_suc_sp,
                                   start = c(a = nls_coff_power[1,1],
                                             b = nls_coff_power[2,1])))
summary(mod_mnfd_mnfunc_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnfd_mnfunc_d_hyperbola = nlmer(mnfd ~ hyperbola_fun(dist = mnfunc_d, a, b) ~ a|field/f_p,
                                       dat_suc_sp,
                                       start = c(a = nls_coff_hyperbola[1,1],
                                                 b = nls_coff_hyperbola[2,1])))
summary(mod_mnfd_mnfunc_d_hyperbola)

anova(mod_lmer, mod_mnfd_mnfunc_d_e, mod_mnfd_mnfunc_d_e2,
      mod_mnfd_mnfunc_d_power, mod_mnfd_mnfunc_d_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)



####### mnabfd ~ mnfunc_d ######
model_e1 = nlsLM(mnabfd ~ max_e(a, b, dist = mnfunc_d),
                 start=list(a = 1, b = 0.01),
                 data = dat_suc_sp)
model_e2 = nlsLM(mnabfd ~ max_e_2(b, dist = mnfunc_d),
                 start=list(b = 0.01),
                 data = dat_suc_sp)
model_power = nlsLM(mnabfd ~ power(a, b, dist = mnfunc_d),
                    start=list(a = 1, b = 0.01),
                    data = dat_suc_sp)
model_hyperbola = nlsLM(mnabfd ~ hyperbola(a, b, dist = mnfunc_d),
                        start=list(a = 1, b = 0.01),
                        data = dat_suc_sp)
model_logistic = nlsLM(mnabfd ~ logistic(a, b, c,
                                           dist = mnfunc_d),
                       start=list(a = 2, b = 1, c = 5),
                       data = dat_suc_sp)

nls_coff_e1 = summary(model_e1)$coefficients
nls_coff_e2 = summary(model_e2)$coefficients
nls_coff_power = summary(model_power)$coefficients
nls_coff_hyperbola = summary(model_hyperbola)$coefficients

## lm
(mod_mnabfd_mnfunc_dlmer = lmer(mnabfd ~ mnfunc_d + (1|field/f_p),
                                  data = dat_suc_sp, REML = TRUE))
summary(mod_mnabfd_mnfunc_dlmer)
ggpredict(mod_mnabfd_mnfunc_dlmer, terms = 'mnfunc_d')

(mod_mnabfd_mnfunc_d_pglmm = pglmm(mnabfd ~ mnfunc_d + (1|species) + (1|f_p) + (1|field), data = dat_suc_sp,
                                     family = "gaussian", cov_ranef = list(species = tree),
                                     bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                 config = TRUE),
                                                          quantiles=c(0.025,0.5,0.975)),
                                     bayes = T))
### nabfd correlate with nfunc_d positively ###

#### predictive curve for nabfd ~ nfunc_d
lincombs.data.nabfd.nfunc_d = data.frame(mnfunc_d=seq(0.0001,
                                                        max(dat_suc_sp$mnfunc_d),length=100))

lincombs.matrix.nabfd.nfunc_d=model.matrix(~mnfunc_d,
                                             data=lincombs.data.nabfd.nfunc_d)
lincombs.matrix.nabfd.nfunc_d=as.data.frame(lincombs.matrix.nabfd.nfunc_d)
lincombs.nabfd.nfunc_d=inla.make.lincombs(lincombs.matrix.nabfd.nfunc_d)

inla.model_lincombs.nabfd.nfunc_d = pglmm(mnabfd ~  mnfunc_d+(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.nabfd.nfunc_d,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.nabfd.nfunc_d$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.nabfd.nfunc_d$predicted.value=inla.model_lincombs.nabfd.nfunc_d$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.nabfd.nfunc_d$lower=inla.model_lincombs.nabfd.nfunc_d$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.nabfd.nfunc_d$upper=inla.model_lincombs.nabfd.nfunc_d$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.nabfd.nfunc_d

save(lincombs.data.nabfd.nfunc_d, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.nabfd.nfunc_d.rdata')

## model_e
max_eform = ~a*(1-exp(-(b*dist)))
max_e_fun = deriv(max_eform,namevec=c("a","b"),
                  function.arg=c("dist","a","b"))
(mod_mnabfd_mnfunc_d_e = nlmer(mnabfd ~ max_e_fun(dist = mnfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                                 dat_suc_sp,
                                 start = c(a = nls_coff_e1[1,1], b=nls_coff_e1[2,1])))
summary(mod_mnabfd_mnfunc_d_e)

## model_e2
max_e2form = ~(1-exp(-(b*dist)))
max_e2_fun = deriv(max_e2form,namevec=c("b"),
                   function.arg=c("dist", "b"))
(mod_mnabfd_mnfunc_d_e2 = nlmer(mnabfd ~ max_e2_fun(dist = mnfunc_d, b) ~ (b|field/f_p) + (b|species),
                                  dat_suc_sp,
                                  start = c(b = nls_coff_e2[1,1])))
summary(mod_mnabfd_mnfunc_d_e2)

## model_power
power_form = ~ a*(dist^b)
power_fun = deriv(power_form,namevec=c("a", "b"),
                  function.arg=c("dist", "a", "b"))
(mod_mnabfd_mnfunc_d_power = nlmer(mnabfd ~ power_fun(dist = mnfunc_d, a, b) ~ (a|field/f_p) + (a|species),
                                     dat_suc_sp,
                                     start = c(a = nls_coff_power[1,1],
                                               b = nls_coff_power[2,1])))
summary(mod_mnabfd_mnfunc_d_power)

## model_hyperbola
hyperbola_form = ~ dist/(a+b*dist)
hyperbola_fun = deriv(hyperbola_form,namevec=c("a", "b"),
                      function.arg=c("dist", "a", "b"))
(mod_mnabfd_mnfunc_d_hyperbola = nlmer(mnabfd ~ hyperbola_fun(dist = mnfunc_d, a, b) ~ a|field/f_p,
                                         dat_suc_sp,
                                         start = c(a = nls_coff_hyperbola[1,1],
                                                   b = nls_coff_hyperbola[2,1])))
summary(mod_mnabfd_mnfunc_d_hyperbola)

anova(mod_lmer, mod_mnabfd_mnfunc_d_e, mod_mnabfd_mnfunc_d_e2,
      mod_mnabfd_mnfunc_d_power, mod_mnabfd_mnfunc_d_hyperbola, test="Chisq") ## mod_lmer the best!
AIC(mod_lmer, mod_e, mod_e2, mod_power, mod_hyperbola)

### Bayesian inla methods to inculde phylogenetic independence for nlmer
pc_prior = list(prec = list("pc.prec",param=c(0.1,0.01)))
dat_suc_sp$species_1 = factor(dat_suc_sp$species)
sp_name = unique(dat_suc_sp$species)

### Calculate the original functional distance
tree = read.tree('data/original data/phylo_tree332.txt')
sp_md = unique(dat_suc_sp$species)
tree_md = keep.tip(tree, sp_md)
vcv_tree_md = ape::vcv(tree_md, model = "Brownian", corr = FALSE)
vcv_tree_md_sparse = inla.as.sparse(solve(vcv_tree_md))

### lmer ###
formula = mnd.a ~ mpd.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)

mod_mnd.a_mpd.a = inla(formula,
                       control.compute = list(dic=T, waic=T, cpo=T),
                       quantiles=c(0.025,0.5,0.975), data=dat_suc_sp)

summary(mod_mnd.a_mpd.a)
mod_mnd.a_mpd.a$summary.fixed

### nlmer ###
lik_e1 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mnd.a ~ a*(1-exp(-(b*mpd.a)))+ field + f_p
                       + species
                       #+ species_1
)
lik_e2 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mnd.a ~ (1-exp(-(b*mpd.a)))+ field + f_p 
                       + species + species_1)
lik_power = inlabru::like(family = "gaussian", data = dat_suc_sp,
                          formula = mnd.a ~ a*(mpd.a^b)+ field + f_p 
                          + species + species_1)
lik_hyperbola = inlabru::like(family = "gaussian", data = dat_suc_sp,
                              formula = mnd.a ~ mpd.a/(a+b*mpd.a)+ field + f_p 
                              + species + species_1)
cmp = ~ a(1) + b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)#+
#species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

cmp_e2 = ~ b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)+
  species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

fit_e1 = bru(components = cmp, lik_e1, 
             options = list(bru_initial = list(a = nls_coff_e1[1,1],
                                               b = nls_coff_e1[2,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e1)
fit_e1$summary.fixed
plot(mpd.a, max_e(fit_e1$summary.fixed[1,1],fit_e1$summary.fixed[2,1],mpd.a))

fit_e2 = bru(components = cmp, lik_e2, 
             options = list(bru_initial = list(b = nls_coff_e2[1,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e2)
fit_e2$summary.fixed
plot(mpd.a, max_e_2(fit_e2$summary.fixed[1,1],mpd.a))

fit_power = bru(components = cmp, lik_power, 
                options = list(bru_initial = list(a = nls_coff_power[1,1], 
                                                  b = nls_coff_power[2,1]),
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975)))
summary(fit_power)
fit_power$summary.fixed
plot(mpd.a, power(fit_power$summary.fixed[2,1], fit_power$summary.fixed[1,1], mpd.a))

fit_hyperbola = bru(components = cmp, lik_hyperbola, 
                    options = list(bru_initial = list(a = nls_coff_hyperbola[1,1], 
                                                      b = nls_coff_hyperbola[2,1]),
                                   control.compute = list(dic=T, waic=T, cpo=T),
                                   quantiles=c(0.025,0.5,0.975)))
summary(fit_hyperbola)
plot(mpd.a, hyperbola(fit_hyperbola$summary.fixed[2,1], fit_hyperbola$summary.fixed[1,1], mpd.a))

c(mod_mnd.a_mpd.a$waic$waic, fit_e1$waic$waic, fit_e2$waic$waic,
  fit_power$waic$waic, fit_hyperbola$waic$waic) ### lmer is the best after considering phy_independece


####### mnd.a ~ single trait differnces ######
# mgrowth.a
formula = mnd.a ~ mgrowth.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mgrowth.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mgrowth.a)
mod_mnd.a_mgrowth.a$summary.fixed

# mspan.a
formula = mnd.a ~ mspan.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mspan.a = inla(formula,
                         control.compute = list(dic=T, waic=T, cpo=T),
                         quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mspan.a)
mod_mnd.a_mspan.a$summary.fixed

# mpollination.a
formula = mnd.a ~ mpollination.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mpollination.a = inla(formula,
                                control.compute = list(dic=T, waic=T, cpo=T),
                                quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mpollination.a)
mod_mnd.a_mpollination.a$summary.fixed

# mdispersal.a
formula = mnd.a ~ mdispersal.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mdispersal.a = inla(formula,
                              control.compute = list(dic=T, waic=T, cpo=T),
                              quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mdispersal.a)
mod_mnd.a_mdispersal.a$summary.fixed

# mclonality.a
formula = mnd.a ~ mclonality.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mclonality.a = inla(formula,
                              control.compute = list(dic=T, waic=T, cpo=T),
                              quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mclonality.a)
mod_mnd.a_mclonality.a$summary.fixed

# mheight.a
formula = mnd.a ~ mheight.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mheight.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mheight.a)
mod_mnd.a_mheight.a$summary.fixed

# mldmc.a
formula = mnd.a ~ mldmc.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mldmc.a = inla(formula,
                         control.compute = list(dic=T, waic=T, cpo=T),
                         quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mldmc.a)
mod_mnd.a_mldmc.a$summary.fixed

# msla.a
formula = mnd.a ~ msla.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_msla.a = inla(formula,
                        control.compute = list(dic=T, waic=T, cpo=T),
                        quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_msla.a)
mod_mnd.a_msla.a$summary.fixed

# mseedmass.a
formula = mnd.a ~ mseedmass.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mnd.a_mseedmass.a = inla(formula,
                             control.compute = list(dic=T, waic=T, cpo=T),
                             quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mnd.a_mseedmass.a)
mod_mnd.a_mseedmass.a$summary.fixed

dat_trait_mnd.a = data.frame(
  Trait = names(dat_suc_sps)[28:36],
  Value = c(mod_mnd.a_mgrowth.a$summary.fixed$mean[2],
            mod_mnd.a_mspan.a$summary.fixed$mean[2],
            mod_mnd.a_mpollination.a$summary.fixed$mean[2],
            mod_mnd.a_mdispersal.a$summary.fixed$mean[2],
            mod_mnd.a_mclonality.a$summary.fixed$mean[2],
            mod_mnd.a_mheight.a$summary.fixed$mean[2],
            mod_mnd.a_mldmc.a$summary.fixed$mean[2],
            mod_mnd.a_msla.a$summary.fixed$mean[2],
            mod_mnd.a_mseedmass.a$summary.fixed$mean[2]),
  lw = c(mod_mnd.a_mgrowth.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mspan.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mpollination.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mdispersal.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mclonality.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mheight.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mldmc.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_msla.a$summary.fixed$'0.025quant'[2],
         mod_mnd.a_mseedmass.a$summary.fixed$'0.025quant'[2]),
  uw = c(mod_mnd.a_mgrowth.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mspan.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mpollination.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mdispersal.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mclonality.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mheight.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mldmc.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_msla.a$summary.fixed$'0.975quant'[2],
         mod_mnd.a_mseedmass.a$summary.fixed$'0.975quant'[2])
)
dat_trait_mnd.a = arrange(dat_trait_mnd.a, dat_trait_mnd.a$Trait)

(p_mnd.a_traits =  ggplot(dat_trait_mnd.a, aes(x = Trait, y = Value, group = 1)) +
    geom_polygon(aes(y = uw), fill = "grey50", alpha = 0.5) +
    geom_polygon(aes(y = lw), fill = "grey99", alpha = 0.7) +
    geom_polygon(fill = NA, colour = "purple", size = 2) +
    coord_polar() +
    geom_hline(yintercept = c(-0.3, 0.3), colour = "black", size = 2) +
    geom_vline(xintercept = c(0:8), colour = "black", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 2)  +
    scale_x_discrete(labels=c("clonality","dispersal","growth","height",
                              "ldmc","pollination","seedmass","sla",
                              "span")) +
    scale_y_continuous(limits = c(-0.3, 0.3),
                       breaks = seq(-0.3, 0.3, by = 0.3),
                       labels = seq(-0.3, 0.3, by = 0.3))+
    theme_custom() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid  = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(-1, "lines"),
          axis.ticks.margin = unit(1.3,"lines"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.line=element_line(),
          axis.line.x=element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0, size = 20))+
    labs(x = "", y = ""))



####### mfd.a ~ mconti_func_d.a ######
mconti_func_d.a = dat_suc_sp$mconti_func_d.a
plot(mconti_func_d.a, dat_suc_sp$mfd.a)
plot(dat_suc_sp$mgrowth.a, dat_suc_sp$mfd.a)
plot(dat_suc_sp$mspan.a, dat_suc_sp$mfd.a)
plot(dat_suc_sp$mpollination.a, dat_suc_sp$mfd.a)
plot(dat_suc_sp$mdispersal.a, dat_suc_sp$mfd.a)
plot(dat_suc_sp$mclonality.a, dat_suc_sp$mfd.a)
plot(dat_suc_sp$mheight.a, dat_suc_sp$mfd.a)
plot(dat_suc_sp$mldmc.a, dat_suc_sp$mfd.a)
plot(dat_suc_sp$msla.a, dat_suc_sp$mfd.a)

#### normal lmer 
(mod_lmer = lmer(mfd.a ~ mconti_func_d.a + (1|field/f_p) + (1|species),
                 data = dat_suc_sp, REML = TRUE))
summary(mod_lmer)

### Bayesian inla methods to inculde phylogenetic independence
pc_prior = list(prec = list("pc.prec",param=c(0.1,0.01)))
dat_suc_sp$species_1 = factor(dat_suc_sp$species)
sp_name = unique(dat_suc_sp$species)

numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

### Calculate the original functional distance
tree = read.tree('data/original data/phylo_tree332.txt')
sp_md = unique(dat_suc_sp$species)
tree_md = keep.tip(tree, sp_md)
vcv_tree_md = ape::vcv(tree_md, model = "Brownian", corr = FALSE)
vcv_tree_md_sparse = inla.as.sparse(solve(vcv_tree_md))

newdata_for_trait = data.frame(mconti_func_d.a=seq(min(dat_suc_sp$mconti_func_d.a),
                                                   max(dat_suc_sp$mconti_func_d.a),length=100))
trait_matrix=model.matrix(~mconti_func_d.a-1,data=newdata_for_trait)
trait_matrix=as.data.frame(trait_matrix)
trait_matrix_lincombs=inla.make.lincombs(trait_matrix)

### lmer ###
formula = mfd.a ~ mconti_func_d.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_mconti_func_d.a = inla(formula,
                                   control.compute = list(dic=T, waic=T, cpo=T),
                                   quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_mconti_func_d.a)
mod_mfd.a_mconti_func_d.a$summary.fixed

mod_mfd.a_mconti_func_d.a_combs = inla(formula,
                                         control.compute = list(dic=T, waic=T, cpo=T),
                                         quantiles=c(0.025,0.5,0.975), data=dat_suc_sps,
                                         lincomb=phy_matrix_lincombs)

# mgrowth.a
formula = mfd.a ~ mgrowth.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_mgrowth.a = inla(formula,
                             control.compute = list(dic=T, waic=T, cpo=T),
                             quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_mgrowth.a)
mod_mfd.a_mgrowth.a$summary.fixed

# mspan.a
formula = mfd.a ~ mspan.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_mspan.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_mspan.a)
mod_mfd.a_mspan.a$summary.fixed

# mpollination.a
formula = mfd.a ~ mpollination.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_mpollination.a = inla(formula,
                                  control.compute = list(dic=T, waic=T, cpo=T),
                                  quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_mpollination.a)
mod_mfd.a_mpollination.a$summary.fixed

# mdispersal.a
formula = mfd.a ~ mdispersal.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_mdispersal.a = inla(formula,
                                control.compute = list(dic=T, waic=T, cpo=T),
                                quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_mdispersal.a)
mod_mfd.a_mdispersal.a$summary.fixed

# mclonality.a
formula = mfd.a ~ mclonality.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_mclonality.a = inla(formula,
                                control.compute = list(dic=T, waic=T, cpo=T),
                                quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_mclonality.a)
mod_mfd.a_mclonality.a$summary.fixed

# mheight.a
formula = mfd.a ~ mheight.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_mheight.a = inla(formula,
                             control.compute = list(dic=T, waic=T, cpo=T),
                             quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_mheight.a)
mod_mfd.a_mheight.a$summary.fixed

# mldmc.a
formula = mfd.a ~ mldmc.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_mldmc.a = inla(formula,
                           control.compute = list(dic=T, waic=T, cpo=T),
                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_mldmc.a)
mod_mfd.a_mldmc.a$summary.fixed

# msla.a
formula = mfd.a ~ msla.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_msla.a = inla(formula,
                          control.compute = list(dic=T, waic=T, cpo=T),
                          quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_msla.a)
mod_mfd.a_msla.a$summary.fixed

# mseedmass.a
formula = mfd.a ~ mseedmass.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mfd.a_mseedmass.a = inla(formula,
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mfd.a_mseedmass.a)
mod_mfd.a_mseedmass.a$summary.fixed

dat_trait_mfd.a = data.frame(
  Trait = names(dat_suc_sps)[28:36],
  Value = c(mod_mfd.a_mgrowth.a$summary.fixed$mean[2],
            mod_mfd.a_mspan.a$summary.fixed$mean[2],
            mod_mfd.a_mpollination.a$summary.fixed$mean[2],
            mod_mfd.a_mdispersal.a$summary.fixed$mean[2],
            mod_mfd.a_mclonality.a$summary.fixed$mean[2],
            mod_mfd.a_mheight.a$summary.fixed$mean[2],
            mod_mfd.a_mldmc.a$summary.fixed$mean[2],
            mod_mfd.a_msla.a$summary.fixed$mean[2],
            mod_mfd.a_mseedmass.a$summary.fixed$mean[2]),
  lw = c(mod_mfd.a_mgrowth.a$summary.fixed$'0.025quant'[2],
         mod_mfd.a_mspan.a$summary.fixed$'0.025quant'[2],
         mod_mfd.a_mpollination.a$summary.fixed$'0.025quant'[2],
         mod_mfd.a_mdispersal.a$summary.fixed$'0.025quant'[2],
         mod_mfd.a_mclonality.a$summary.fixed$'0.025quant'[2],
         mod_mfd.a_mheight.a$summary.fixed$'0.025quant'[2],
         mod_mfd.a_mldmc.a$summary.fixed$'0.025quant'[2],
         mod_mfd.a_msla.a$summary.fixed$'0.025quant'[2],
         mod_mfd.a_mseedmass.a$summary.fixed$'0.025quant'[2]),
  uw = c(mod_mfd.a_mgrowth.a$summary.fixed$'0.975quant'[2],
         mod_mfd.a_mspan.a$summary.fixed$'0.975quant'[2],
         mod_mfd.a_mpollination.a$summary.fixed$'0.975quant'[2],
         mod_mfd.a_mdispersal.a$summary.fixed$'0.975quant'[2],
         mod_mfd.a_mclonality.a$summary.fixed$'0.975quant'[2],
         mod_mfd.a_mheight.a$summary.fixed$'0.975quant'[2],
         mod_mfd.a_mldmc.a$summary.fixed$'0.975quant'[2],
         mod_mfd.a_msla.a$summary.fixed$'0.975quant'[2],
         mod_mfd.a_mseedmass.a$summary.fixed$'0.975quant'[2])
)
dat_trait_mfd.a = arrange(dat_trait_mfd.a, dat_trait_mfd.a$Trait)

(p_mfd.a_traits =  ggplot(dat_trait_mfd.a, aes(x = Trait, y = Value, group = 1)) +
    geom_polygon(aes(y = uw), fill = "grey50", alpha = 0.5) +
    geom_polygon(aes(y = lw), fill = "grey99", alpha = 0.7) +
    geom_polygon(fill = NA, colour = "purple", size = 2) +
    coord_polar() +
    geom_hline(yintercept = c(-0.3, 0.3), colour = "black", size = 2) +
    geom_vline(xintercept = c(0:8), colour = "black", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 2)  +
    scale_x_discrete(labels=c("clonality","dispersal","growth","height",
                              "ldmc","pollination","seedmass","sla",
                              "span")) +
    scale_y_continuous(limits = c(-0.3, 0.3),
                       breaks = seq(-0.3, 0.3, by = 0.3),
                       labels = seq(-0.3, 0.3, by = 0.3))+
    theme_custom() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid  = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(-1, "lines"),
          axis.ticks.margin = unit(1.3,"lines"),
          axis.text.x = element_text(size = 23),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.line=element_line(),
          axis.line.x=element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0, size = 20))+
    labs(x = "", y = ""))


####### mabfd.a ~ mconti_func_d.a ######
mconti_func_d.a = dat_suc_sp$mconti_func_d.a
plot(mconti_func_d.a, dat_suc_sp$mabfd.a)
plot(dat_suc_sp$mgrowth.a, dat_suc_sp$mabfd.a)
plot(dat_suc_sp$mspan.a, dat_suc_sp$mabfd.a)
plot(dat_suc_sp$mpollination.a, dat_suc_sp$mabfd.a)
plot(dat_suc_sp$mdispersal.a, dat_suc_sp$mabfd.a)
plot(dat_suc_sp$mclonality.a, dat_suc_sp$mabfd.a)
plot(dat_suc_sp$mheight.a, dat_suc_sp$mabfd.a)
plot(dat_suc_sp$mldmc.a, dat_suc_sp$mabfd.a)
plot(dat_suc_sp$msla.a, dat_suc_sp$mabfd.a)

#### normal lmer 
(mod_lmer = lmer(mabfd.a ~ mconti_func_d.a + (1|field/f_p) + (1|species),
                 data = dat_suc_sp, REML = TRUE))
summary(mod_lmer)

### Bayesian inla methods to inculde phylogenetic independence
pc_prior = list(prec = list("pc.prec",param=c(0.1,0.01)))
dat_suc_sp$species_1 = factor(dat_suc_sp$species)
sp_name = unique(dat_suc_sp$species)

numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

### Calculate the original functional distance
tree = read.tree('data/original data/phylo_tree332.txt')
sp_md = unique(dat_suc_sp$species)
tree_md = keep.tip(tree, sp_md)
vcv_tree_md = ape::vcv(tree_md, model = "Brownian", corr = FALSE)
vcv_tree_md_sparse = inla.as.sparse(solve(vcv_tree_md))

newdata_for_trait = data.frame(mconti_func_d.a=seq(min(dat_suc_sp$mconti_func_d.a),
                                                   max(dat_suc_sp$mconti_func_d.a),length=100))
trait_matrix=model.matrix(~mconti_func_d.a-1,data=newdata_for_trait)
trait_matrix=as.data.frame(trait_matrix)
trait_matrix_lincombs=inla.make.lincombs(trait_matrix)

### lmer ###
formula = mabfd.a ~ mconti_func_d.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_mconti_func_d.a = inla(formula,
                                     control.compute = list(dic=T, waic=T, cpo=T),
                                     quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_mconti_func_d.a)
mod_mabfd.a_mconti_func_d.a$summary.fixed

mod_mabfd.a_mconti_func_d.a_combs = inla(formula,
                                           control.compute = list(dic=T, waic=T, cpo=T),
                                           quantiles=c(0.025,0.5,0.975), data=dat_suc_sps,
                                           lincomb=phy_matrix_lincombs)

# mgrowth.a
formula = mabfd.a ~ mgrowth.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_mgrowth.a = inla(formula,
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_mgrowth.a)
mod_mabfd.a_mgrowth.a$summary.fixed

# mspan.a
formula = mabfd.a ~ mspan.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_mspan.a = inla(formula,
                             control.compute = list(dic=T, waic=T, cpo=T),
                             quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_mspan.a)
mod_mabfd.a_mspan.a$summary.fixed

# mpollination.a
formula = mabfd.a ~ mpollination.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_mpollination.a = inla(formula,
                                    control.compute = list(dic=T, waic=T, cpo=T),
                                    quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_mpollination.a)
mod_mabfd.a_mpollination.a$summary.fixed

# mdispersal.a
formula = mabfd.a ~ mdispersal.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_mdispersal.a = inla(formula,
                                  control.compute = list(dic=T, waic=T, cpo=T),
                                  quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_mdispersal.a)
mod_mabfd.a_mdispersal.a$summary.fixed

# mclonality.a
formula = mabfd.a ~ mclonality.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_mclonality.a = inla(formula,
                                  control.compute = list(dic=T, waic=T, cpo=T),
                                  quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_mclonality.a)
mod_mabfd.a_mclonality.a$summary.fixed

# mheight.a
formula = mabfd.a ~ mheight.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_mheight.a = inla(formula,
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_mheight.a)
mod_mabfd.a_mheight.a$summary.fixed

# mldmc.a
formula = mabfd.a ~ mldmc.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_mldmc.a = inla(formula,
                             control.compute = list(dic=T, waic=T, cpo=T),
                             quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_mldmc.a)
mod_mabfd.a_mldmc.a$summary.fixed

# msla.a
formula = mabfd.a ~ msla.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_msla.a = inla(formula,
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_msla.a)
mod_mabfd.a_msla.a$summary.fixed

# mseedmass.a
formula = mabfd.a ~ mseedmass.a +
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_mseedmass.a = inla(formula,
                                 control.compute = list(dic=T, waic=T, cpo=T),
                                 quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_mseedmass.a)
mod_mabfd.a_mseedmass.a$summary.fixed

## all traits
formula = mabfd.a ~ mgrowth.a + mspan.a + mpollination.a + mdispersal.a + 
  mclonality.a + mheight.a + mldmc.a + msla.a + mseedmass.a + 
  f(species, model="iid", hyper = pc_prior)+
  f(field, model="iid", hyper = pc_prior)+
  f(f_p, model="iid", hyper = pc_prior) + 
  f(species_1, model="generic0",
    Cmatrix=vcv_tree_md_sparse, values = sp_md, 
    hyper = pc_prior)
mod_mabfd.a_alltraits = inla(formula,
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975), data=dat_suc_sps)
summary(mod_mabfd.a_alltraits)
mod_mabfd.a_alltraits$summary.fixed

dat_trait_mabfd.a = data.frame(
  Trait = names(dat_suc_sps)[28:36],
  Value = c(mod_mabfd.a_mgrowth.a$summary.fixed$mean[2],
            mod_mabfd.a_mspan.a$summary.fixed$mean[2],
            mod_mabfd.a_mpollination.a$summary.fixed$mean[2],
            mod_mabfd.a_mdispersal.a$summary.fixed$mean[2],
            mod_mabfd.a_mclonality.a$summary.fixed$mean[2],
            mod_mabfd.a_mheight.a$summary.fixed$mean[2],
            mod_mabfd.a_mldmc.a$summary.fixed$mean[2],
            mod_mabfd.a_msla.a$summary.fixed$mean[2],
            mod_mabfd.a_mseedmass.a$summary.fixed$mean[2]),
  lw = c(mod_mabfd.a_mgrowth.a$summary.fixed$'0.025quant'[2],
         mod_mabfd.a_mspan.a$summary.fixed$'0.025quant'[2],
         mod_mabfd.a_mpollination.a$summary.fixed$'0.025quant'[2],
         mod_mabfd.a_mdispersal.a$summary.fixed$'0.025quant'[2],
         mod_mabfd.a_mclonality.a$summary.fixed$'0.025quant'[2],
         mod_mabfd.a_mheight.a$summary.fixed$'0.025quant'[2],
         mod_mabfd.a_mldmc.a$summary.fixed$'0.025quant'[2],
         mod_mabfd.a_msla.a$summary.fixed$'0.025quant'[2],
         mod_mabfd.a_mseedmass.a$summary.fixed$'0.025quant'[2]),
  uw = c(mod_mabfd.a_mgrowth.a$summary.fixed$'0.975quant'[2],
         mod_mabfd.a_mspan.a$summary.fixed$'0.975quant'[2],
         mod_mabfd.a_mpollination.a$summary.fixed$'0.975quant'[2],
         mod_mabfd.a_mdispersal.a$summary.fixed$'0.975quant'[2],
         mod_mabfd.a_mclonality.a$summary.fixed$'0.975quant'[2],
         mod_mabfd.a_mheight.a$summary.fixed$'0.975quant'[2],
         mod_mabfd.a_mldmc.a$summary.fixed$'0.975quant'[2],
         mod_mabfd.a_msla.a$summary.fixed$'0.975quant'[2],
         mod_mabfd.a_mseedmass.a$summary.fixed$'0.975quant'[2])
)
dat_trait_mabfd.a = arrange(dat_trait_mabfd.a, dat_trait_mabfd.a$Trait)
require(ggplot2)
(p_mabfd.a_traits =  ggplot(dat_trait_mabfd.a, aes(x = Trait, y = Value, group = 1)) +
    geom_polygon(aes(y = uw), fill = "grey50", alpha = 0.5) +
    geom_polygon(aes(y = lw), fill = "grey99", alpha = 0.7) +
    geom_polygon(fill = NA, colour = "purple", size = 2) +
    coord_polar() +
    geom_hline(yintercept = c(-0.3, 0.3), colour = "black", size = 2) +
    geom_vline(xintercept = c(0:8), colour = "black", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 2)  +
    scale_x_discrete(labels=c("clonality","dispersal","growth","height",
                              "ldmc","pollination","seedmass","sla",
                              "span")) +
    scale_y_continuous(limits = c(-0.3, 0.3),
                       breaks = seq(-0.3, 0.3, by = 0.3),
                       labels = seq(-0.3, 0.3, by = 0.3))+
    theme_custom() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid  = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(-1, "lines"),
          axis.ticks.margin = unit(1.3,"lines"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.line=element_line(),
          axis.line.x=element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0, size = 20))+
    labs(x = "", y = ""))


## Merge plots mpd.a results for mnd.a and mabfd.a ##
gap = ggplot(NULL)+theme_void()
require(ggpubr)
plo_mnd.a_mabfd.a_traits = ggarrange(p_mnd.a_traits, p_mabfd.a_traits,
                                       nrow = 1, ncol = 2,
                                       widths = c(1, 1), labels = c('a)', 'b)'),
                                       hjust = -3, #vjust = 0.1,
                                       font.label = list(size = 100))

ggsave(plot = plo_mnd.a_mabfd.a_traits,
       'results/figures_ages1_35_top50_aposi/plo_mnd.a_mabfd.a_traits.jpg',
       width = 120, height = 40,dpi = 300, units = 'cm',
       limitsize = F)


### nlmer ###
lik_e1 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mnd.a ~ a*(1-exp(-(b*mconti_func_d.a)))+ field + f_p
                       + species + species_1)
lik_e2 = inlabru::like(family = "gaussian", data = dat_suc_sp,
                       formula = mnd.a ~ (1-exp(-(b*mconti_func_d.a)))+ field + f_p 
                       + species + species_1)
lik_power = inlabru::like(family = "gaussian", data = dat_suc_sp,
                          formula = mnd.a ~ a*(mconti_func_d.a^b)+ field + f_p 
                          + species + species_1)
lik_hyperbola = inlabru::like(family = "gaussian", data = dat_suc_sp,
                              formula = mnd.a ~ mconti_func_d.a/(a+b*mconti_func_d.a)+ field + f_p 
                              + species + species_1)
cmp = ~ a(1) + b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)+
  species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

cmp_e2 = ~ b(1) + 
  field(1, model="iid", hyper = pc_prior)+
  f_p(1, model="iid", hyper = pc_prior) +
  species(1, model="iid", hyper = pc_prior)+
  species_1(1, model="generic0", Cmatrix=vcv_tree_md_sparse, values = sp_md, hyper = pc_prior) 

fit_e1 = bru(components = cmp, lik_e1, 
             options = list(bru_initial = list(a = nls_coff_e1[1,1],
                                               b = nls_coff_e1[2,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e1)
fit_e1$summary.fixed
plot(mconti_func_d.a, max_e(fit_e1$summary.fixed[1,1],fit_e1$summary.fixed[2,1],mconti_func_d.a))

fit_e2 = bru(components = cmp, lik_e2, 
             options = list(bru_initial = list(b = nls_coff_e2[1,1]),
                            control.compute = list(dic=T, waic=T, cpo=T),
                            quantiles=c(0.025,0.5,0.975)))
summary(fit_e2)
fit_e2$summary.fixed
plot(mconti_func_d.a, max_e_2(fit_e2$summary.fixed[1,1],mconti_func_d.a))

fit_power = bru(components = cmp, lik_power, 
                options = list(bru_initial = list(a = nls_coff_power[1,1], 
                                                  b = nls_coff_power[2,1]),
                               control.compute = list(dic=T, waic=T, cpo=T),
                               quantiles=c(0.025,0.5,0.975)))
summary(fit_power)
fit_power$summary.fixed
plot(mconti_func_d.a, power(fit_power$summary.fixed[2,1], fit_power$summary.fixed[1,1], mconti_func_d.a))

fit_hyperbola = bru(components = cmp, lik_hyperbola, 
                    options = list(bru_initial = list(a = nls_coff_hyperbola[1,1], 
                                                      b = nls_coff_hyperbola[2,1]),
                                   control.compute = list(dic=T, waic=T, cpo=T),
                                   quantiles=c(0.025,0.5,0.975)))
summary(fit_hyperbola)
plot(mconti_func_d.a, hyperbola(fit_hyperbola$summary.fixed[2,1],
                                fit_hyperbola$summary.fixed[1,1],
                                mconti_func_d.a))

c(mod_mnd.a_mconti_func_d.a$waic$waic, fit_e1$waic$waic, fit_e2$waic$waic,
  fit_power$waic$waic, fit_hyperbola$waic$waic) ### lmer is the best after considering phy_independece

mod_nd_pd_fd_fixed = mod_nd_pd_fd$summary.fixed
mod_abfd_pd_fd_fixed = mod_abfd_pd_fd$summary.fixed

require(ggplot2)
mnd.a_predicted.phy = mod_mnd.a_mconti_func_d.a_combs$summary.lincomb.derived
mnd.a_newdata.lmm_phy = data.frame(newdata_for_phy,
                                   mean=mnd.a_predicted.phy$mean,
                                   lower=mnd.a_predicted.phy$'0.025quant',
                                   upper=mnd.a_predicted.phy$'0.975quant')

plo_mnd.a_mconti_func_d.a = ggplot() + 
  geom_point(data = dat_suc_sp, aes(x = mconti_func_d.a, y = mnd.a),
             shape = 1, alpha = 1, color = 'grey', size = 16) + 
  
  geom_line(aes(x = seq(0, max(mconti_func_d.a), length.out = 100),
                y = max_e(mod_e_coff[1,1],mod_e_coff[2,1],seq(0, max(mconti_func_d.a),
                                                              length.out = 100))), 
            col = 'blue', linewidth = 4) + 
  geom_line(data=mnd.a_newdata.lmm_phy, aes(x=mconti_func_d.a,y=mean), linetype = "dashed",
            size = 12, color = 'blue')+
  geom_ribbon(data=mnd.a_newdata.lmm_phy,
              aes(x=mconti_func_d.a,ymin=lower,ymax=upper),fill="purple",alpha=0.3)+
  theme_custom() 

#### Fast start for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for fitted species ####
setwd("D:/R projects/BSS")
require(dplyr)
library(phyr)
library(tibble)
library(lme4)
require(ape)
load('code/results_analyzing/analysing_ages1_35_top50_aposi_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

## Check the co-linearity
### only continuous functional trait distance
car::vif(
  glmer(estab ~ mnd + mfd + mpd + mconti_func_d + (1|species) + (1|f_p),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnd.a + mfd.a + mpd.a + mconti_func_d.a + (1|species) + (1|f_p),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnnd + mnfd + mntd + mnconti_func_d + (1|species) + (1|f_p),
        family=binomial,data=dat_suc_sps)
)

### all functional trait distance
car::vif(
  glmer(estab ~ mnd + mfd + mpd + mfunc_d + (1|species) + (1|f_p),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnd.a + mfd.a + mpd.a + mfunc_d.a + (1|species) + (1|f_p),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnnd + mnfd + mntd + mnfunc_d + (1|species) + (1|f_p),
        family=binomial,data=dat_suc_sps)
)
###no co-linearity problem

##### estab #####

pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

### only continuous trait distance 
estab_model_all_md_conti_func_d = pglmm(estab~mnd+mfd+mpd+mconti_func_d+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sps,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975)),
                                        bayes = T)

estab_model_all_md.a_conti_func_d = pglmm(estab~mnd.a+mfd.a+mpd.a+mconti_func_d.a+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_suc_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)

estab_model_all_mnd_conti_func_d = pglmm(estab~mnnd+mnfd+mntd+mnconti_func_d+(1|species) + 
                                           (1|f_p) + (1|field), data = dat_suc_sps,
                                         family = "binomial", 
                                         cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,
                                                                                     waic=T,
                                                                                     cpo=T,
                                                                                     config = TRUE),
                                                              quantiles=c(0.025,0.5,0.975)),
                                         bayes = T)

c(estab_model_all_md_conti_func_d$WAIC,
  estab_model_all_md.a_conti_func_d$WAIC,
  estab_model_all_mnd_conti_func_d$WAIC) 
summary(estab_model_all_md_conti_func_d)
summary(estab_model_all_md.a_conti_func_d)
summary(estab_model_all_mnd_conti_func_d)

### all functional trait distance
estab_model_all_md_func_d = pglmm(estab~mnd+mfd+mpd+mfunc_d+(1|species) + 
                                    (1|f_p) + (1|field), data = dat_suc_sps,
                                  family = "binomial", cov_ranef = list(species = tree),
                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                              config = TRUE),
                                                       quantiles=c(0.025,0.5,0.975)),
                                  bayes = T)

estab_model_all_md.a_func_d = pglmm(estab~mnd.a+mfd.a+mpd.a+mfunc_d.a+(1|species) + 
                                      (1|f_p) + (1|field), data = dat_suc_sps,
                                    family = "binomial", cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975)),
                                    bayes = T)

estab_model_all_mnd_func_d = pglmm(estab~mnnd+mnfd+mntd+mnfunc_d+(1|species) + 
                                     (1|f_p) + (1|field), data = dat_suc_sps,
                                   family = "binomial", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T)

c(estab_model_all_md_func_d$WAIC, estab_model_all_md.a_func_d$WAIC,
  estab_model_all_mnd_func_d$WAIC) 
summary(estab_model_all_md_func_d)
summary(estab_model_all_md.a_func_d)
summary(estab_model_all_mnd_func_d)

# plot 
estab_data.inla.all.md_conti_func_d_intercept1 = estab_model_all_md_conti_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.all.md_conti_func_d_intercept = estab_data.inla.all.md_conti_func_d_intercept1%>%
  mutate(rowname=c("ND","RFD","PD","FD"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

estab_data.inla.all.md.a_conti_func_d_intercept1 = estab_model_all_md.a_conti_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.all.md.a_conti_func_d_intercept = estab_data.inla.all.md.a_conti_func_d_intercept1%>%
  mutate(rowname=c("mnd.ab","mfitness_d.ab","mpd.ab","mconti_func_d.ab"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

estab_data.inla.all.mnd_conti_func_d_intercept1 = estab_model_all_mnd_conti_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.all.mnd_conti_func_d_intercept = estab_data.inla.all.mnd_conti_func_d_intercept1%>%
  mutate(rowname=c("mnnd","mnfitness_d","mnpd","mnconti_func_d"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

estab_data.inla.all.md_func_d_intercept1 = estab_model_all_md_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.all.md_func_d_intercept = estab_data.inla.all.md_func_d_intercept1%>%
  mutate(rowname=c("MND","MRFD","MPD","MFD"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

estab_data.inla.all.md.a_func_d_intercept1 = estab_model_all_md.a_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.all.md.a_func_d_intercept = estab_data.inla.all.md.a_func_d_intercept1%>%
  mutate(rowname=c("MND.ab","MRFD.ab","MPD.ab","MFD.ab"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

estab_data.inla.all.mnd_func_d_intercept1 = estab_model_all_mnd_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

estab_data.inla.all.mnd_func_d_intercept = estab_data.inla.all.mnd_func_d_intercept1%>%
  mutate(rowname=c("MNND","MNRFD","MNTD","MNFD"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

### Effect size plot for mean differences
# point + effect size
(estab_md.varied.intercept.plot =
    ggplot(data=estab_data.inla.all.md_func_d_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('establishment')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD','MPD', 'MRFD', 'MND'))+
    scale_x_continuous(limits=c(-1.5,1))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())
# R2%
(estab_md.varied.intercept.R2.plot = 
    ggplot(data=estab_data.inla.all.md_func_d_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MFD','MPD', 'MRFD', 'MND'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
established_all_md = ggarrange(estab_md.varied.intercept.plot,
                               estab_md.varied.intercept.R2.plot,
                               widths=c(2,1.2))


### Effect size plot for abundance weighted mean differences
# point + effect size
(estab_md.a.varied.intercept.plot =
    ggplot(data=estab_data.inla.all.md.a_func_d_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('establishment')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab','MPD.ab', 'MRFD.ab', 'MND.ab'))+
    scale_x_continuous(limits=c(-1.5,1))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())
# R2%
(estab_md.a.varied.intercept.R2.plot = 
    ggplot(data=estab_data.inla.all.md.a_func_d_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MFD.ab','MPD.ab', 'MRFD.ab', 'MND.ab'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
established_all_md.a = ggarrange(estab_md.a.varied.intercept.plot,
                                 estab_md.a.varied.intercept.R2.plot,
                                 widths=c(2,1.2))

### Effect size plot for mean nearest differences
# point + effect size
(estab_mnd.varied.intercept.plot =
    ggplot(data=estab_data.inla.all.mnd_func_d_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('establishment')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNFD','MNTD', 'MNRFD', 'MNND'))+
    scale_x_continuous(limits=c(-1.5,1))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())
# R2%
(estab_mnd.varied.intercept.R2.plot = 
    ggplot(data=estab_data.inla.all.mnd_func_d_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MNFD','MNTD', 'MNRFD', 'MNND'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
established_all_mnd = ggarrange(estab_mnd.varied.intercept.plot,
                                estab_mnd.varied.intercept.R2.plot,
                                widths=c(2,1.2))

###### predictive curves for mnd
lincombs.data.estab.mnd = data.frame(mnd=seq(min(dat_suc_sp$mnd),max(dat_suc_sp$mnd),length=100),
                                     mfd=mean(dat_suc_sp$mfd),
                                     mpd = mean(dat_suc_sp$mpd),
                                     mconti_func_d = mean(dat_suc_sp$mconti_func_d))

lincombs.matrix.estab.mnd=model.matrix(~mnd+mfd+mpd+mconti_func_d,
                                       data=lincombs.data.estab.mnd)
lincombs.matrix.estab.mnd=as.data.frame(lincombs.matrix.estab.mnd)
lincombs.estab.mnd=inla.make.lincombs(lincombs.matrix.estab.mnd)

inla.model_lincombs.estab.mnd = pglmm(estab ~ mnd+mfd+mpd+mconti_func_d+(1|species) + 
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
#lincombs.data.estab.mnd

save(lincombs.data.estab.mnd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnd.rdata")
(estab.mnd.partial.logistic=ggplot(data=lincombs.data.estab.mnd,
                                   aes(x=mnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1',size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd, y=estab,
                                    color = estab),shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x='Introduced-native niche difference',
         y="Establishment probability")+
    annotate(geom="text",x=c(0.33,0.33),y=c(0.80,0.70),
             label=c("italic(β)['I-N ND'] == -2.25",
                     "'95%CI' == '[-3.10, -1.40]'"),parse=T,size=3.5)
)


# Draw logistic curve
lincombs.data.estab.mfd = data.frame(mfd=seq(min(dat_suc_sp$mfd),max(dat_suc_sp$mfd),length=100),
                                       mnd=mean(dat_suc_sp$mnd),
                                       mpd = mean(dat_suc_sp$mpd),
                                       mconti_func_d = mean(dat_suc_sp$mconti_func_d))

lincombs.matrix.estab.mfd=model.matrix(~mnd+mfd+mpd+mconti_func_d,
                                         data=lincombs.data.estab.mfd)
lincombs.matrix.estab.mfd=as.data.frame(lincombs.matrix.estab.mfd)
lincombs.estab.mfd=inla.make.lincombs(lincombs.matrix.estab.mfd)

inla.model_lincombs.estab.mfd = pglmm(estab ~ mnd+mfd+mpd+mconti_func_d+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.estab.mfd,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)

lincombs.posterior.estab.mfd = inla.model_lincombs.estab.mfd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mfd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mfd$predicted.value=unlist(lapply(lincombs.posterior.estab.mfd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mfd$lower=unlist(lapply(lincombs.posterior.estab.mfd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mfd$upper=unlist(lapply(lincombs.posterior.estab.mfd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mfd
save(lincombs.data.estab.mfd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic prediction curve
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfd.rdata")
(estab.mfd.partial.logistic=ggplot(data=lincombs.data.estab.mfd,aes(x=mfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1',size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfd, y=estab, color = estab),
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Introduced-native fitness difference', y=NULL)+
    annotate(geom="text",x=c(0,0),y=c(0.80,0.70),
             label=c("italic(β)['I-N RFD'] == 2.36", "'95%CI' == '[1.51, 3.21]'"),
             parse=T,size=3.5)
)



##### domin #####
dat_dom_sp = dat_suc_sp %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c('establish', 'dominant'))

## Check the co-linearity
car::vif(
  glmer(domin ~ mnd + mfd + mpd + mconti_func_d + (1|species) + (1|f_p) + (1|field),
        family=binomial,data=dat_dom_sps)
)

car::vif(
  glmer(domin ~ mnd.a + mfd.a + mpd.a + mconti_func_d.a + (1|species) + (1|f_p) + (1|field),
        family=binomial,data=dat_dom_sps)
)

car::vif(
  glmer(domin ~ mnnd + mnfd + mntd + mnconti_func_d + (1|species) + (1|f_p) + (1|field),
        family=binomial,data=dat_dom_sps)
)
summary()

### no problem

### Only continuous functional trait distance 
domin_model_all_md_conti_func_d = pglmm(domin~mnd+mfd+mpd+mconti_func_d+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_dom_sps,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975)),
                                        bayes = T)

domin_model_all_md.a_conti_func_d = pglmm(domin~mnd.a+mfd.a+mpd.a+mconti_func_d.a+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_dom_sps,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975)),
                                          bayes = T)

domin_model_all_mnd_conti_func_d = pglmm(domin~mnnd+mnfd+mntd+mnconti_func_d+(1|species) + 
                                           (1|f_p) + (1|field), data = dat_dom_sps,
                                         family = "binomial", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = TRUE),
                                                              quantiles=c(0.025,0.5,0.975)),
                                         bayes = T)

c(domin_model_all_md_conti_func_d$WAIC, domin_model_all_md.a_conti_func_d$WAIC,
  domin_model_all_mnd_conti_func_d$WAIC) 
summary(domin_model_all_md_conti_func_d)
summary(domin_model_all_md.a_conti_func_d)
summary(domin_model_all_mnd_conti_func_d)

### all functional trait distance
domin_model_all_md_func_d = pglmm(domin~mnd+mfd+mpd+mfunc_d+(1|species) + 
                                    (1|f_p) + (1|field), data = dat_dom_sps,
                                  family = "binomial", cov_ranef = list(species = tree),
                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                              config = TRUE),
                                                       quantiles=c(0.025,0.5,0.975)),
                                  bayes = T)

domin_model_all_md.a_func_d = pglmm(domin~mnd.a+mfd.a+mpd.a+mfunc_d.a+(1|species) + 
                                      (1|f_p) + (1|field), data = dat_dom_sps,
                                    family = "binomial", cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975)),
                                    bayes = T)

domin_model_all_mnd_func_d = pglmm(domin~mnnd+mnfd+mntd+mnfunc_d+(1|species) + 
                                     (1|f_p) + (1|field), data = dat_dom_sps,
                                   family = "binomial", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T)

c(domin_model_all_md_func_d$WAIC, domin_model_all_md.a_func_d$WAIC,
  domin_model_all_mnd_func_d$WAIC) ### 
summary(domin_model_all_md_func_d)
summary(domin_model_all_md.a_func_d)
summary(domin_model_all_mnd_func_d)

# plot 
domin_data.inla.all.md_conti_func_d_intercept1 = domin_model_all_md_conti_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.all.md_conti_func_d_intercept = domin_data.inla.all.md_conti_func_d_intercept1%>%
  mutate(rowname=c("ND","RFD","PD","FD"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

domin_data.inla.all.md.a_conti_func_d_intercept1 = domin_model_all_md.a_conti_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.all.md.a_conti_func_d_intercept = domin_data.inla.all.md.a_conti_func_d_intercept1%>%
  mutate(rowname=c("mnd.ab","mfitness_d.ab","mpd.ab","mconti_func_d.ab"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

domin_data.inla.all.mnd_conti_func_d_intercept1 = domin_model_all_mnd_conti_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.all.mnd_conti_func_d_intercept = domin_data.inla.all.mnd_conti_func_d_intercept1%>%
  mutate(rowname=c("mnnd","mnfitness_d","mnpd","mnconti_func_d"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

domin_data.inla.all.md_func_d_intercept1 = domin_model_all_md_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.all.md_func_d_intercept = domin_data.inla.all.md_func_d_intercept1%>%
  mutate(rowname=c("MND","MRFD","MPD","MFD"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

domin_data.inla.all.md.a_func_d_intercept1 = domin_model_all_md.a_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.all.md.a_func_d_intercept = domin_data.inla.all.md.a_func_d_intercept1%>%
  mutate(rowname=c("MND.ab","MRFD.ab","MPD.ab","MFD.ab"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

domin_data.inla.all.mnd_func_d_intercept1 = domin_model_all_mnd_func_d$inla.model$summary.fixed%>%
  rownames_to_column()%>%.[-1,]%>%rename(lower="0.025quant",median="0.5quant",upper="0.975quant")

domin_data.inla.all.mnd_func_d_intercept = domin_data.inla.all.mnd_func_d_intercept1%>%
  mutate(rowname=c("MNND","MNRFD","MNTD","MNFD"))%>%
  mutate(percent=abs(mean)/sum(abs(mean)))

### Effect size plot for mean differences
# point + effect size
(domin_md.varied.intercept.plot =
    ggplot(data=domin_data.inla.all.md_func_d_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('dominance')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD','MPD', 'MRFD', 'MND'))+
    scale_x_continuous(limits=c(-1.5,1))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())
# R2%
(domin_md.varied.intercept.R2.plot = 
    ggplot(data=domin_data.inla.all.md_func_d_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MFD','MPD', 'MRFD', 'MND'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
dominance_all_md = ggarrange(domin_md.varied.intercept.plot,
                             domin_md.varied.intercept.R2.plot,
                             widths=c(2,1.2))


### Effect size plot for abundance weighted mean differences
# point + effect size
(domin_md.a.varied.intercept.plot =
    ggplot(data=domin_data.inla.all.md.a_func_d_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('dominance')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab','MPD.ab', 'MRFD.ab', 'MND.ab'))+
    scale_x_continuous(limits=c(-1.5,1))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())
# R2%
(domin_md.a.varied.intercept.R2.plot = 
    ggplot(data=domin_data.inla.all.md.a_func_d_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MFD.ab','MPD.ab', 'MRFD.ab', 'MND.ab'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
dominance_all_md.a = ggarrange(domin_md.a.varied.intercept.plot,
                               domin_md.a.varied.intercept.R2.plot,
                               widths=c(2,1.2))

### Effect size plot for mean nearest differences
# point + effect size
(domin_mnd.varied.intercept.plot =
    ggplot(data=domin_data.inla.all.mnd_func_d_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('dominance')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNFD','MNTD', 'MNRFD', 'MNND'))+
    scale_x_continuous(limits=c(-1.5,1.1))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())
# R2%
(domin_mnd.varied.intercept.R2.plot = 
    ggplot(data=domin_data.inla.all.mnd_func_d_intercept,aes(percent,rowname,fill=rowname))+
    geom_bar(stat="identity",width=0.5)+
    geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),
              hjust=-0.1, size = 3)+
    theme_void()+
    scale_y_discrete(limits=c('MNFD','MNTD', 'MNRFD', 'MNND'))+
    scale_fill_viridis_d()+
    theme(plot.margin=unit(c(1.5,0,2.1,-0.3),units="lines"))+
    xlim(0,0.8)+
    guides(fill="none"))

# Merge effect size + R2 
dominance_all_mnd = ggarrange(domin_mnd.varied.intercept.plot,
                              domin_mnd.varied.intercept.R2.plot,
                              widths=c(2,1.2))


###### predictive curves
lincombs.data.domin.nd = data.frame(mnd=seq(min(dat_dom_sp$mnd),max(dat_dom_sp$mnd),length=100),
                                    mfd=mean(dat_dom_sp$mfd),
                                    mpd = mean(dat_dom_sp$mpd),
                                    mconti_func_d = mean(dat_dom_sp$mconti_func_d))

lincombs.matrix.domin.nd=model.matrix(~mnd+mfd+mpd+mconti_func_d,
                                      data=lincombs.data.domin.nd)
lincombs.matrix.domin.nd=as.data.frame(lincombs.matrix.domin.nd)
lincombs.domin.nd=inla.make.lincombs(lincombs.matrix.domin.nd)

inla.model_lincombs.domin.nd = pglmm(domin ~ mnd+mfd+mpd+mconti_func_d+(1|species) + 
                                       (1|f_p) + (1|field), data = dat_dom_sp,
                                     family = "binomial", cov_ranef = list(species = tree),
                                     bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                 config = TRUE),
                                                          quantiles=c(0.025,0.5,0.975),
                                                          lincomb=lincombs.domin.nd,
                                                          control.predictor=list(compute=T)),
                                     bayes = T)

lincombs.posterior.domin.nd = inla.model_lincombs.domin.nd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.nd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.nd$predicted.value=unlist(lapply(lincombs.posterior.domin.nd,
                                                     function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.nd$lower=unlist(lapply(lincombs.posterior.domin.nd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.nd$upper=unlist(lapply(lincombs.posterior.domin.nd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.nd

save(lincombs.data.domin.nd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.nd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve
#load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.nd.rdata")
(domin.nd.partial.logistic=ggplot(data=lincombs.data.domin.nd,
                                  aes(x=mnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1',size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd, y=domin,
                                    color = domin),shape=1,
               alpha = 0.5,
               position=position_jitter(height=0.02))+
    labs(x='Established-native niche difference',
         y="Dominance probability")+
    annotate(geom="text",x=c(0.33,0.33),y=c(0.80,0.70),
             label=c("italic(β)['I-N ND'] == -2.21",
                     "'95%CI' == '[-3.73, -0.68]'"),parse=T,size=3.5)
)


# Draw logistic curve
lincombs.data.domin.fd = data.frame(mfd=seq(min(dat_dom_sp$mfd),max(dat_dom_sp$mfd),length=100),
                                      mnd=mean(dat_dom_sp$mnd),
                                      mpd = mean(dat_dom_sp$mpd),
                                      mconti_func_d = mean(dat_dom_sp$mconti_func_d))

lincombs.matrix.domin.fd=model.matrix(~mnd+mfd+mpd+mconti_func_d,
                                        data=lincombs.data.domin.fd)
lincombs.matrix.domin.fd=as.data.frame(lincombs.matrix.domin.fd)
lincombs.domin.fd=inla.make.lincombs(lincombs.matrix.domin.fd)

inla.model_lincombs.domin.fd = pglmm(domin ~ mnd+mfd+mpd+mconti_func_d+(1|species) + 
                                         (1|f_p) + (1|field), data = dat_dom_sp,
                                       family = "binomial", cov_ranef = list(species = tree),
                                       bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                   config = TRUE),
                                                            quantiles=c(0.025,0.5,0.975),
                                                            lincomb=lincombs.domin.fd,
                                                            control.predictor=list(compute=T)),
                                       bayes = T)

lincombs.posterior.domin.fd = inla.model_lincombs.domin.fd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.fd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.fd$predicted.value=unlist(lapply(lincombs.posterior.domin.fd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.fd$lower=unlist(lapply(lincombs.posterior.domin.fd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.fd$upper=unlist(lapply(lincombs.posterior.domin.fd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.fd
save(lincombs.data.domin.fd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.fd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic prediction curve
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.fd.rdata")
(domin.fd.partial.logistic=ggplot(data=lincombs.data.domin.fd,aes(x=mfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill='grey1',alpha=0.2)+
    geom_line(color='grey1',size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfd, y=domin, color = domin),
               shape=1,
               alpha = 0.5,
               position=position_jitter(height=0.02))+
    labs(x='Established-native fitness difference', y=NULL)+
    annotate(geom="text",x=c(0,0),y=c(0.80,0.70),
             label=c("italic(β)['I-N RFD'] == 2.84", "'95%CI' == '[1.32, 4.36]'"),
             parse=T,size=3.5)
)

#### Fig. S1 estab/domin probability ~ different kind of mean differences 
gap = ggplot(NULL)+theme_void()
library(ggpubr)
estab_domin_md = ggarrange(gap, established_all_md,
                           gap, dominance_all_md,
                           nrow = 4, ncol = 1,
                           labels = c('','a)','','b)'),
                           heights = c(0.1, 1, 0.1, 1))

estab_domin_md.a = ggarrange(gap, established_all_md.a,
                             gap, dominance_all_md.a,
                             nrow = 4, ncol = 1,
                             labels = c('','a)','','b)'),
                             heights = c(0.1, 1, 0.1, 1))

estab_domin_mnd = ggarrange(gap, established_all_mnd,
                            gap, dominance_all_mnd,
                            nrow = 4, ncol = 1,
                            labels = c('','a)','','b)'),
                            heights = c(0.1, 1, 0.1, 1))

ggsave(plot = estab_domin_md,
       'results/figures_ages1_35_top50_aposi/fitted_species_d/estab_domin_md.svg',
       width = 10, height = 15,
       dpi = 300, units = 'cm',
       limitsize = F)

ggsave(plot = estab_domin_md.a,
       'results/figures_ages1_35_top50_aposi/fitted_species_d/estab_domin_md.a.svg',
       width = 10, height = 15,
       dpi = 300, units = 'cm',
       limitsize = F)

ggsave(plot = estab_domin_mnd,
       'results/figures_ages1_35_top50_aposi/fitted_species_d/estab_domin_mnd.svg',
       width = 10, height = 15,
       dpi = 300, units = 'cm',
       limitsize = F)




#### Single regression results for analyzing invasion success probability ~  mnd+mfd+mpd+mfunc_d for all species ####
##### establishment predictive curves for md.ab #####
setwd("D:/R projects/BSS")
library(phyr)
library(tibble)
library(lme4)
require(ape)
library(scales)
library(ggthemes)
load('code/results_analyzing/analysing_ages1_35_top50_aposi_data/dat_suc_sp.rdata')
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
inla.model_lincombs.estab.mnd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(3)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnd.a_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mnd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnd.a_single$lower=unlist(lapply(lincombs.posterior.estab.mnd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnd.a_single$upper=unlist(lapply(lincombs.posterior.estab.mnd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnd.a_single
save(lincombs.data.estab.mnd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnd.a
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnd.a_single.rdata")
(estab.mnd.a_single.logistic=ggplot(data=lincombs.data.estab.mnd.a_single,aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1,
              linetype = 2)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd.a, y=estab),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",x=c(-30,-30),y=c(0.35,0.25),
             label=c("italic(β)['ND'[ab]] == '0.03'",
                     "'95%CI' == '[-0.02, 0.08]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mnd.a       -2.4900    -2.8393    -2.1408



### get establishment predict data for mfd.a
lincombs.data.estab.mfd.a_single = data.frame(mfd.a=seq(min(dat_suc_sp$mfd.a),
                                                            max(dat_suc_sp$mfd.a),
                                                            length=100))

lincombs.matrix.estab.mfd.a_single=model.matrix(~mfd.a,
                                                  data=lincombs.data.estab.mfd.a_single)
lincombs.matrix.estab.mfd.a_single=as.data.frame(lincombs.matrix.estab.mfd.a_single)
lincombs.estab.mfd.a_single=inla.make.lincombs(lincombs.matrix.estab.mfd.a_single)

inla.model_lincombs.estab.mfd.a_single = pglmm(estab ~ mfd.a+#(1|species) + 
                                                   (1|f_p) + (1|field), data = dat_suc_sp,
                                                 family = "binomial", cov_ranef = list(species = tree),
                                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                             config = TRUE),
                                                                      quantiles=c(0.025,0.5,0.975),
                                                                      lincomb=lincombs.estab.mfd.a_single,
                                                                      control.predictor=list(compute=T)),
                                                 bayes = T)

lincombs.posterior.estab.mfd.a_single = inla.model_lincombs.estab.mfd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mfd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mfd.a_single$predicted.value=unlist(lapply(lincombs.posterior.estab.mfd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mfd.a_single$lower=unlist(lapply(lincombs.posterior.estab.mfd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mfd.a_single$upper=unlist(lapply(lincombs.posterior.estab.mfd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mfd.a_single
save(lincombs.data.estab.mfd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfd.a
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfd.a_single.rdata")
(estab.mfd.a_single.logistic=ggplot(data=lincombs.data.estab.mfd.a_single,aes(x=mfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfd.a, y=estab),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y=' ')+
    annotate(geom="text",x=c(-125,-125),y=c(0.35,0.25),
             label=c("italic(β)['RFD'[ab]] == '-0.04'",
                     "'95%CI' == '[-0.06, -0.03]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mfd.a      5.8046     5.2584     6.3508


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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mpd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mpd.a_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mpd.a_single.rdata")
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
    annotate(geom="text",x=c(200,200),y=c(0.90,0.80),
             label=c("italic(β)['MPD'[ab]] == '-0.0175'",
                     "'95%CI' == '[-0.0190, -0.0159]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mpd.a_all   -0.0175    -0.0190    -0.0159


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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfunc_d.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfunc_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfunc_d.a_single.rdata")
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
             label=c("italic(β)['MFD'[ab]] == '-2.7069'",
                     "'95%CI' == '[-3.5855, -1.8281]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mfunc_d.a_all -2.7069    -3.5855    -1.8281


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
inla.model_lincombs.domin.mnd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnd.a_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mnd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnd.a_single$lower=unlist(lapply(lincombs.posterior.domin.mnd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnd.a_single$upper=unlist(lapply(lincombs.posterior.domin.mnd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnd.a_single
save(lincombs.data.domin.mnd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnd.a
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnd.a_single.rdata")
(domin.mnd.a_single.logistic=ggplot(data=lincombs.data.domin.mnd.a_single,aes(x=mnd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1],alpha=0.2)+
    geom_line(color=colors_4d[1],size=1,
              linetype = 2)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd.a, y=domin),
               color = colors_4d[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Niche difference (ND)', y='Dominance probability')+
    annotate(geom="text",x=c(-30,-30),y=c(0.35,0.25),
             label=c("italic(β)['ND'[ab]] == '-0.02'",
                     "'95%CI' == '[-0.06, 0.02]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mnd.a        0.0024    -0.0367     0.0415



### get dominance predict data for mfd.a
lincombs.data.domin.mfd.a_single = data.frame(mfd.a=seq(min(dat_dom_sp$mfd.a),
                                                            max(dat_dom_sp$mfd.a),
                                                            length=100))

lincombs.matrix.domin.mfd.a_single=model.matrix(~mfd.a,
                                                  data=lincombs.data.domin.mfd.a_single)
lincombs.matrix.domin.mfd.a_single=as.data.frame(lincombs.matrix.domin.mfd.a_single)
lincombs.domin.mfd.a_single=inla.make.lincombs(lincombs.matrix.domin.mfd.a_single)

inla.model_lincombs.domin.mfd.a_single = pglmm(domin ~ mfd.a+#(1|species) + 
                                                   (1|f_p) + (1|field), data = dat_dom_sp,
                                                 family = "binomial", cov_ranef = list(species = tree),
                                                 bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                             config = TRUE),
                                                                      quantiles=c(0.025,0.5,0.975),
                                                                      lincomb=lincombs.domin.mfd.a_single,
                                                                      control.predictor=list(compute=T)),
                                                 bayes = T)

lincombs.posterior.domin.mfd.a_single = inla.model_lincombs.domin.mfd.a_single$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mfd.a_single$inla.model$summary.fixed[c(1,3,5)]%>%round(4)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mfd.a_single$predicted.value=unlist(lapply(lincombs.posterior.domin.mfd.a_single,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mfd.a_single$lower=unlist(lapply(lincombs.posterior.domin.mfd.a_single,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mfd.a_single$upper=unlist(lapply(lincombs.posterior.domin.mfd.a_single,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mfd.a_single
save(lincombs.data.domin.mfd.a_single, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfd.a
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfd.a_single.rdata")
(domin.mfd.a_single.logistic=ggplot(data=lincombs.data.domin.mfd.a_single,aes(x=mfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1,
              linetype = 1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfd.a, y=domin),
               color = colors_4d[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='Relative fitness difference (RFD)', y='')+
    annotate(geom="text",x=c(-125,-125),y=c(0.35,0.25),
             label=c("italic(β)['RFD'[ab]] == '-0.06'",
                     "'95%CI' == '[-0.09, -0.02]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mfd.a       -0.0427    -0.0591    -0.0263



### get dominance predict data for mpd.a_all
lincombs.data.domin.mpd.a_single = data.frame(mpd.a_all=seq(min(dat_dom_sp$mpd.a_all),
                                                            max(dat_dom_sp$mpd.a_all),
                                                            length=100))

lincombs.matrix.domin.mpd.a_single=model.matrix(~mpd.a_all,
                                                data=lincombs.data.domin.mpd.a_single)
lincombs.matrix.domin.mpd.a_single=as.data.frame(lincombs.matrix.domin.mpd.a_single)
lincombs.domin.mpd.a_single=inla.make.lincombs(lincombs.matrix.domin.mpd.a_single)

inla.model_lincombs.domin.mpd.a_single = pglmm(domin ~ mpd.a_all+(1|species) + 
                                                 (1|f_p) + (1|field), data = dat_dom_sp,
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mpd.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mpd.a_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mpd.a_single.rdata")
(domin.mpd.a_single.logistic=ggplot(data=lincombs.data.domin.mpd.a_single,aes(x=mpd.a_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[3],alpha=0.2)+
    geom_line(color=colors_4d[3],size=1,
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
    annotate(geom="text",x=c(200,200),y=c(0.90,0.80),
             label=c("italic(β)['MPD'[ab]] == '-0.0030'",
                     "'95%CI' == '[-0.0101, 0.0042]'"),
             parse=T,size=3.5)+
    theme_regular()
)
#mpd.a_all   -0.0030    -0.0101     0.0042

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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfunc_d.a_single.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfunc_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfunc_d.a_single.rdata")
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
             label=c("italic(β)['MFD'[ab]] == '5.8854'",
                     "'95%CI' == '[4.1421, 7.6156]'"),
             parse=T,size=3.5)+
    theme_regular()
)
# mfunc_d.a_all  5.8854     4.1421     7.6156


#### Fast start for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for all species ####
setwd("D:/R projects/BSS")
library(phyr)
library(tibble)
library(lme4)
require(ape)
library(scales)
library(ggthemes)
load('code/results_analyzing/analysing_ages1_35_top50_aposi_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

## Check the co-linearity
### only continuous functional trait distance
car::vif(
  glmer(estab ~ mnd + mfd + mpd_all + mconti_func_d_all + (1|species) + (1|f_p) + (1|field),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnd.a + mfd.a + mpd.a_all + mconti_func_d.a_all + (1|species) + (1|f_p)+ (1|field),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnnd + mnfd + mntd_all + mnconti_func_d_all + (1|species) + (1|f_p)+ (1|field),
        family=binomial,data=dat_suc_sps)
)

### all functional trait distance
car::vif(
  glmer(estab ~ mnd + mfd + mpd_all + mfunc_d_all + (1|species) + (1|f_p)+ (1|field),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnd.a + mfd.a + mpd.a_all + mfunc_d.a_all + (1|species) + (1|f_p)+ (1|field),
        family = binomial, data = dat_suc_sps)
)

car::vif(
  glmer(estab ~ mnnd + mnfd + mntd_all + mnfunc_d_all + (1|species) + (1|f_p)+ (1|field),
        family=binomial,data=dat_suc_sps)
)
###no co-linearity problem, all VIF < 3

##### Estab #####
pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))

### all functional trait distance
estab_model_md_func_d_all = pglmm(estab~mnd+mfd+mpd_all+mfunc_d_all+#(1|species) + 
                                    (1|f_p) + (1|field), data = dat_suc_sps,
                                  family = "binomial", cov_ranef = list(species = tree),
                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                              config = TRUE),
                                                       quantiles=c(0.025,0.5,0.975)),
                                  bayes = T)

estab_model_md.a_func_d_all = pglmm(estab~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+#(1|species) + 
                                      (1|f_p) + (1|field), data = dat_suc_sps,
                                    family = "binomial", cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975)),
                                    bayes = T)

estab_model_mnd_func_d_all = pglmm(estab~mnnd+mnfd+mntd_all+mnfunc_d_all+#(1|species) + 
                                     (1|f_p) + (1|field), data = dat_suc_sps,
                                   family = "binomial", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T)

summary(estab_model_md_func_d_all)
summary(estab_model_md.a_func_d_all)
summary(estab_model_mnd_func_d_all)

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


### Effect size plot for mean differences
# point + effect size
require(ggplot2)
require(ggpubr)
(estab_nd_rfd.all.varied.intercept.plot =
    ggplot((data=estab_data.inla.md_func_d.all_intercept %>% 
              filter(rowname %in% c('MND','MRFD'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 2)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD', 'MND'))+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND" = "MND",
                                  "MRFD" = "MRFD"),
                       name = ' ')+
    annotate(geom="text", x=-1.4, y=2.5, family = 'Arial',
             label='Establishment')+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

(estab_pd_fd.all.varied.intercept.plot =
    ggplot(data=(data=estab_data.inla.md_func_d.all_intercept %>% 
                   filter(rowname %in% c('MFD_all','MPD_all'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 2)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD_all','MPD_all'),
                     label=c('MFD','MPD'))+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD_all" = "MPD_all",
                                  "MFD_all" = "MFD_all"),
                       name = ' ')+
    annotate(geom="text", x=-1.4, y=2.5, family = 'Arial',
             label='Establishment')+
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
(estab_md.a.all.varied.intercept.plot =
    ggplot(data=estab_data.inla.md.a_func_d.all_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all', 'MRFD.ab', 'MND.ab'))+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_viridis_d()+
    annotate(geom="text", x=-1.4, y=2.5, family = 'Arial',
             label='Establishment', fontface = 'bold')+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())

(estab_nd_rfd.a.all.varied.intercept.plot =
    ggplot((data=estab_data.inla.md.a_func_d.all_intercept %>% 
              filter(rowname %in% c('MRFD.ab', 'MND.ab'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 3)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD.ab', 'MND.ab'),
                     label = c(expression(RFD['ab']), expression(ND['ab'])))+
    scale_x_continuous(limits=c(-0.5,0.58))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND.ab" = "MND.ab",
                                  "MRFD.ab" = "MRFD.ab"),
                       name = ' ')+
    annotate(geom="text", x=0.31, y=2.5, family = 'Arial',
             label='Establishment', fontface = 'bold')+
    labs(x = '  ', y = '  ')+
    guides(color="none")+
    theme_regular())

(estab_pd_fd.a.all.varied.intercept.plot =
    ggplot(data=(data=estab_data.inla.md.a_func_d.all_intercept %>% 
                   filter(rowname %in% c('MFD.ab_all','MPD.ab_all'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 3)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all'),
                     label = c(expression(MFD['ab']), expression(MPD['ab'])))+
    scale_x_continuous(limits=c(-1.05,1.05))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD.ab_all" = "MPD.ab_all",
                                  "MFD.ab_all" = "MFD.ab_all"),
                       name = ' ')+
    annotate(geom="text", x=0.6, y=2.5, family = 'Arial',
             label='Establishment')+
    labs(x = '  ', y = ' ')+
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
(estab_mnd.all.varied.intercept.plot =
    ggplot(data=estab_data.inla.mnd_func_d.all_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('Establishment')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNFD_all','MNTD_all', 'MNRFD', 'MNND'))+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())

(estab_nnd_nrfd.all.varied.intercept.plot =
    ggplot((data=estab_data.inla.mnd_func_d.all_intercept %>% 
              filter(rowname %in% c('MNRFD', 'MNND'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 2)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNRFD', 'MNND'))+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MNND" = "MNND",
                                  "MNRFD" = "MNRFD"),
                       name = ' ')+
    annotate(geom="text", x=-1.4, y=2.5, family = 'Arial',
             label='Establishment')+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

(estab_npd_nfd.all.varied.intercept.plot =
    ggplot(data=(data=estab_data.inla.mnd_func_d.all_intercept %>% 
                   filter(rowname %in% c('MNFD_all','MNTD_all'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 2)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNFD_all','MNTD_all'),
                     label = c('MNFD','MNTD'))+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MNTD_all" = "MNTD_all",
                                  "MNFD_all" = "MNFD_all"),
                       name = ' ')+
    annotate(geom="text", x=-1.4, y=2.5, family = 'Arial',
             label='Establishment')+
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
                                     mfd=mean(dat_suc_sp$mfd),
                                     mpd_all = mean(dat_suc_sp$mpd_all),
                                     mfunc_d_all = mean(dat_suc_sp$mfunc_d_all))

lincombs.matrix.estab.mnd=model.matrix(~mnd+mfd+mpd_all+mfunc_d_all,
                                       data=lincombs.data.estab.mnd)
lincombs.matrix.estab.mnd=as.data.frame(lincombs.matrix.estab.mnd)
lincombs.estab.mnd=inla.make.lincombs(lincombs.matrix.estab.mnd)

inla.model_lincombs.estab.mnd = pglmm(estab ~ mnd+mfd+mpd_all+mfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment~mnd
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnd.rdata")

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


### get establishment predict data for mfd
lincombs.data.estab.mfd = data.frame(mfd=seq(min(dat_suc_sp$mfd),max(dat_suc_sp$mfd),length=100),
                                       mnd=mean(dat_suc_sp$mnd),
                                       mpd_all = mean(dat_suc_sp$mpd_all),
                                       mfunc_d_all = mean(dat_suc_sp$mfunc_d_all))

lincombs.matrix.estab.mfd=model.matrix(~mnd+mfd+mpd_all+mfunc_d_all,
                                         data=lincombs.data.estab.mfd)
lincombs.matrix.estab.mfd=as.data.frame(lincombs.matrix.estab.mfd)
lincombs.estab.mfd=inla.make.lincombs(lincombs.matrix.estab.mfd)

inla.model_lincombs.estab.mfd = pglmm(estab ~ mnd+mfd+mpd_all+mfunc_d_all+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.estab.mfd,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)

lincombs.posterior.estab.mfd = inla.model_lincombs.estab.mfd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mfd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mfd$predicted.value=unlist(lapply(lincombs.posterior.estab.mfd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mfd$lower=unlist(lapply(lincombs.posterior.estab.mfd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mfd$upper=unlist(lapply(lincombs.posterior.estab.mfd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mfd
save(lincombs.data.estab.mfd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfd
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfd.rdata")
(estab.mfd.partial.logistic=ggplot(data=lincombs.data.estab.mfd,aes(x=mfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfd, y=estab),
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
                                         mfd = mean(dat_suc_sp$mfd),
                                         mfunc_d_all = mean(dat_suc_sp$mfunc_d_all))

lincombs.matrix.estab.mpd_all=model.matrix(~mnd+mfd+mpd_all+mfunc_d_all,
                                           data=lincombs.data.estab.mpd_all)
lincombs.matrix.estab.mpd_all=as.data.frame(lincombs.matrix.estab.mpd_all)
lincombs.estab.mpd_all=inla.make.lincombs(lincombs.matrix.estab.mpd_all)

inla.model_lincombs.estab.mpd_all = pglmm(estab ~ mnd+mfd+mpd_all+mfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mpd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mpd_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mpd_all.rdata")
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
                                             mfd = mean(dat_suc_sp$mfd),
                                             mfunc_d_all = seq(min(dat_suc_sp$mfunc_d_all),
                                                               max(dat_suc_sp$mfunc_d_all),length=100))

lincombs.matrix.estab.mfunc_d_all=model.matrix(~mnd+mfd+mpd_all+mfunc_d_all,
                                               data=lincombs.data.estab.mfunc_d_all)
lincombs.matrix.estab.mfunc_d_all=as.data.frame(lincombs.matrix.estab.mfunc_d_all)
lincombs.estab.mfunc_d_all=inla.make.lincombs(lincombs.matrix.estab.mfunc_d_all)

inla.model_lincombs.estab.mfunc_d_all = pglmm(estab ~ mnd+mfd+mpd_all+mfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfunc_d_all.rdata")
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
                                      mnfd=mean(dat_suc_sp$mnfd),
                                      mntd_all = mean(dat_suc_sp$mntd_all),
                                      mnfunc_d_all = mean(dat_suc_sp$mnfunc_d_all))

lincombs.matrix.estab.mnnd=model.matrix(~mnnd+mnfd+mntd_all+mnfunc_d_all,
                                        data=lincombs.data.estab.mnnd)
lincombs.matrix.estab.mnnd=as.data.frame(lincombs.matrix.estab.mnnd)
lincombs.estab.mnnd=inla.make.lincombs(lincombs.matrix.estab.mnnd)

inla.model_lincombs.estab.mnnd = pglmm(estab ~ mnnd+mnfd+mntd_all+mnfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment~mnnd
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnnd.rdata")

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


### get establishment predict data for mnfd
lincombs.data.estab.mnfd = data.frame(mnfd=seq(min(dat_suc_sp$mnfd),max(dat_suc_sp$mnfd),length=100),
                                        mnnd=mean(dat_suc_sp$mnnd),
                                        mntd_all = mean(dat_suc_sp$mntd_all),
                                        mnfunc_d_all = mean(dat_suc_sp$mnfunc_d_all))

lincombs.matrix.estab.mnfd=model.matrix(~mnnd+mnfd+mntd_all+mnfunc_d_all,
                                          data=lincombs.data.estab.mnfd)
lincombs.matrix.estab.mnfd=as.data.frame(lincombs.matrix.estab.mnfd)
lincombs.estab.mnfd=inla.make.lincombs(lincombs.matrix.estab.mnfd)

inla.model_lincombs.estab.mnfd = pglmm(estab ~ mnnd+mnfd+mntd_all+mnfunc_d_all+(1|species) + 
                                           (1|f_p) + (1|field), data = dat_suc_sp,
                                         family = "binomial", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = TRUE),
                                                              quantiles=c(0.025,0.5,0.975),
                                                              lincomb=lincombs.estab.mnfd,
                                                              control.predictor=list(compute=T)),
                                         bayes = T)

lincombs.posterior.estab.mnfd = inla.model_lincombs.estab.mnfd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mnfd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnfd$predicted.value=unlist(lapply(lincombs.posterior.estab.mnfd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnfd$lower=unlist(lapply(lincombs.posterior.estab.mnfd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnfd$upper=unlist(lapply(lincombs.posterior.estab.mnfd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mnfd
save(lincombs.data.estab.mnfd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnfd
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnfd.rdata")
(estab.mnfd.partial.logistic=ggplot(data=lincombs.data.estab.mnfd,aes(x=mnfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnfd, y=estab),
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
                                          mnfd = mean(dat_suc_sp$mnfd),
                                          mnfunc_d_all = mean(dat_suc_sp$mnfunc_d_all))

lincombs.matrix.estab.mntd_all=model.matrix(~mnnd+mnfd+mntd_all+mnfunc_d_all,
                                            data=lincombs.data.estab.mntd_all)
lincombs.matrix.estab.mntd_all=as.data.frame(lincombs.matrix.estab.mntd_all)
lincombs.estab.mntd_all=inla.make.lincombs(lincombs.matrix.estab.mntd_all)

inla.model_lincombs.estab.mntd_all = pglmm(estab ~ mnnd+mnfd+mntd_all+mnfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mntd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mntd_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mntd_all.rdata")
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
                                              mnfd = mean(dat_suc_sp$mnfd),
                                              mnfunc_d_all = seq(min(dat_suc_sp$mnfunc_d_all),
                                                                 max(dat_suc_sp$mnfunc_d_all),length=100))

lincombs.matrix.estab.mnfunc_d_all=model.matrix(~mnnd+mnfd+mntd_all+mnfunc_d_all,
                                                data=lincombs.data.estab.mnfunc_d_all)
lincombs.matrix.estab.mnfunc_d_all=as.data.frame(lincombs.matrix.estab.mnfunc_d_all)
lincombs.estab.mnfunc_d_all=inla.make.lincombs(lincombs.matrix.estab.mnfunc_d_all)

inla.model_lincombs.estab.mnfunc_d_all = pglmm(estab ~ mnnd+mnfd+mntd_all+mnfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnfunc_d_all.rdata")
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
                                       mfd.a=mean(dat_suc_sp$mfd.a),
                                       mpd.a_all = mean(dat_suc_sp$mpd.a_all),
                                       mfunc_d.a_all = mean(dat_suc_sp$mfunc_d.a_all))

lincombs.matrix.estab.mnd.a=model.matrix(~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all,
                                         data=lincombs.data.estab.mnd.a)
lincombs.matrix.estab.mnd.a=as.data.frame(lincombs.matrix.estab.mnd.a)
lincombs.estab.mnd.a=inla.make.lincombs(lincombs.matrix.estab.mnd.a)

inla.model_lincombs.estab.mnd.a = pglmm(estab ~ mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnd.a.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment~mnd.ab
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mnd.a.rdata")
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


### get establishment predict data for mfd.ab
lincombs.data.estab.mfd.a = data.frame(mfd.a=seq(min(dat_suc_sp$mfd.a),max(dat_suc_sp$mfd.a),length=100),
                                         mnd.a=mean(dat_suc_sp$mnd.a),
                                         mpd.a_all = mean(dat_suc_sp$mpd.a_all),
                                         mfunc_d.a_all = mean(dat_suc_sp$mfunc_d.a_all))

lincombs.matrix.estab.mfd.a=model.matrix(~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all,
                                           data=lincombs.data.estab.mfd.a)
lincombs.matrix.estab.mfd.a=as.data.frame(lincombs.matrix.estab.mfd.a)
lincombs.estab.mfd.a=inla.make.lincombs(lincombs.matrix.estab.mfd.a)

inla.model_lincombs.estab.mfd.a = pglmm(estab ~ mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_suc_sp,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.estab.mfd.a,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)

lincombs.posterior.estab.mfd.a = inla.model_lincombs.estab.mfd.a$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.estab.mfd.a$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mfd.a$predicted.value=unlist(lapply(lincombs.posterior.estab.mfd.a,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mfd.a$lower=unlist(lapply(lincombs.posterior.estab.mfd.a,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mfd.a$upper=unlist(lapply(lincombs.posterior.estab.mfd.a,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.estab.mfd.a
save(lincombs.data.estab.mfd.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfd.a.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfd.ab
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfd.a.rdata")
(estab.mfd.a.partial.logistic=ggplot(data=lincombs.data.estab.mfd.a,aes(x=mfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfd.a, y=estab),
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
                                           mfd.a = mean(dat_suc_sp$mfd.a),
                                           mfunc_d.a_all = mean(dat_suc_sp$mfunc_d.a_all))

lincombs.matrix.estab.mpd.a_all=model.matrix(~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all,
                                             data=lincombs.data.estab.mpd.a_all)
lincombs.matrix.estab.mpd.a_all=as.data.frame(lincombs.matrix.estab.mpd.a_all)
lincombs.estab.mpd.a_all=inla.make.lincombs(lincombs.matrix.estab.mpd.a_all)

inla.model_lincombs.estab.mpd.a_all = pglmm(estab ~ mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mpd.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mpd.a_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mpd.a_all.rdata")
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
                                               mfd.a = mean(dat_suc_sp$mfd.a),
                                               mfunc_d.a_all = seq(min(dat_suc_sp$mfunc_d.a_all),
                                                                   max(dat_suc_sp$mfunc_d.a_all),length=100))

lincombs.matrix.estab.mfunc_d.a_all=model.matrix(~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all,
                                                 data=lincombs.data.estab.mfunc_d.a_all)
lincombs.matrix.estab.mfunc_d.a_all=as.data.frame(lincombs.matrix.estab.mfunc_d.a_all)
lincombs.estab.mfunc_d.a_all=inla.make.lincombs(lincombs.matrix.estab.mfunc_d.a_all)

inla.model_lincombs.estab.mfunc_d.a_all = pglmm(estab ~ mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfunc_d.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mfunc_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.estab.mfunc_d.a_all.rdata")
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
  glmer(domin ~ mnd + mfd + mpd_all + mfunc_d_all + (1|species) + (1|f_p) + (1|field),
        family=binomial,data=dat_dom_sps)
)

car::vif(
  glmer(domin ~ mnd.a + mfd.a + mpd.a_all + mfunc_d.a_all + (1|species) + (1|f_p)+ (1|field),
        family=binomial,data=dat_dom_sps)
)

car::vif(
  glmer(domin ~ mnnd + mnfd + mntd_all + mnfunc_d_all + (1|species) + (1|f_p)+ (1|field),
        family=binomial,data=dat_dom_sps)
)
summary(glmer(domin ~ mnd + mfd + mpd_all + mfunc_d_all + 
                (1|species) + (1|f_p) + (1|field),
              family=binomial,data=dat_dom_sps))
summary(glmer(domin ~ mpd_all + mfunc_d_all + 
                (1|species) + (1|f_p) + (1|field),
              family=binomial,data=dat_dom_sps))
### no colinearity problem

### all functional trait distance
domin_model_md_func_d_all = pglmm(domin~mnd+mfd+mpd_all+mfunc_d_all+#(1|species) + 
                                    (1|f_p) + (1|field), data = dat_dom_sps,
                                  family = "binomial", cov_ranef = list(species = tree),
                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                              config = TRUE),
                                                       quantiles=c(0.025,0.5,0.975)),
                                  bayes = T)

domin_model_md.a_func_d_all = pglmm(domin~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+#(1|species) + 
                                      (1|f_p) + (1|field), data = dat_dom_sps,
                                    family = "binomial", cov_ranef = list(species = tree),
                                    bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                config = TRUE),
                                                         quantiles=c(0.025,0.5,0.975)),
                                    bayes = T)

domin_model_mnd_func_d_all = pglmm(domin~mnnd+mnfd+mntd_all+mnfunc_d_all+#(1|species) + 
                                     (1|f_p) + (1|field), data = dat_dom_sps,
                                   family = "binomial", cov_ranef = list(species = tree),
                                   bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                               config = TRUE),
                                                        quantiles=c(0.025,0.5,0.975)),
                                   bayes = T)
summary(domin_model_md_func_d_all)
summary(domin_model_md.a_func_d_all)
summary(domin_model_mnd_func_d_all)

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


### Effect size plot for mean differences
# point + effect size
(domin_md.all.varied.intercept.plot =
    ggplot(data=domin_data.inla.md_func_d.all_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('dominance')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD_all','MPD_all', 'MRFD', 'MND'))+
    scale_x_continuous(limits=c(-1.5,1.5), breaks = seq(-1.5, 1.5, 0.5))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())

(domin_nd_rfd.all.varied.intercept.plot =
    ggplot((data=domin_data.inla.md_func_d.all_intercept %>% 
              filter(rowname %in% c('MND','MRFD'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 2)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD', 'MND'))+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND" = "MND",
                                  "MRFD" = "MRFD"),
                       name = ' ')+
    annotate(geom="text", x=-1.4, y=2.5, family = 'Arial',
             label='Dominance')+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_regular())

(domin_pd_fd.all.varied.intercept.plot =
    ggplot(data=(data=domin_data.inla.md_func_d.all_intercept %>% 
                   filter(rowname %in% c('MFD_all','MPD_all'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 2)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD_all','MPD_all'),
                     label =c('MFD','MPD') )+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD_all" = "MPD_all",
                                  "MFD_all" = "MFD_all"),
                       name = ' ')+
    annotate(geom="text", x=-1.4, y=2.5, family = 'Arial',
             label='Dominance')+
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
(domin_md.a.all.varied.intercept.plot =
    ggplot(data=domin_data.inla.md.a_func_d.all_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('dominance')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all', 'MRFD.ab', 'MND.ab'))+
    scale_x_continuous(limits=c(-1.5,1.5), breaks = seq(-1.5, 1.5, 0.5))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())

(domin_nd_rfd.a.all.varied.intercept.plot =
    ggplot((data=domin_data.inla.md.a_func_d.all_intercept %>% 
              filter(rowname %in% c('MRFD.ab', 'MND.ab'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 3)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MRFD.ab', 'MND.ab'),
                     label = c(expression(RFD['ab']), expression(ND['ab'])))+
    scale_x_continuous(limits=c(-0.5,0.55))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MND.ab" = "MND.ab",
                                  "MRFD.ab" = "MRFD.ab"),
                       name = ' ')+
    annotate(geom="text", x=0.26, y=2.5, family = 'Arial',
             label='Dominance', fontface = 'bold')+
    labs(x = 'Standardized effects', y = '  ')+
    guides(color="none")+
    theme_regular())

(domin_pd_fd.a.all.varied.intercept.plot =
    ggplot(data=(data=domin_data.inla.md.a_func_d.all_intercept %>% 
                   filter(rowname %in% c('MFD.ab_all','MPD.ab_all'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 3)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MFD.ab_all','MPD.ab_all'),
                     label = c(expression(MFD['ab']), expression(MPD['ab'])))+
    scale_x_continuous(limits=c(-1.05,1.05))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MPD.ab_all" = "MPD.ab_all",
                                  "MFD.ab_all" = "MFD.ab_all"),
                       name = ' ')+
    annotate(geom="text", x=0.6, y=2.5, family = 'Arial',
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
(domin_mnd.all.varied.intercept.plot =
    ggplot(data=domin_data.inla.mnd_func_d.all_intercept,aes(x=mean,y=rowname,color=rowname))+
    geom_point()+
    ggtitle('dominance')+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNFD_all','MNTD_all', 'MNRFD', 'MNND'))+
    scale_x_continuous(limits=c(-1.5,1.5), breaks = seq(-1.5, 1.5, 0.5))+
    scale_color_viridis_d()+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())

(domin_nnd_nrfd.all.varied.intercept.plot =
    ggplot((data=domin_data.inla.mnd_func_d.all_intercept %>% 
              filter(rowname %in% c('MNRFD', 'MNND'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 2)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNRFD', 'MNND'))+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_manual(values = c(colors_4d[1],
                                  colors_4d[2]),
                       labels = c("MNND" = "MNND",
                                  "MNRFD" = "MNRFD"),
                       name = ' ')+
    annotate(geom="text", x=-1.4, y=2.5, family = 'Arial',
             label='Dominance')+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())

(domin_npd_nfd.all.varied.intercept.plot =
    ggplot(data=(data=domin_data.inla.mnd_func_d.all_intercept %>% 
                   filter(rowname %in% c('MNFD_all','MNTD_all'))),
           aes(x=mean,y=rowname,color=rowname))+
    geom_point(size = 2)+
    geom_linerange(aes(xmin=lower,xmax=upper))+
    geom_vline(xintercept=0,linetype=2,color="grey40")+
    scale_y_discrete(limits=c('MNFD_all','MNTD_all'),
                     label =c('MNFD','MNTD') )+
    scale_x_continuous(limits=c(-1.5,1.5))+
    scale_color_manual(values = c(colors_4d[4],
                                  colors_4d[3]),
                       labels = c("MNTD_all" = "MNTD_all",
                                  "MNFD_all" = "MNFD_all"),
                       name = ' ')+
    annotate(geom="text", x=-1.4, y=2.5, family = 'Arial',
             label='Dominance')+
    labs(x = 'Standardized effects', y = NULL)+
    guides(color="none")+
    theme_test())

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
                                     mfd=mean(dat_dom_sp$mfd),
                                     mpd_all = mean(dat_dom_sp$mpd_all),
                                     mfunc_d_all = mean(dat_dom_sp$mfunc_d_all))

lincombs.matrix.domin.mnd=model.matrix(~mnd+mfd+mpd_all+mfunc_d_all,
                                       data=lincombs.data.domin.mnd)
lincombs.matrix.domin.mnd=as.data.frame(lincombs.matrix.domin.mnd)
lincombs.domin.mnd=inla.make.lincombs(lincombs.matrix.domin.mnd)

inla.model_lincombs.domin.mnd = pglmm(domin ~ mnd+mfd+mpd_all+mfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance~mnd
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnd.rdata")

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


### get dominance predict data for mfd
lincombs.data.domin.mfd = data.frame(mfd=seq(min(dat_dom_sp$mfd),max(dat_dom_sp$mfd),length=100),
                                       mnd=mean(dat_dom_sp$mnd),
                                       mpd_all = mean(dat_dom_sp$mpd_all),
                                       mfunc_d_all = mean(dat_dom_sp$mfunc_d_all))

lincombs.matrix.domin.mfd=model.matrix(~mnd+mfd+mpd_all+mfunc_d_all,
                                         data=lincombs.data.domin.mfd)
lincombs.matrix.domin.mfd=as.data.frame(lincombs.matrix.domin.mfd)
lincombs.domin.mfd=inla.make.lincombs(lincombs.matrix.domin.mfd)

inla.model_lincombs.domin.mfd = pglmm(domin ~ mnd+mfd+mpd_all+mfunc_d_all+(1|species) + 
                                          (1|f_p) + (1|field), data = dat_dom_sp,
                                        family = "binomial", cov_ranef = list(species = tree),
                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                    config = TRUE),
                                                             quantiles=c(0.025,0.5,0.975),
                                                             lincomb=lincombs.domin.mfd,
                                                             control.predictor=list(compute=T)),
                                        bayes = T)

lincombs.posterior.domin.mfd = inla.model_lincombs.domin.mfd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mfd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mfd$predicted.value=unlist(lapply(lincombs.posterior.domin.mfd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mfd$lower=unlist(lapply(lincombs.posterior.domin.mfd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mfd$upper=unlist(lapply(lincombs.posterior.domin.mfd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mfd
save(lincombs.data.domin.mfd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfd
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfd.rdata")
(domin.mfd.partial.logistic=ggplot(data=lincombs.data.domin.mfd,aes(x=mfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfd, y=domin),
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
                                         mfd = mean(dat_dom_sp$mfd),
                                         mfunc_d_all = mean(dat_dom_sp$mfunc_d_all))

lincombs.matrix.domin.mpd_all=model.matrix(~mnd+mfd+mpd_all+mfunc_d_all,
                                           data=lincombs.data.domin.mpd_all)
lincombs.matrix.domin.mpd_all=as.data.frame(lincombs.matrix.domin.mpd_all)
lincombs.domin.mpd_all=inla.make.lincombs(lincombs.matrix.domin.mpd_all)

inla.model_lincombs.domin.mpd_all = pglmm(domin ~ mnd+mfd+mpd_all+mfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mpd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mpd_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mpd_all.rdata")
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
                                             mfd = mean(dat_dom_sp$mfd),
                                             mfunc_d_all = seq(min(dat_dom_sp$mfunc_d_all),
                                                               max(dat_dom_sp$mfunc_d_all),length=100))

lincombs.matrix.domin.mfunc_d_all=model.matrix(~mnd+mfd+mpd_all+mfunc_d_all,
                                               data=lincombs.data.domin.mfunc_d_all)
lincombs.matrix.domin.mfunc_d_all=as.data.frame(lincombs.matrix.domin.mfunc_d_all)
lincombs.domin.mfunc_d_all=inla.make.lincombs(lincombs.matrix.domin.mfunc_d_all)

inla.model_lincombs.domin.mfunc_d_all = pglmm(domin ~ mnd+mfd+mpd_all+mfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfunc_d_all.rdata")
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
                                      mnfd=mean(dat_dom_sp$mnfd),
                                      mntd_all = mean(dat_dom_sp$mntd_all),
                                      mnfunc_d_all = mean(dat_dom_sp$mnfunc_d_all))

lincombs.matrix.domin.mnnd=model.matrix(~mnnd+mnfd+mntd_all+mnfunc_d_all,
                                        data=lincombs.data.domin.mnnd)
lincombs.matrix.domin.mnnd=as.data.frame(lincombs.matrix.domin.mnnd)
lincombs.domin.mnnd=inla.make.lincombs(lincombs.matrix.domin.mnnd)

inla.model_lincombs.domin.mnnd = pglmm(domin ~ mnnd+mnfd+mntd_all+mnfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance~mnnd
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnnd.rdata")

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


### get dominance predict data for mnfd
lincombs.data.domin.mnfd = data.frame(mnfd=seq(min(dat_dom_sp$mnfd),max(dat_dom_sp$mnfd),length=100),
                                        mnnd=mean(dat_dom_sp$mnnd),
                                        mntd_all = mean(dat_dom_sp$mntd_all),
                                        mnfunc_d_all = mean(dat_dom_sp$mnfunc_d_all))

lincombs.matrix.domin.mnfd=model.matrix(~mnnd+mnfd+mntd_all+mnfunc_d_all,
                                          data=lincombs.data.domin.mnfd)
lincombs.matrix.domin.mnfd=as.data.frame(lincombs.matrix.domin.mnfd)
lincombs.domin.mnfd=inla.make.lincombs(lincombs.matrix.domin.mnfd)

inla.model_lincombs.domin.mnfd = pglmm(domin ~ mnnd+mnfd+mntd_all+mnfunc_d_all+(1|species) + 
                                           (1|f_p) + (1|field), data = dat_dom_sp,
                                         family = "binomial", cov_ranef = list(species = tree),
                                         bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                     config = TRUE),
                                                              quantiles=c(0.025,0.5,0.975),
                                                              lincomb=lincombs.domin.mnfd,
                                                              control.predictor=list(compute=T)),
                                         bayes = T)

lincombs.posterior.domin.mnfd = inla.model_lincombs.domin.mnfd$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mnfd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mnfd$predicted.value=unlist(lapply(lincombs.posterior.domin.mnfd,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mnfd$lower=unlist(lapply(lincombs.posterior.domin.mnfd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mnfd$upper=unlist(lapply(lincombs.posterior.domin.mnfd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mnfd
save(lincombs.data.domin.mnfd, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnfd
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnfd.rdata")
(domin.mnfd.partial.logistic=ggplot(data=lincombs.data.domin.mnfd,aes(x=mnfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    theme_test()+
    theme(axis.title = element_text(face="bold"),
          legend.position = 'none')+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnfd, y=domin),
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
                                          mnfd = mean(dat_dom_sp$mnfd),
                                          mnfunc_d_all = mean(dat_dom_sp$mnfunc_d_all))

lincombs.matrix.domin.mntd_all=model.matrix(~mnnd+mnfd+mntd_all+mnfunc_d_all,
                                            data=lincombs.data.domin.mntd_all)
lincombs.matrix.domin.mntd_all=as.data.frame(lincombs.matrix.domin.mntd_all)
lincombs.domin.mntd_all=inla.make.lincombs(lincombs.matrix.domin.mntd_all)

inla.model_lincombs.domin.mntd_all = pglmm(domin ~ mnnd+mnfd+mntd_all+mnfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mntd_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mntd_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mntd_all.rdata")
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
                                              mnfd = mean(dat_dom_sp$mnfd),
                                              mnfunc_d_all = seq(min(dat_dom_sp$mnfunc_d_all),
                                                                 max(dat_dom_sp$mnfunc_d_all),length=100))

lincombs.matrix.domin.mnfunc_d_all=model.matrix(~mnnd+mnfd+mntd_all+mnfunc_d_all,
                                                data=lincombs.data.domin.mnfunc_d_all)
lincombs.matrix.domin.mnfunc_d_all=as.data.frame(lincombs.matrix.domin.mnfunc_d_all)
lincombs.domin.mnfunc_d_all=inla.make.lincombs(lincombs.matrix.domin.mnfunc_d_all)

inla.model_lincombs.domin.mnfunc_d_all = pglmm(domin ~ mnnd+mnfd+mntd_all+mnfunc_d_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnfunc_d_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnfunc_d_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnfunc_d_all.rdata")
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
                                       mfd.a=mean(dat_dom_sp$mfd.a),
                                       mpd.a_all = mean(dat_dom_sp$mpd.a_all),
                                       mfunc_d.a_all = mean(dat_dom_sp$mfunc_d.a_all))

lincombs.matrix.domin.mnd.a=model.matrix(~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all,
                                         data=lincombs.data.domin.mnd.a)
lincombs.matrix.domin.mnd.a=as.data.frame(lincombs.matrix.domin.mnd.a)
lincombs.domin.mnd.a=inla.make.lincombs(lincombs.matrix.domin.mnd.a)

inla.model_lincombs.domin.mnd.a = pglmm(domin ~ mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnd.a.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance~mnd.ab
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mnd.a.rdata")
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


### get dominance predict data for mfd.ab
lincombs.data.domin.mfd.a = data.frame(mfd.a=seq(min(dat_dom_sp$mfd.a),max(dat_dom_sp$mfd.a),length=100),
                                         mnd.a=mean(dat_dom_sp$mnd.a),
                                         mpd.a_all = mean(dat_dom_sp$mpd.a_all),
                                         mfunc_d.a_all = mean(dat_dom_sp$mfunc_d.a_all))

lincombs.matrix.domin.mfd.a=model.matrix(~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all,
                                           data=lincombs.data.domin.mfd.a)
lincombs.matrix.domin.mfd.a=as.data.frame(lincombs.matrix.domin.mfd.a)
lincombs.domin.mfd.a=inla.make.lincombs(lincombs.matrix.domin.mfd.a)

inla.model_lincombs.domin.mfd.a = pglmm(domin ~ mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
                                            (1|f_p) + (1|field), data = dat_dom_sp,
                                          family = "binomial", cov_ranef = list(species = tree),
                                          bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                      config = TRUE),
                                                               quantiles=c(0.025,0.5,0.975),
                                                               lincomb=lincombs.domin.mfd.a,
                                                               control.predictor=list(compute=T)),
                                          bayes = T)

lincombs.posterior.domin.mfd.a = inla.model_lincombs.domin.mfd.a$inla.model$marginals.lincomb.derived[c(1:100)]
inla.model_lincombs.domin.mfd.a$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.domin.nd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.domin.mfd.a$predicted.value=unlist(lapply(lincombs.posterior.domin.mfd.a,function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.domin.mfd.a$lower=unlist(lapply(lincombs.posterior.domin.mfd.a,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.domin.mfd.a$upper=unlist(lapply(lincombs.posterior.domin.mfd.a,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))
#lincombs.data.domin.mfd.a
save(lincombs.data.domin.mfd.a, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfd.a.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfd.ab
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfd.a.rdata")
(domin.mfd.a.partial.logistic=ggplot(data=lincombs.data.domin.mfd.a,aes(x=mfd.a, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(color=colors_4d[2],size=1)+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfd.a, y=domin),
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
                                           mfd.a = mean(dat_dom_sp$mfd.a),
                                           mfunc_d.a_all = mean(dat_dom_sp$mfunc_d.a_all))

lincombs.matrix.domin.mpd.a_all=model.matrix(~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all,
                                             data=lincombs.data.domin.mpd.a_all)
lincombs.matrix.domin.mpd.a_all=as.data.frame(lincombs.matrix.domin.mpd.a_all)
lincombs.domin.mpd.a_all=inla.make.lincombs(lincombs.matrix.domin.mpd.a_all)

inla.model_lincombs.domin.mpd.a_all = pglmm(domin ~ mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mpd.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mpd.a_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mpd.a_all.rdata")
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
                                               mfd.a = mean(dat_dom_sp$mfd.a),
                                               mfunc_d.a_all = seq(min(dat_dom_sp$mfunc_d.a_all),
                                                                   max(dat_dom_sp$mfunc_d.a_all),length=100))

lincombs.matrix.domin.mfunc_d.a_all=model.matrix(~mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all,
                                                 data=lincombs.data.domin.mfunc_d.a_all)
lincombs.matrix.domin.mfunc_d.a_all=as.data.frame(lincombs.matrix.domin.mfunc_d.a_all)
lincombs.domin.mfunc_d.a_all=inla.make.lincombs(lincombs.matrix.domin.mfunc_d.a_all)

inla.model_lincombs.domin.mfunc_d.a_all = pglmm(domin ~ mnd.a+mfd.a+mpd.a_all+mfunc_d.a_all+(1|species) + 
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
     file = 'code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfunc_d.a_all.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mfunc_d.a_all
load("code/results_analyzing/analysing_ages1_35_top50_aposi_data/lincombs.data.domin.mfunc_d.a_all.rdata")
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



#### Fig.S1：Demonstration of how Spaak's ND/RFD affects successful colonization and dominance of species, including presentation of raw data and biased regression results ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Define dimensions for the plots and rows
plot_width = 1 / 3  # The width for plots based on the widest row (3 plots)
normal_height = plot_width  # Maintain width-height equality
first_row_height = 1.5 * normal_height  # The first row height is 1.5 times other rows

# Calculate the total height and determine the scaling factor
total_height = first_row_height + 2 * normal_height
scale_factor = 1 / total_height  # Adjust total height to fit into a unit height of 1

# Adjust each plot’s height and width according to the scaling factor
adjusted_first_row_height = first_row_height * scale_factor
adjusted_normal_height = normal_height * scale_factor
adjusted_plot_width = plot_width * scale_factor  # Adjust width to keep proportions

# Use ggdraw and draw_plot to create the graphical layout
Fig.S1.all = ggdraw() +
  draw_plot(Fig.S1_compare_estab_original_pie,
            x = 0, y = 1 - adjusted_first_row_height, width = adjusted_plot_width * 1.5, height = adjusted_first_row_height) +
  draw_plot(Fig.S1_compare_domin_original_pie,
            x = adjusted_plot_width * 1.5, y = 1 - adjusted_first_row_height, width = adjusted_plot_width * 1.5, height = adjusted_first_row_height) +
  draw_plot(estab.mnd.a_single.logistic,
            x = 0, y = 1 - adjusted_first_row_height - adjusted_normal_height, width = adjusted_plot_width, height = adjusted_normal_height) +
  draw_plot(estab.mfd.a_single.logistic,
            x = adjusted_plot_width, y = 1 - adjusted_first_row_height - adjusted_normal_height, width = adjusted_plot_width, height = adjusted_normal_height) +
  draw_plot(estab_nd_rfd.a.all.varied.intercept.plot,
            x = 2 * adjusted_plot_width, y = 1 - adjusted_first_row_height - adjusted_normal_height, width = adjusted_plot_width, height = adjusted_normal_height) +
  draw_plot(domin.mnd.a_single.logistic,
            x = 0, y = 1 - adjusted_first_row_height - 1.9 * adjusted_normal_height, width = adjusted_plot_width, height = adjusted_normal_height) +
  draw_plot(domin_nd_rfd.a.all.varied.intercept.plot,
            x = 2 * adjusted_plot_width, y = 1 - adjusted_first_row_height - 1.9 * adjusted_normal_height, width = adjusted_plot_width, height = adjusted_normal_height)+
  draw_plot(domin.mfd.a_single.logistic,
            x = adjusted_plot_width, y = 1 - adjusted_first_row_height - 1.9 * adjusted_normal_height, width = adjusted_plot_width, height = adjusted_normal_height) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"),
                  x = c(0, 1.5 * adjusted_plot_width,
                        0, adjusted_plot_width, 2 * adjusted_plot_width,
                        0, adjusted_plot_width, 2 * adjusted_plot_width),
                  y = c(1, 1,
                        1 - adjusted_first_row_height,
                        1 - adjusted_first_row_height,
                        1 - adjusted_first_row_height,
                        1 - 2.4 * adjusted_normal_height,
                        1 - 2.4 * adjusted_normal_height,
                        1 - 2.4 * adjusted_normal_height),
                  hjust = -0.1, vjust = 1, size = 15)

Fig.S1.all
emf('results/figures_ages1_35_top50_aposi/Fig.S1.all_pie.emf',
    width = 25.2, height = 25.2, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S1.all
dev.off() #turn off device and finalize file


#### Fig.2：Demonstration of how PD/FD affects successful colonization and dominance of species, including presentation of raw data and biased regression results ####
library(ggpubr)
library(export)
library(devEMF)
library(cowplot)

# Adjust grid settings
num_rows <- 2
num_cols <- 3
plot_width <- 1 / num_cols
plot_height <- 1 / num_rows

# Create a ggdraw object with all plots
Fig.2.all = ggdraw() +
  draw_plot(estab.mpd.a_single.logistic, x = 0, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab.mfunc_d.a_single.logistic, x = plot_width, y = 0.5,
            width = plot_width, height = plot_height) +
  draw_plot(estab_pd_fd.a.all.varied.intercept.plot, x = 2 * plot_width, y = 0.5,
            width = plot_width,height = plot_height) +
  draw_plot(domin.mpd.a_single.logistic, x = 0, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot(domin.mfunc_d.a_single.logistic, x = plot_width, y = 0.05, 
            width = plot_width, height = plot_height) +
  draw_plot(domin_pd_fd.a.all.varied.intercept.plot, x = 2 * plot_width, y = 0.05,
            width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
    x = c(0, 1 * plot_width, 2 * plot_width, 0, 1 * plot_width, 2 * plot_width),
    y = c(1, 1, 1, 0.55, 0.55, 0.55),
    hjust = 0, vjust = 1.1, size = 15
  )


Fig.2.all
emf('results/figures_ages1_35_top50_aposi/Fig.2.all.emf',
    width = 25.2, height = 16.2, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.2.all
dev.off() #turn off device and finalize file

Fig.2.md.a = ggdraw() +
  draw_plot(estab.mpd.a_all.partial.logistic, 0.02, 0.5, .48, .46) +
  draw_plot(estab.mfunc_d.a_all.partial.logistic, 0.5, 0.5, .48, .46) +
  draw_plot(domin.mpd.a_all.partial.logistic, 0.02, 0.08, .48, .46) +
  draw_plot(domin.mfunc_d.a_all.partial.logistic, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12
  )
emf('results/figures_ages1_35_top50_aposi/all_species_d/Fig.2.md.a.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.2.md.a
dev.off() #turn off device and finalize file

Fig.2.md.a_single = ggdraw() +
  draw_plot(estab.mpd.a_single.logistic, 0.02, 0.5, .48, .46) +
  draw_plot(estab.mfunc_d.a_single.logistic, 0.5, 0.5, .48, .46) +
  draw_plot(domin.mpd.a_single.logistic, 0.02, 0.08, .48, .46) +
  draw_plot(domin.mfunc_d.a_single.logistic, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12
  )
emf('results/figures_ages1_35_top50_aposi/all_species_d/Fig.2.md.a_single.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.2.md.a_single
dev.off() #turn off device and finalize file


Fig.2.md = ggdraw() +
  draw_plot(estab.mpd_all.partial.logistic, 0.02, 0.5, .48, .46) +
  draw_plot(estab.mfunc_d_all.partial.logistic, 0.5, 0.5, .48, .46) +
  draw_plot(domin.mpd_all.partial.logistic, 0.02, 0.08, .48, .46) +
  draw_plot(domin.mfunc_d_all.partial.logistic, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12
  )
emf('results/figures_ages1_35_top50_aposi/all_species_d/Fig.2.md.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.2.md
dev.off() #turn off device and finalize file

#### SEM: Success ~ mND mRFD ~ mfd mpd for fitted species ####
### Original model ###
require(glmmTMB)
library(piecewiseSEM)
require(optimx)
require(dplyr)

### SEM for establishment
setwd("D:/R projects/BSS")
load('code/results_analyzing/analysing_ages1_35_top50_aposi_data/dat_suc_sp.rdata')

numcols = grep("^m.",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)


estab_sem1_md = psem(
  glmmTMB(estab~mnd+mfd+mpd+mfunc_d+
            (1|species)+
            (1|field)
          #+(1|f_p)
          , family=binomial, data=dat_suc_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))), 
  glmmTMB(mnd~mpd+
            #mfunc_d+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  glmmTMB(mfd~#mpd+
            mfunc_d+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  mfd %~~% mnd,
  #estab %~~% mpd,
  #estab %~~% mfunc_d,
  data=dat_suc_sps)

summary(estab_sem1_md)
print(fisherC(estab_sem1_md), digits = 9)
AIC(estab_sem1_md)
coefs(estab_sem1_md)

estab_sem1_md.a = psem(
  glmmTMB(estab~mnd.a+mfd.a+mpd.a+mfunc_d.a+
            (1|species)+
            (1|field)
          #+(1|f_p)
          , family=binomial, data=dat_suc_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a+mfunc_d.a+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data = dat_suc_sps),
  glmmTMB(mfd.a~mpd.a+
            #mfunc_d.a+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  mfd.a %~~% mnd.a,
  #estab %~~% mpd,
  #estab %~~% mfunc_d,
  data=dat_suc_sps)

summary(estab_sem1_md.a)
print(fisherC(estab_sem1_md.a), digits = 9)
AIC(estab_sem1_md.a)

estab_sem1_mnd = psem(
  glmmTMB(estab~mnnd+mnfd+
            mntd+
            mnfunc_d+
            (1|species)+
            (1|field)
          #+(1|f_p)
          , family=binomial, data=dat_suc_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))
  ),
  glmmTMB(mnnd~#mntd+
            mnfunc_d+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  glmmTMB(mnfd~mntd+
            mnfunc_d+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  mnfd %~~% mnnd,
  #estab %~~% mpd,
  #estab %~~% mfunc_d,
  data=dat_suc_sps)

summary(estab_sem1_mnd)
print(fisherC(estab_sem1_mnd), digits = 9)
AIC(estab_sem1_mnd)


### SEM for Dominance ###
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c("establish", "dominant"))

domin_sem1_md = psem(
  glmmTMB(domin~#mnd+
            mfd+
            #mpd+
            #mfunc_d+
            (1|species)+
            (1|field)
          +(1|f_p)
          , family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd~mpd+
            #mfunc_d+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_dom_sps),
  glmmTMB(mfd~#mpd+
            mfunc_d+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_dom_sps),
  mfd %~~% mnd,
  #domin %~~% mpd,
  #domin %~~% mfunc_d,
  data=dat_dom_sps)

summary(domin_sem1_md)
print(fisherC(domin_sem1_md), digits = 9)
AIC(domin_sem1_md)

domin_sem1_md.a = psem(
  glmmTMB(domin~mnd.a+mfd.a+mpd.a+mfunc_d.a+
            (1|species)+
            (1|field)+
            (1|f_p), family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a+#mfunc_d.a+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data = dat_dom_sps),
  glmmTMB(mfd.a~#mpd.a+
            mfunc_d.a+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_dom_sps),
  mfd.a %~~% mnd.a,
  #domin %~~% mpd,
  #domin %~~% mfunc_d,
  data=dat_dom_sps)

summary(domin_sem1_md.a)
print(fisherC(domin_sem1_md.a), digits = 9)
AIC(domin_sem1_md.a)

domin_sem1_mnd = psem(
  glmmTMB(domin~mnnd+mnfd+
            mntd+
            mnfunc_d+
            (1|species)+
            (1|field)+
            (1|f_p), family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnnd~mntd+
            #mnfunc_d+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_dom_sps),
  glmmTMB(mnfd~#mntd+
            mnfunc_d+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_dom_sps),
  mnfd %~~% mnnd,
  #domin %~~% mpd,
  #domin %~~% mfunc_d,
  data=dat_dom_sps)

summary(domin_sem1_mnd)
print(fisherC(domin_sem1_mnd), digits = 9)
AIC(domin_sem1_mnd)

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
sd.y_estab_sem1_md = sqrt(predict(estab_sem1_md[[1]],re.form=NA)%>%var+sum(as.numeric(unlist(VarCorr(estab_sem1_md[[1]]))))+pi^2/3)
sd.y_estab_sem1_md
std.estab_sem1_md=unlist(fixef(estab_sem1_md[[1]]))[-1]*
  c(sd(dat_suc_sps$mpd),
    sd(dat_suc_sps$mfunc_d),
    sd(dat_suc_sps$mnd),
    sd(dat_suc_sps$mfd))/sd.y_estab_sem1_md
std.estab_sem1_md

# merge the new corrected cofficients into the original cofficients data frame
coefs_estab_sem1_md=coefs(estab_sem1_md)# original coefficients of piecewiseSEM
coefs_estab_sem1_md[1:4,8]=std.estab_sem1_md

### Dominance corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_domin_sem1_md=sqrt(predict(domin_sem1_md[[1]],re.form=NA)%>%var+sum(as.numeric(unlist(VarCorr(domin_sem1_md[[1]]))))+pi^2/3)
sd.y_domin_sem1_md
std.domin_sem1_md=unlist(fixef(domin_sem1_md[[1]]))[-1]*
  c(sd(dat_dom_sps$mpd),
    sd(dat_dom_sps$mconti_func_d),
    sd(dat_dom_sps$mnd),
    sd(dat_dom_sps$mfd))/sd.y_domin_sem1_md
std.domin_sem1_md

# merge the new corrected cofficients into the original cofficients data frame
coefs_domin_sem1_md=coefs(domin_sem1_md)# original coefficients of piecewiseSEM
coefs_domin_sem1_md[1:4,8]=std.domin_sem1_md
coefs_domin_sem1_md

### Calculate direct and indirect effects of SEM following Xu Meng (2022) GCB
### Estab
sem_estab_effects_md = data.frame(predictor = rep(unique(coefs_estab_sem1_md$Predictor)[1:4], 3),
                                  type = rep(c('direct', 'indirect', 'total'), each = 4),
                                  value = rep(0, 12),
                                  significant = rep('yes', 12))
coefs_estab_sem1_md = coefs_estab_sem1_md[-nrow(coefs_domin_sem1_md),]
coefs_estab_sem1_md_c = coefs_estab_sem1_md[coefs_estab_sem1_md$P.Value < 0.05,]
for (i in 1:length(unique(coefs_estab_sem1_md_c$Predictor)[1:nrow(coefs_estab_sem1_md_c)])) {
  #i = 4
  predictor_1 = unique(coefs_estab_sem1_md_c$Predictor)[1:4][i]
  dat = coefs_estab_sem1_md_c[coefs_estab_sem1_md_c$Predictor == predictor_1,]
  dat_2 = coefs_estab_sem1_md_c[coefs_estab_sem1_md_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_estab_effects_md[sem_estab_effects_md$predictor == predictor_1&
                           sem_estab_effects_md$type == 'direct',]$value = dat_3$Std.Estimate
    sem_estab_effects_md[sem_estab_effects_md$predictor == predictor_1&
                           sem_estab_effects_md$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    sem_estab_effects_md[sem_estab_effects_md$predictor == predictor_1&
                           sem_estab_effects_md$type == 'direct',]$value = dat_3[dat_3$Response == 'estab' &
                                                                                   dat_3$Predictor == predictor_1,]$Std.Estimate
    sem_estab_effects_md[sem_estab_effects_md$predictor == predictor_1&
                           sem_estab_effects_md$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'estab',]$Std.Estimate*dat_3[dat_3$Response == 'estab' &
                                                                                                                                        dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_estab_effects_md[sem_estab_effects_md$predictor == predictor_1&
                           sem_estab_effects_md$type == 'total',]$value  = sem_estab_effects_md[sem_estab_effects_md$predictor == predictor_1&
                                                                                                  sem_estab_effects_md$type == 'direct',]$value + sem_estab_effects_md[sem_estab_effects_md$predictor == predictor_1&
                                                                                                                                                                         sem_estab_effects_md$type == 'indirect',]$value
  }
}

### Domin
sem_domin_effects_md = data.frame(predictor = rep(unique(coefs_domin_sem1_md$Predictor)[1:4], 3),
                                  type = rep(c('direct', 'indirect', 'total'), each = 4),
                                  value = rep(0, 12))
coefs_domin_sem1_md = coefs_domin_sem1_md[-nrow(coefs_domin_sem1_md),]
coefs_domin_sem1_md_c = coefs_domin_sem1_md[coefs_domin_sem1_md$P.Value < 0.05,]
for (i in 1:length(unique(coefs_domin_sem1_md_c$Predictor)[1:nrow(coefs_domin_sem1_md_c)])) {
  #i = 4
  predictor_1 = unique(coefs_domin_sem1_md_c$Predictor)[1:nrow(coefs_domin_sem1_md_c)][i]
  dat = coefs_domin_sem1_md_c[coefs_domin_sem1_md_c$Predictor == predictor_1,]
  dat_2 = coefs_domin_sem1_md_c[coefs_domin_sem1_md_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_domin_effects_md[sem_domin_effects_md$predictor == predictor_1&
                           sem_domin_effects_md$type == 'direct',]$value = dat_3$Std.Estimate
    sem_domin_effects_md[sem_domin_effects_md$predictor == predictor_1&
                           sem_domin_effects_md$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    sem_domin_effects_md[sem_domin_effects_md$predictor == predictor_1&
                           sem_domin_effects_md$type == 'direct',]$value = dat_3[dat_3$Response == 'domin' &
                                                                                   dat_3$Predictor == predictor_1,]$Std.Estimate
    sem_domin_effects_md[sem_domin_effects_md$predictor == predictor_1&
                           sem_domin_effects_md$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'domin',]$Std.Estimate*dat_3[dat_3$Response == 'domin' &
                                                                                                                                        dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_domin_effects_md[sem_domin_effects_md$predictor == predictor_1&
                           sem_domin_effects_md$type == 'total',]$value  = sem_domin_effects_md[sem_domin_effects_md$predictor == predictor_1&
                                                                                                  sem_domin_effects_md$type == 'direct',]$value + sem_domin_effects_md[sem_domin_effects_md$predictor == predictor_1&
                                                                                                                                                                         sem_domin_effects_md$type == 'indirect',]$value
  }
}

#-----------Draw plots for SEMs' relative total effects, direct and indirect effects 
### Estab
require(ggthemes)
require(ggplot2)
require(ggpubr)

sem_estab_effects_md_total = sem_estab_effects_md %>% filter(type == 'total')
sem_estab_effects_md.percent = sem_estab_effects_md_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_estab_effects_md.longer = sem_estab_effects_md %>% filter(type != 'total' 
                                                              #& value != 0
)

# Relative total effect 
sem_estab_effects_md.percent.plot = 
  ggplot(data=sem_estab_effects_md.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mfunc_d", "mpd", "mfd", "mnd"),
                   labels=c("I-N MFD", "I-N MPD", "I-N MRFD", "I-N MND"),
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
sem_estab_effects_md.dir_indir.mpd.plot = 
  ggplot(data=sem_estab_effects_md.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects","indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mfunc_d", "mpd", "mfd", "mnd"),
                   position="right")+
  scale_x_continuous(breaks = c(seq(-0.4, 0.1, 0.2)))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_estab_md_space = ggplot()+
  facet_wrap(~"Establishment")+
  theme_test() + 
  theme(strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.015,1.2,0.2),units="lines"))

# Merge three plots
sem_estab_md_plot_all = ggarrange(sem_estab_md_space,
                                  sem_estab_effects_md.dir_indir.mpd.plot,
                                  sem_estab_effects_md.percent.plot,
                                  nrow = 1, ncol = 3,
                                  labels = c('(a)', '', ''),
                                  vjust = 1.8)

emf('results/figures_ages1_35_top50_aposi/fitted_species_d/sem_estab_md_plot_all.emf',
    width=9, height=3, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
sem_estab_md_plot_all
dev.off() #turn off device and finalize file


### domin
sem_domin_effects_md_total = sem_domin_effects_md %>% filter(type == 'total')
sem_domin_effects_md.percent = sem_domin_effects_md_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_domin_effects_md.longer = sem_domin_effects_md %>% filter(type != 'total' 
                                                              #&value != 0
)

# Relative total effect
sem_domin_effects_md.percent.plot = 
  ggplot(data=sem_domin_effects_md.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mfunc_d", "mpd", "mfd", "mnd"),
                   labels=c("E-N MFD", "E-N MPD", "E-N MRFD", "E-N MND"),
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
sem_domin_effects_md.dir_indir.mpd.plot = 
  ggplot(data=sem_domin_effects_md.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects",
                                          "indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mfunc_d", "mpd", "mfd", "mnd"),
                   position="right")+
  scale_x_continuous(breaks = seq(-0.1, 0.1, 0.1),
                     limits = c(-0.1, 0.15))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_domin_md_space = ggplot()+
  facet_wrap(~"Dominance")+
  theme_test() + 
  theme(strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.015,1.2,0.2),units="lines"))

# Merge three plots
sem_domin_md_plot_all = ggarrange(sem_domin_md_space,
                                  sem_domin_effects_md.dir_indir.mpd.plot,
                                  sem_domin_effects_md.percent.plot,
                                  nrow = 1, ncol = 3,
                                  labels = c('(b)', '', ''),
                                  vjust = 1.8)


ggsave(plot = sem_domin_md_plot_all,
       "results/figures_ages1_35_top50_aposi/fitted_species_d/sem_domin_md_plot_all.svg",
       width=9,height=3)


##### Draw plot for md.a
### Establishment corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_estab_sem1_md.a = sqrt(predict(estab_sem1_md.a[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(estab_sem1_md.a[[1]])))+pi^2/3)
sd.y_estab_sem1_md.a
std.estab_sem1_md.a=unlist(fixef(estab_sem1_md.a[[1]]))[-1]*
  c(sd(dat_suc_sps$mpd.a),
    sd(dat_suc_sps$mfunc_d.a),
    sd(dat_suc_sps$mnd.a),
    sd(dat_suc_sps$mfd.a))/sd.y_estab_sem1_md.a
std.estab_sem1_md.a

# merge the new corrected cofficients into the original cofficients data frame
coefs_estab_sem1_md.a=coefs(estab_sem1_md.a)# original coefficients of piecewiseSEM
coefs_estab_sem1_md.a[1:4,8]=std.estab_sem1_md.a

### Dominance corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_domin_sem1_md.a=sqrt(predict(domin_sem1_md.a[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(domin_sem1_md.a[[1]])))+pi^2/3)
sd.y_domin_sem1_md.a
std.domin_sem1_md.a=unlist(fixef(domin_sem1_md.a[[1]]))[-1]*
  c(sd(dat_dom_sps$mpd.a),
    sd(dat_dom_sps$mfunc_d.a),
    sd(dat_dom_sps$mnd.a),
    sd(dat_dom_sps$mfd.a))/sd.y_domin_sem1_md.a
std.domin_sem1_md.a

# merge the new corrected cofficients into the original cofficients data frame
coefs_domin_sem1_md.a=coefs(domin_sem1_md.a)# original coefficients of piecewiseSEM
coefs_domin_sem1_md.a[1:4,8]=std.domin_sem1_md.a
coefs_domin_sem1_md.a

### Calculate direct and indirect effects of SEM following Xu Meng (2022) GCB
### Estab
sem_estab_effects_md.a = data.frame(predictor = rep(unique(coefs_estab_sem1_md.a$Predictor)[1:4], 3),
                                    type = rep(c('direct', 'indirect', 'total'), each = 4),
                                    value = rep(0, 12),
                                    significant = rep('yes', 12))
coefs_estab_sem1_md.a = coefs_estab_sem1_md.a[-nrow(coefs_estab_sem1_md.a),]
coefs_estab_sem1_md.a_c = coefs_estab_sem1_md.a[coefs_estab_sem1_md.a$P.Value < 0.05,]
for (i in 1:length(unique(coefs_estab_sem1_md.a_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_estab_sem1_md.a_c$Predictor)[i]
  dat = coefs_estab_sem1_md.a_c[coefs_estab_sem1_md.a_c$Predictor == predictor_1,]
  dat_2 = coefs_estab_sem1_md.a_c[coefs_estab_sem1_md.a_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_estab_effects_md.a[sem_estab_effects_md.a$predictor == predictor_1&
                             sem_estab_effects_md.a$type == 'direct',]$value = dat_3$Std.Estimate
    sem_estab_effects_md.a[sem_estab_effects_md.a$predictor == predictor_1&
                             sem_estab_effects_md.a$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    if (nrow(dat_3[dat_3$Response == 'estab' &
                   dat_3$Predictor == predictor_1,]) > 0) {
      sem_estab_effects_md.a[sem_estab_effects_md.a$predictor == predictor_1&
                               sem_estab_effects_md.a$type == 'direct',]$value = dat_3[dat_3$Response == 'estab' &
                                                                                         dat_3$Predictor == predictor_1,]$Std.Estimate
    }
    
    sem_estab_effects_md.a[sem_estab_effects_md.a$predictor == predictor_1&
                             sem_estab_effects_md.a$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'estab',]$Std.Estimate*dat_3[dat_3$Response == 'estab' &
                                                                                                                                            dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_estab_effects_md.a[sem_estab_effects_md.a$predictor == predictor_1&
                             sem_estab_effects_md.a$type == 'total',]$value  = sem_estab_effects_md.a[sem_estab_effects_md.a$predictor == predictor_1&
                                                                                                        sem_estab_effects_md.a$type == 'direct',]$value + sem_estab_effects_md.a[sem_estab_effects_md.a$predictor == predictor_1&
                                                                                                                                                                                   sem_estab_effects_md.a$type == 'indirect',]$value
  }
}

### Domin
sem_domin_effects_md.a = data.frame(predictor = rep(unique(coefs_domin_sem1_md.a$Predictor)[1:4], 3),
                                    type = rep(c('direct', 'indirect', 'total'), each = 4),
                                    value = rep(0, 12))
coefs_domin_sem1_md.a = coefs_domin_sem1_md.a[-nrow(coefs_domin_sem1_md.a),]
coefs_domin_sem1_md.a_c = coefs_domin_sem1_md.a[coefs_domin_sem1_md.a$P.Value < 0.05,]
for (i in 1:length(unique(coefs_domin_sem1_md.a_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_domin_sem1_md.a_c$Predictor)[i]
  dat = coefs_domin_sem1_md.a_c[coefs_domin_sem1_md.a_c$Predictor == predictor_1,]
  dat_2 = coefs_domin_sem1_md.a_c[coefs_domin_sem1_md.a_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_domin_effects_md.a[sem_domin_effects_md.a$predictor == predictor_1&
                             sem_domin_effects_md.a$type == 'direct',]$value = dat_3$Std.Estimate
    sem_domin_effects_md.a[sem_domin_effects_md.a$predictor == predictor_1&
                             sem_domin_effects_md.a$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    sem_domin_effects_md.a[sem_domin_effects_md.a$predictor == predictor_1&
                             sem_domin_effects_md.a$type == 'direct',]$value = dat_3[dat_3$Response == 'domin' &
                                                                                       dat_3$Predictor == predictor_1,]$Std.Estimate
    sem_domin_effects_md.a[sem_domin_effects_md.a$predictor == predictor_1&
                             sem_domin_effects_md.a$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'domin',]$Std.Estimate*dat_3[dat_3$Response == 'domin' &
                                                                                                                                            dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_domin_effects_md.a[sem_domin_effects_md.a$predictor == predictor_1&
                             sem_domin_effects_md.a$type == 'total',]$value  = sem_domin_effects_md.a[sem_domin_effects_md.a$predictor == predictor_1&
                                                                                                        sem_domin_effects_md.a$type == 'direct',]$value + sem_domin_effects_md.a[sem_domin_effects_md.a$predictor == predictor_1&
                                                                                                                                                                                   sem_domin_effects_md.a$type == 'indirect',]$value
  }
}

#-----------Draw plots for SEMs' relative total effects, direct and indirect effects 
### Estab
require(ggthemes)
sem_estab_effects_md.a_total = sem_estab_effects_md.a %>% filter(type == 'total')
sem_estab_effects_md.a.percent = sem_estab_effects_md.a_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_estab_effects_md.a.longer = sem_estab_effects_md.a %>% filter(type != 'total' 
                                                                  #& value != 0
)

# Relative total effect 
sem_estab_effects_md.a.percent.plot = 
  ggplot(data=sem_estab_effects_md.a.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mfunc_d.a", "mpd.a", "mfd.a", "mnd.a"),
                   labels=c("FD", "PD", "RFD", "ND"),
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
sem_estab_effects_md.a.dir_indir.mpd.plot = 
  ggplot(data=sem_estab_effects_md.a.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects","indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mfunc_d.a", "mpd.a", "mfd.a", "mnd.a"),
                   position="right")+
  scale_x_continuous(breaks = c(seq(-0.4, 0.1, 0.2)))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
library(cowplot)
sem_estab_md.a_space = ggplot()+
  theme_nothing() + ## 
  theme(strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.015,1.2,0.2),units="lines"))

# Merge three plots
sem_estab_md.a_plot_all = ggarrange(sem_estab_md.a_space,
                                    sem_estab_effects_md.a.dir_indir.mpd.plot,
                                    sem_estab_effects_md.a.percent.plot,
                                    nrow = 1, ncol = 3,
                                    labels = c('(a)', '', ''),
                                    vjust = 1.8)


emf('results/figures_ages1_35_top50_aposi/fitted_species_d/sem_estab_md.a_plot_all.emf',
    width=18, height=6, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
sem_estab_md.a_plot_all
dev.off() #turn off device and finalize file



### domin
sem_domin_effects_md.a_total = sem_domin_effects_md.a %>% filter(type == 'total')
sem_domin_effects_md.a.percent = sem_domin_effects_md.a_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_domin_effects_md.a.longer = sem_domin_effects_md.a %>% filter(type != 'total' 
                                                                  #&value != 0
)

# Relative total effect
sem_domin_effects_md.a.percent.plot = 
  ggplot(data=sem_domin_effects_md.a.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mfunc_d.a", "mpd.a", "mfd.a", "mnd.a"),
                   labels=c("FD", "PD", "RFD", "ND"),
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
sem_domin_effects_md.a.dir_indir.mpd.plot = 
  ggplot(data=sem_domin_effects_md.a.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects",
                                          "indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mfunc_d.a", "mpd.a", "mfd.a", "mnd.a"),
                   position="right")+
  scale_x_continuous(breaks = seq(-0.1, 0.1, 0.1),
                     limits = c(-0.1, 0.3))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_domin_md.a_space = ggplot()+
  facet_wrap(~"Dominance")+
  theme_test() + 
  theme(strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.015,1.2,0.2),units="lines"))

# Merge three plots
sem_domin_md.a_plot_all = ggarrange(sem_domin_md.a_space,
                                    sem_domin_effects_md.a.dir_indir.mpd.plot,
                                    sem_domin_effects_md.a.percent.plot,
                                    nrow = 1, ncol = 3,
                                    labels = c('(b)', '', ''),
                                    vjust = 1.8)

emf('results/figures_ages1_35_top50_aposi/fitted_species_d/sem_domin_md.a_plot_all.emf',
    width=18, height=6, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
sem_domin_md.a_plot_all
dev.off() #turn off device and finalize file



#### Draw plot for mnd  
### Establishment corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_estab_sem1_mnd = sqrt(predict(estab_sem1_mnd[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(estab_sem1_mnd[[1]])))+pi^2/3)
sd.y_estab_sem1_mnd
std.estab_sem1_mnd=unlist(fixef(estab_sem1_mnd[[1]]))[-1]*
  c(sd(dat_suc_sps$mntd),
    sd(dat_suc_sps$mnfunc_d),
    sd(dat_suc_sps$mnnd),
    sd(dat_suc_sps$mnfd))/sd.y_estab_sem1_mnd
std.estab_sem1_mnd

# merge the new corrected cofficients into the original cofficients data frame
coefs_estab_sem1_mnd=coefs(estab_sem1_mnd)# original coefficients of piecewiseSEM
coefs_estab_sem1_mnd[1:4,8]=std.estab_sem1_mnd

### Dominance corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_domin_sem1_mnd=sqrt(predict(domin_sem1_mnd[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(domin_sem1_mnd[[1]])))+pi^2/3)
sd.y_domin_sem1_mnd
std.domin_sem1_mnd=unlist(fixef(domin_sem1_mnd[[1]]))[-1]*
  c(sd(dat_dom_sps$mntd),
    sd(dat_dom_sps$mnfunc_d),
    sd(dat_dom_sps$mnnd),
    sd(dat_dom_sps$mnfd))/sd.y_domin_sem1_mnd
std.domin_sem1_mnd

# merge the new corrected cofficients into the original cofficients data frame
coefs_domin_sem1_mnd=coefs(domin_sem1_mnd)# original coefficients of piecewiseSEM
coefs_domin_sem1_mnd[1:4,8]=std.domin_sem1_mnd
coefs_domin_sem1_mnd

### Calculate direct and indirect effects of SEM following Xu Meng (2022) GCB
### Estab
sem_estab_effects_mnd = data.frame(predictor = rep(unique(coefs_estab_sem1_mnd$Predictor)[1:4], 3),
                                   type = rep(c('direct', 'indirect', 'total'), each = 4),
                                   value = rep(0, 12),
                                   significant = rep('yes', 12))
coefs_estab_sem1_mnd = coefs_estab_sem1_mnd[-nrow(coefs_estab_sem1_mnd), ]
coefs_estab_sem1_mnd_c = coefs_estab_sem1_mnd[coefs_estab_sem1_mnd$P.Value < 0.05,]
for (i in 1:length(unique(coefs_estab_sem1_mnd_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_estab_sem1_mnd_c$Predictor)[i]
  dat = coefs_estab_sem1_mnd_c[coefs_estab_sem1_mnd_c$Predictor == predictor_1,]
  dat_2 = coefs_estab_sem1_mnd_c[coefs_estab_sem1_mnd_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_estab_effects_mnd[sem_estab_effects_mnd$predictor == predictor_1&
                            sem_estab_effects_mnd$type == 'direct',]$value = dat_3$Std.Estimate
    sem_estab_effects_mnd[sem_estab_effects_mnd$predictor == predictor_1&
                            sem_estab_effects_mnd$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    if (nrow(dat_3[dat_3$Response == 'estab' &
                   dat_3$Predictor == predictor_1,]) > 0){
      sem_estab_effects_mnd[sem_estab_effects_mnd$predictor == predictor_1&
                              sem_estab_effects_mnd$type == 'direct',]$value = dat_3[dat_3$Response == 'estab' &
                                                                                       dat_3$Predictor == predictor_1,]$Std.Estimate
    }
    sem_estab_effects_mnd[sem_estab_effects_mnd$predictor == predictor_1&
                            sem_estab_effects_mnd$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'estab',]$Std.Estimate*dat_3[dat_3$Response == 'estab' &
                                                                                                                                          dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_estab_effects_mnd[sem_estab_effects_mnd$predictor == predictor_1&
                            sem_estab_effects_mnd$type == 'total',]$value  = sem_estab_effects_mnd[sem_estab_effects_mnd$predictor == predictor_1&
                                                                                                     sem_estab_effects_mnd$type == 'direct',]$value + sem_estab_effects_mnd[sem_estab_effects_mnd$predictor == predictor_1&
                                                                                                                                                                              sem_estab_effects_mnd$type == 'indirect',]$value
  }
}

### Domin
sem_domin_effects_mnd = data.frame(predictor = rep(unique(coefs_domin_sem1_mnd$Predictor)[1:4], 3),
                                   type = rep(c('direct', 'indirect', 'total'), each = 4),
                                   value = rep(0, 12))
coefs_domin_sem1_mnd = coefs_domin_sem1_mnd[-nrow(coefs_domin_sem1_mnd),]
coefs_domin_sem1_mnd_c = coefs_domin_sem1_mnd[coefs_domin_sem1_mnd$P.Value < 0.05,]
for (i in 1:length(unique(coefs_domin_sem1_mnd_c$Predictor))) {
  #i = 4
  predictor_1 = unique(coefs_domin_sem1_mnd_c$Predictor)[i]
  dat = coefs_domin_sem1_mnd_c[coefs_domin_sem1_mnd_c$Predictor == predictor_1,]
  dat_2 = coefs_domin_sem1_mnd_c[coefs_domin_sem1_mnd_c$Predictor %in% dat$Response,]
  dat_3 = rbind(dat, dat_2)
  if (nrow(dat_3) == 1) {
    sem_domin_effects_mnd[sem_domin_effects_mnd$predictor == predictor_1&
                            sem_domin_effects_mnd$type == 'direct',]$value = dat_3$Std.Estimate
    sem_domin_effects_mnd[sem_domin_effects_mnd$predictor == predictor_1&
                            sem_domin_effects_mnd$type == 'total',]$value = dat_3$Std.Estimate
  } else if(nrow(dat_3) != 1) {
    if (nrow(dat_3[dat_3$Response == 'domin' &
                   dat_3$Predictor == predictor_1,]) > 0) {
      sem_domin_effects_mnd[sem_domin_effects_mnd$predictor == predictor_1&
                              sem_domin_effects_mnd$type == 'direct',]$value = dat_3[dat_3$Response == 'domin' &
                                                                                       dat_3$Predictor == predictor_1,]$Std.Estimate
    }
    sem_domin_effects_mnd[sem_domin_effects_mnd$predictor == predictor_1&
                            sem_domin_effects_mnd$type == 'indirect',]$value = sum(dat_3[dat_3$Response != 'domin',]$Std.Estimate*dat_3[dat_3$Response == 'domin' &
                                                                                                                                          dat_3$Predictor != predictor_1,]$Std.Estimate)
    sem_domin_effects_mnd[sem_domin_effects_mnd$predictor == predictor_1&
                            sem_domin_effects_mnd$type == 'total',]$value  = sem_domin_effects_mnd[sem_domin_effects_mnd$predictor == predictor_1&
                                                                                                     sem_domin_effects_mnd$type == 'direct',]$value + sem_domin_effects_mnd[sem_domin_effects_mnd$predictor == predictor_1&
                                                                                                                                                                              sem_domin_effects_mnd$type == 'indirect',]$value
  }
}

#-----------Draw plots for SEMs' relative total effects, direct and indirect effects 
### Estab
require(ggthemes)
sem_estab_effects_mnd_total = sem_estab_effects_mnd %>% filter(type == 'total')
sem_estab_effects_mnd.percent = sem_estab_effects_mnd_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_estab_effects_mnd.longer = sem_estab_effects_mnd %>% filter(type != 'total' 
                                                                #& value != 0
)

# Relative total effect 
sem_estab_effects_mnd.percent.plot = 
  ggplot(data=sem_estab_effects_mnd.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mnfunc_d", "mntd", "mnfd", "mnnd"),
                   labels=c("I-N MNFD", "I-N MNTD", "I-N MNRFD", "I-N MNND"),
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
sem_estab_effects_mnd.dir_indir.mpd.plot = 
  ggplot(data=sem_estab_effects_mnd.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects","indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mnfunc_d", "mntd", "mnfd", "mnnd"),
                   position="right")+
  scale_x_continuous(breaks = c(seq(-0.4, 0.1, 0.2)))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_estab_mnd_space = ggplot()+
  facet_wrap(~"Establishment")+
  theme_test() + 
  theme(strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.015,1.2,0.2),units="lines"))

# Merge three plots
sem_estab_mnd_plot_all = ggarrange(sem_estab_mnd_space,
                                   sem_estab_effects_mnd.dir_indir.mpd.plot,
                                   sem_estab_effects_mnd.percent.plot,
                                   nrow = 1, ncol = 3,
                                   labels = c('(a)', '', ''),
                                   vjust = 1.8)


ggsave(plot = sem_estab_mnd_plot_all,
       "results/figures_ages1_35_top50_aposi/fitted_species_d/sem_estab_mnd_plot_all.svg",
       width=9,height=3)

### domin
sem_domin_effects_mnd_total = sem_domin_effects_mnd %>% filter(type == 'total')
sem_domin_effects_mnd.percent = sem_domin_effects_mnd_total%>%
  mutate(percent=abs(value)/sum(abs(value)),value=NULL)
sem_domin_effects_mnd.longer = sem_domin_effects_mnd %>% filter(type != 'total' 
                                                                #&value != 0
)

# Relative total effect
sem_domin_effects_mnd.percent.plot = 
  ggplot(data=sem_domin_effects_mnd.percent,aes(percent,predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_text(aes(label=paste(sprintf("%.1f",percent*100),"%",sep="")),hjust=-0.1,size=3)+
  facet_wrap(~"Relative total effect")+
  scale_y_discrete(limits=c("mnfunc_d", "mntd", "mnfd", "mnnd"),
                   labels=c("E-N MNFD", "E-N MNTD", "E-N MNRFD", "E-N MNND"),
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
sem_domin_effects_mnd.dir_indir.mpd.plot = 
  ggplot(data=sem_domin_effects_mnd.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects",
                                          "indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mnfunc_d", "mntd", "mnfd", "mnnd"),
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
sem_domin_mnd_space = ggplot()+
  facet_wrap(~"Dominance")+
  theme_test() + 
  theme(strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.015,1.2,0.2),units="lines"))

# Merge three plots
sem_domin_mnd_plot_all = ggarrange(sem_domin_mnd_space,
                                   sem_domin_effects_mnd.dir_indir.mpd.plot,
                                   sem_domin_effects_mnd.percent.plot,
                                   nrow = 1, ncol = 3,
                                   labels = c('(b)', '', ''),
                                   vjust = 1.8)


ggsave(plot = sem_domin_mnd_plot_all,
       "results/figures_ages1_35_top50_aposi/fitted_species_d/sem_domin_mnd_plot_all.svg",
       width=9,height=3)




#### SEM: Success ~ mND mRFD ~ mfd mpd for all species ####
### Original model ###
require(glmmTMB)
library(piecewiseSEM)
require(optimx)
require(dplyr)

### SEM for establishment
setwd("D:/R projects/BSS")
load('code/results_analyzing/analysing_ages1_35_top50_aposi_data/dat_suc_sp.rdata')

numcols = grep("^m.",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)

estab_sem1_md_all = psem(
  glmmTMB(estab~mnd+mfd+mpd_all+mfunc_d_all+
            (1|species)+
            (1|field)
          #+(1|f_p)
          , family=binomial, data=dat_suc_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd~mpd_all+
            mfunc_d_all+
            (1|species)+
            (1|field)
          +(1|f_p)
          , family=gaussian, data=dat_suc_sps),
  glmmTMB(mfd~mpd_all+
            mfunc_d_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  mfd %~~% mnd,
  #estab %~~% mpd,
  #estab %~~% mfunc_d,
  data=dat_suc_sps)
summary(estab_sem1_md_all)

estab_sem1_md.a_all = psem(
  glmmTMB(estab~mnd.a+mfd.a
          +mpd.a_all
          +mfunc_d.a_all
          #+(1|species)
          +(1|field)
          #+(1|f_p)
          , family=binomial, data=dat_suc_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a_all+
            mfunc_d.a_all+
            #(1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data = dat_suc_sps),
  glmmTMB(mfd.a~mpd.a_all+
            #mfunc_d.a_all+
            #(1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  mfd.a %~~% mnd.a,
  #estab %~~% mpd,
  #estab %~~% mfunc_d,
  data=dat_suc_sps)
summary(estab_sem1_md.a_all)

estab_sem2_md.a_all = psem(
  glmmTMB(estab~mnd.a+mfd.a
          +mpd.a_all
          +mfunc_d.a_all
          +(1|species)
          +(1|field)
          #+(1|f_p)
          , family=binomial, data=dat_suc_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~mpd.a_all+
            mfunc_d.a_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data = dat_suc_sps),
  glmmTMB(mfd.a~mpd.a_all+
            #mfunc_d.a_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  mfd.a %~~% mnd.a,
  #estab %~~% mpd,
  #estab %~~% mfunc_d,
  data=dat_suc_sps)
summary(estab_sem2_md.a_all)

estab_sem1_mnd_all = psem(
  glmmTMB(estab~mnnd+mnfd+
            mntd_all+
            mnfunc_d_all+
            (1|species)+
            (1|field)
          #+(1|f_p)
          , family=binomial, data=dat_suc_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))
  ),
  glmmTMB(mnnd~#mntd_all+
            mnfunc_d_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  glmmTMB(mnfd~mntd_all+
            mnfunc_d_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_suc_sps),
  mnfd %~~% mnnd,
  #estab %~~% mpd,
  #estab %~~% mfunc_d,
  data=dat_suc_sps)
summary(estab_sem1_mnd_all)


### SEM for Dominance ###
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c("establish", "dominant"))

domin_sem1_md_all = psem(
  glmmTMB(domin~#mnd+
            mfd+mpd_all+mfunc_d_all+
            #(1|species)+
            (1|field)
          #+(1|f_p)
          , family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd~#mpd_all+
            mfunc_d_all+
            #(1|species)+
            (1|field)
          #+(1|f_p)
          , family=gaussian, data=dat_dom_sps),
  glmmTMB(mfd~mpd_all+
            mfunc_d_all+
            #(1|species)+
            (1|field)
          #+(1|f_p)
          , family=gaussian, data=dat_dom_sps),
  mfd %~~% mnd,
  #domin %~~% mpd_all,
  #domin %~~% mfunc_d_all,
  data=dat_dom_sps)
summary(domin_sem1_md_all)

domin_sem1_md.a_all = psem(
  glmmTMB(domin~#mnd.a+
            mfd.a+
            mpd.a_all+
            mfunc_d.a_all+
            #(1|species)+
            (1|field)
          #+(1|f_p)
          , family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnd.a~#mpd.a_all+
            mfunc_d.a_all+
            #(1|species)+
            (1|field)
          #+(1|f_p)
          , family=gaussian, data = dat_dom_sps),
  glmmTMB(mfd.a~mpd.a_all+
            mfunc_d.a_all+
            #(1|species)+
            (1|field)
          +(1|f_p)
          , family=gaussian, data=dat_dom_sps),
  mfd.a %~~% mnd.a,
  #domin %~~% mpd,
  #domin %~~% mfunc_d,
  data=dat_dom_sps)
summary(domin_sem1_md.a_all)

domin_sem1_mnd_all = psem(
  glmmTMB(domin~mnnd+mnfd+
            mntd_all+
            mnfunc_d_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=binomial, data=dat_dom_sps,
          control=glmmTMBControl(optimizer=optim,
                                 optArgs=list(method="BFGS"))),
  glmmTMB(mnnd~mntd_all+
            #mnfunc_d_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_dom_sps),
  glmmTMB(mnfd~#mntd_all+
            mnfunc_d_all+
            (1|species)+
            (1|field)+
            (1|f_p), family=gaussian, data=dat_dom_sps),
  mnfd %~~% mnnd,
  #domin %~~% mpd,
  #domin %~~% mfunc_d,
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
    sd(dat_suc_sps$mfd))/sd.y_estab_sem1_md_all
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
    sd(dat_dom_sps$mfd))/sd.y_domin_sem1_md_all
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
  scale_y_discrete(limits=c("mfunc_d_all", "mpd_all", "mfd", "mnd"),
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
  scale_y_discrete(limits=c("mfunc_d_all", "mpd_all", "mfd", "mnd"),
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
       "results/figures_ages1_35_top50_aposi/all_species_d/sem_estab_md_all_plot_all.svg",
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
  scale_y_discrete(limits=c("mfunc_d_all", "mpd_all", "mfd", "mnd"),
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
  scale_y_discrete(limits=c("mfunc_d_all", "mpd_all", "mfd", "mnd"),
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
       "results/figures_ages1_35_top50_aposi/all_species_d/sem_domin_md_all_plot_all.svg",
       width=9,height=3)


##### Draw plot for md.a ##### 
### Establishment corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_estab_sem1_md.a_all = sqrt(predict(estab_sem1_md.a_all[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(estab_sem1_md.a_all[[1]])))+pi^2/3)
sd.y_estab_sem1_md.a_all
std.estab_sem1_md.a_all=unlist(fixef(estab_sem1_md.a_all[[1]]))[-1]*
  c(sd(dat_suc_sps$mpd.a_all),
    sd(dat_suc_sps$mfunc_d.a_all),
    sd(dat_suc_sps$mnd),
    sd(dat_suc_sps$mfd))/sd.y_estab_sem1_md.a_all
std.estab_sem1_md.a_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_estab_sem1_md.a_all=coefs(estab_sem1_md.a_all)# original coefficients of piecewiseSEM
coefs_estab_sem1_md.a_all[1:4,8]=std.estab_sem1_md.a_all

### Dominance corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_domin_sem1_md.a_all=sqrt(predict(domin_sem1_md.a_all[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(domin_sem1_md.a_all[[1]])))+pi^2/3)
sd.y_domin_sem1_md.a_all
std.domin_sem1_md.a_all=unlist(fixef(domin_sem1_md.a_all[[1]]))[-1]*
  c(sd(dat_dom_sps$mpd.a_all),
    sd(dat_dom_sps$mfunc_d.a_all),
    sd(dat_dom_sps$mnd),
    sd(dat_dom_sps$mfd))/sd.y_domin_sem1_md.a_all
std.domin_sem1_md.a_all

# merge the new corrected cofficients into the original cofficients data frame
coefs_domin_sem1_md.a_all=coefs(domin_sem1_md.a_all)# original coefficients of piecewiseSEM
coefs_domin_sem1_md.a_all[1:4,8]=std.domin_sem1_md.a_all
coefs_domin_sem1_md.a_all

### Calculate direct and indirect effects of SEM following Xu Meng (2022) GCB
### Estab
coefs_estab_sem1_md.a_all = coefs_estab_sem1_md.a_all[-nrow(coefs_estab_sem1_md.a_all),]
sem_estab_effects_md.a_all = data.frame(predictor = rep(c('mnd.a', 'mfd.a',
                                                          'mpd.a_all', 'mfunc_d.a_all'), 3),
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
sem_domin_effects_md.a_all = data.frame(predictor = rep(c('mnd.a', 'mfd.a',
                                                          'mpd.a_all', 'mfunc_d.a_all'), 3),
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
  scale_y_discrete(limits=c("mfunc_d.a_all", "mpd.a_all", "mfd.a", "mnd.a"),
                   labels=#c(bquote('I-N MFD'[ab]),bquote('I-N MPD'[ab]),bquote('I-N MRFD'[ab]),bquote('I-N MND'[ab]))
                     c('FD','PD','RFD','ND') ,
                   position="right")+
  scale_fill_stata()+# color palatte in ggthemes
  xlim(0,0.5)+
  xlab("")+
  #theme_test()+
  theme(axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.4,0.4,0.2,-0.2),units="lines"))+
  guides(fill="none")

# Direct and indirect effects
sem_estab_effects_md.a_all.dir_indir.mpd_all.plot = 
  ggplot(data=sem_estab_effects_md.a_all.longer,aes(x=value,y=predictor,fill=predictor))+
  geom_col(width=0.8)+
  geom_vline(xintercept = 0, linetype = 2)+
  facet_wrap(~type,labeller=as_labeller(c("direct"="Direct effects","indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mfunc_d.a_all", "mpd.a_all", "mfd.a", "mnd.a"),
                   position="right")+
  scale_x_continuous(breaks = c(seq(-0.4, 0.4, 0.2)))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_estab_md.a_all_space1 = ggplot(NULL)+theme_void()
sem_estab_md.a_all_space2 = ggplot(NULL)+theme_void()
# Merge three plots
sem_estab_md.a_all_plot_all = ggarrange(sem_estab_md.a_all_space1, sem_estab_md.a_all_space2,
                                        sem_estab_effects_md.a_all.dir_indir.mpd_all.plot,
                                        sem_estab_effects_md.a_all.percent.plot,
                                        nrow = 2, ncol = 2,
                                        heights = c(0.7, 0.3),
                                        labels = c('(a)', '', ''),
                                        vjust = 1.8)

emf('results/figures_ages1_35_top50_aposi/sem_estab_md.a_all_plot_all.emf',
    width=16, height=12, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
sem_estab_md.a_all_plot_all
dev.off() #turn off device and finalize file

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
  scale_y_discrete(limits=c("mfunc_d.a_all", "mpd.a_all", "mfd.a", "mnd.a"),
                   labels=#c(bquote('E-N MFD'[ab]),bquote('E-N MPD'[ab]),bquote('E-N MRFD'[ab]),bquote('E-N MND'[ab]))
                     c('FD','PD','RFD','ND'),
                   position="right")+
  scale_fill_stata()+# color palatte in ggthemes
  xlim(0,0.5)+
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
                                          "indirect"="Indirect effects")))+
  scale_y_discrete(limits=c("mfunc_d.a_all", "mpd.a_all", "mfd.a", "mnd.a"),
                   position="right")+
  scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2),
                     limits = c(-0.4, 0.4))+
  scale_fill_stata()+
  guides(fill="none")+
  #theme_test()+
  theme(axis.title=element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text=element_text(face="bold"))

### Space for SEM
sem_domin_md.a_all_space1 = ggplot(NULL)+theme_void()
sem_domin_md.a_all_space2 = ggplot(NULL)+theme_void()

# Merge three plots
sem_domin_md.a_all_plot_all = ggarrange(sem_domin_md.a_all_space1,sem_domin_md.a_all_space2,
                                        sem_domin_effects_md.a_all.dir_indir.mpd_all.plot,
                                        sem_domin_effects_md.a_all.percent.plot,
                                        nrow = 2, ncol = 2,
                                        heights = c(0.7, 0.3),
                                        labels = c('(b)', '', ''),
                                        vjust = 1.8)

emf('results/figures_ages1_35_top50_aposi/sem_domin_md.a_all_plot_all.emf',
    width=16, height=12, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
sem_domin_md.a_all_plot_all
dev.off() #turn off device and finalize file

#### Draw plot for mnd_all  
### Establishment corrected std coefficients in the glmer, following nakagawa 2013 MEE
sd.y_estab_sem1_mnd_all = sqrt(predict(estab_sem1_mnd_all[[1]],re.form=NA)%>%var+sum(unlist(VarCorr(estab_sem1_mnd_all[[1]])))+pi^2/3)
sd.y_estab_sem1_mnd_all
std.estab_sem1_mnd_all=unlist(fixef(estab_sem1_mnd_all[[1]]))[-1]*
  c(sd(dat_suc_sps$mntd_all),
    sd(dat_suc_sps$mnfunc_d_all),
    sd(dat_suc_sps$mnnd),
    sd(dat_suc_sps$mnfd))/sd.y_estab_sem1_mnd_all
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
    sd(dat_dom_sps$mnfd))/sd.y_domin_sem1_mnd_all
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
  scale_y_discrete(limits=c("mnfunc_d_all", "mntd_all", "mnfd", "mnnd"),
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
  scale_y_discrete(limits=c("mnfunc_d_all", "mntd_all", "mnfd", "mnnd"),
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
       "results/figures_ages1_35_top50_aposi/all_species_d/sem_estab_mnd_all_plot_all.svg",
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
  scale_y_discrete(limits=c("mnfunc_d_all", "mntd_all", "mnfd", "mnnd"),
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
  scale_y_discrete(limits=c("mnfunc_d_all", "mntd_all", "mnfd", "mnnd"),
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
       "results/figures_ages1_35_top50_aposi/all_species_d/sem_domin_mnd_all_plot_all.svg",
       width=9,height=3)

