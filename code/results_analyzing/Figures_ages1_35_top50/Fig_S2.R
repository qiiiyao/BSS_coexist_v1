############### Fast Start ####################
rm(list = ls())
setwd("D:/R projects/BSS")
load("results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top40/inter_all_c_alltime.rdata")
inter_all_c_alltime_top40 = inter_all_c
load("results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top30/inter_all_c_alltime.rdata")
inter_all_c_alltime_top30 = inter_all_c
rm(inter_all_c)
inter_all_c_alltime_top30$select_std = 'top30'
inter_all_c_alltime_top40$select_std = 'top40'

inter_all_c_alltime = rbind(inter_all_c_alltime_top30, inter_all_c_alltime_top40)

rm(list = c('inter_all_c_alltime_top40', 'inter_all_c_alltime_top30'))

### load packages
library(dplyr)
library(betareg)
library(stringr)
library(data.table)
library(ape)
library(devEMF)
library(reticulate)
inter_all_c = inter_all_c_alltime

## Functions for plot 
theme_for_coe_plot = function(x){
  ggplot2::theme_test() + 
    theme(legend.position = c(0.11, 0.76),
          legend.background = element_blank(),
          legend.key = element_blank(),
          text = element_text(family = 'Arial'),
          axis.title = element_text(face="bold"),
          axis.title.x = element_text(margin = margin(t = 0.5, r = 0,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          axis.title.y = element_text(margin = margin(t = 0, r = 1,
                                                      b = 0, l = 0,
                                                      unit = "cm")),
          plot.margin = margin(t = 0.2, r = 0.2, b = 0.5, l = 0.2, 
                               unit = "cm"))
}

# get the abs value of lgfd
inter_all_c$ablgfd = inter_all_c$lgfd
inter_all_c$ablgfd[inter_all_c$ablgfd < 0] = inter_all_c$ablgfd[inter_all_c$ablgfd < 0]*-1


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

gre = "#55a868ff"
ora = "#dd8452ff"
yel = "#ccb974ff"
blu = "#4c72b0ff"
colors_4d = c(blu,ora,gre,yel)
colors_2d = c(ora, blu)


########## Stage split two parts: estab noestab domin nodomin ##########
#### just mean and sd of raw data ####
inter_all_forestab = inter_all_c %>%
  filter(stage_ij_estab %in% c("intro.estab_native",
                               "intro.noestab_native"))

length(unique(inter_all_forestab$sp_pair))
length(unique(inter_all_forestab$species_i))
length(unique(inter_all_forestab$species_j))

inter_all_forestab_msd = inter_all_forestab %>% 
  group_by(select_std,stage_ij_estab) %>% summarise_at(vars(nd, lgfd),
                                            c(mean), na.rm = T)
inter_all_forestab_outcome = inter_all_forestab %>% 
  group_by(select_std,stage_ij_estab) %>% summarise_at(vars(coexists,
                                                 priority,
                                                 i_exclusion_j,
                                                 j_exclusion_i),
                                            c(probs = 
                                                function(x){
                                                  y = length(x[which(x == 1)])/length(x)}))
inter_all_forestab_outcomes = cbind(inter_all_forestab_msd,
                                    inter_all_forestab_outcome[,3:ncol(inter_all_forestab_outcome)])

inter_all_forestab_outcomes$select_std = factor(inter_all_forestab_outcomes$select_std,
                                           levels = c('top30',
                                                      'top40'))
dat_text_est = data.frame(
  label = c("Establishment", " "),
  select_std = c("top30", "top40")
)

dat_text_est$select_std = factor(dat_text_est$select_std,
                                 levels = c('top30',
                                            'top40'))

BoundaryValues = data.frame(x = seq(0,0.99,l=100), 
                            y1 =  log10(1/ (1-seq(0,0.99,l=100))),
                            y2 = log10((1-seq(0,0.99,l=100))))      

BoundaryValues_P = data.frame(x = seq(0,-1,l=100), 
                              y1 =  log10( 1/ (1-seq(0,-1,l=100))),
                              y2 = log10((1-seq(0,-1,l=100))))      
OneLine =  data.frame(x= c(rev(BoundaryValues$x),BoundaryValues$x),
                      y= c(rev(BoundaryValues$y2),BoundaryValues$y1 ))
rect_data_1 = data.frame(
  rect_xmin = rep(-1,2),
  rect_xmax = rep(1,2),
  rect_ymin = rep(0,2),  # Different ymin and ymax for each facet
  rect_ymax = rep(1,2),
  select_std = c("top30", "top40")  # Faceting variable matching point_data
)

rect_data_2 = data.frame(
  rect_xmin = rep(-1,2),
  rect_xmax = rep(1,2),
  rect_ymin = rep(-1,2),  # Different ymin and ymax for each facet
  rect_ymax = rep(0,2),
  select_std = c("top30", "top40")  # Faceting variable matching point_data
)

#"#4c72b0ff" "#dd8452ff" "#55a868ff" "#ccb974ff"
Fig.S2_compare_estab_original_pie = ggplot()+
  geom_rect(data = rect_data_1,aes(xmin = rect_xmin, xmax = rect_xmax,
                                   ymin = rect_ymin, ymax = rect_ymax),
            fill = '#ccb974ff',
            alpha = 0.4)+
  geom_rect(data = rect_data_2,aes(xmin = rect_xmin, xmax = rect_xmax,
                                   ymin = rect_ymin, ymax = rect_ymax),
            fill = '#55a868ff',
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
                  pie_scale = 35,
                 cols= c('j_exclusion_i_probs', 'i_exclusion_j_probs',
                          'coexists_probs' ,'priority_probs'))+
  coord_fixed(xlim =c(-0.06, 0.49), ylim= c(-0.15, 0.40))+
  geom_point(aes(x = nd, y = lgfd, col = stage_ij_estab), size = 6,
             inter_all_forestab_outcomes)+ 
  geom_path(aes(nd, lgfd), arrow = arrow(length = unit(0.3, "cm"), ends = 'first'), size =1.5,
           data = inter_all_forestab_outcomes)+
  geom_text(data = dat_text_est,
            mapping = aes(x = 0.05, y = 0.4, label = label), fontface = 'bold'
  ) + 
  labs(x = 'Niche difference (ND)', 
       y = 'Relative fitness difference\n(RFD)') +
  scale_colour_manual(values = c(colors_4d[2],colors_4d[1]),
                      labels = c("intro.estab_native" = "Successful",
                                 "intro.noestab_native" = "Failed"),
                      name = ' ') + 
  scale_fill_manual(values = c('#55a868ff', '#ccb974ff', 'grey55','grey85'),
                    labels = c("i_exclusion_j_probs" = "Exotics Win",
                               "j_exclusion_i_probs" = "Natives Win",
                               "coexists_probs" = "Coexistence",
                               "priority_probs" = "Priority Effect"),
                    name = 'Outcome')+
  guides(
    fill = guide_legend(order = 1),  # Put shape legend first
    color = guide_legend(order = 2)   # Put color legend second
  ) +
  facet_wrap(~select_std,
             labeller=as_labeller(c("top30"="Top 15 exotics vs. Top 15 natives",
                                    "top40"="Top 20 exotics vs. Top 20 natives")),
             nrow = 2, ncol = 1)+
  theme_for_coe_plot()+
  theme(legend.position = c(0.45, 0.85),
        legend.direction = "vertical", legend.box = "horizontal",
        #legend.box.background = element_rect(colour = 'black',
        #                                   linewidth = 2),
        #legend.title = element_text(margin = margin(l = 0)),  # create more white space above the legend titles
        legend.box.margin = margin(l = 10)
  ) 
Fig.S2_compare_estab_original_pie

(Fig.S2_compare_estab_original = ggplot(NULL) +
    geom_pointrange(data = inter_all_forestab_msd,
                    mapping = aes(x = nd_mean,y = lgfd_mean,
                                  ymin = lgfd_mean-lgfd_sd,
                                  ymax = lgfd_mean+lgfd_sd,
                                  color = stage,
                                  shape = stage),
                    alpha = 1.2
                    ,fatten = 8,linewidth = 0.8
    )+
    geom_pointrange(data = inter_all_forestab_msd,
                    mapping = aes(x = nd_mean,y = lgfd_mean,
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
    scale_x_continuous(limits=c(-0.015, 0.75)) +
    scale_y_continuous(limits=c(-0.35, 0.35)) +
    geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.3) +
    labs(x = 'Niche difference (ND)', 
         y = 'Relative Fitness difference (RFD)') +
    stat_function(fun=spss.f, colour="black", linewidth = 0.3)+
    stat_function(fun=spss.ff, colour="black", linewidth = 0.3)+
    geom_text(data = dat_text_est,
              mapping = aes(x = 0.17, y = 0.35, label = label),
              fontface = 'bold'
    ) + 
    guides(colour=guide_legend(title=NULL,
                               override.aes = list(linewidth = 0.3),
                               family = 'Arial'),
           shape=guide_legend(title=NULL,
                              family = 'Arial')) +
    facet_wrap(~select_std,
               labeller=as_labeller(c("top30"="Top 15 exotics vs. Top 15 natives",
                                      "top40"="Top 20 exotics vs. Top 20 natives")),
               nrow = 2, ncol = 1)+
    theme_for_coe_plot() + 
    theme(legend.position = c(0.15, 0.6))
)


inter_all_fordomin = inter_all_c %>%
  filter(stage_ij_domin %in% c("estab.domin_native",
                               "estab.nodomin_native"))

length(unique(inter_all_fordomin$sp_pair))
length(unique(inter_all_fordomin$species_i))
length(unique(inter_all_fordomin$species_j))

inter_all_fordomin_msd = inter_all_fordomin %>% 
  group_by(select_std,stage_ij_domin) %>% summarise_at(vars(nd, lgfd),
                                                       c(mean), na.rm = T)
inter_all_fordomin_outcome = inter_all_fordomin %>% 
  group_by(select_std,stage_ij_domin) %>% summarise_at(vars(coexists,
                                                            priority,
                                                            i_exclusion_j,
                                                            j_exclusion_i),
                                                       c(probs = 
                                                           function(x){
                                                             y = length(x[which(x == 1)])/length(x)}))
inter_all_fordomin_outcomes = cbind(inter_all_fordomin_msd,
                                    inter_all_fordomin_outcome[,3:ncol(inter_all_fordomin_outcome)])

inter_all_fordomin_outcomes$select_std = factor(inter_all_fordomin_outcomes$select_std,
                                                levels = c('top30',
                                                           'top40'))
BoundaryValues = data.frame(x = seq(0,0.99,l=100), 
                            y1 =  log10(1/ (1-seq(0,0.99,l=100))),
                            y2 = log10((1-seq(0,0.99,l=100))))      

BoundaryValues_P = data.frame(x = seq(0,-1,l=100), 
                              y1 =  log10( 1/ (1-seq(0,-1,l=100))),
                              y2 = log10((1-seq(0,-1,l=100))))      
OneLine =  data.frame(x= c(rev(BoundaryValues$x),BoundaryValues$x),
                      y= c(rev(BoundaryValues$y2),BoundaryValues$y1 ))
rect_data_1 = data.frame(
  rect_xmin = rep(-1,2),
  rect_xmax = rep(1,2),
  rect_ymin = rep(0,2),  # Different ymin and ymax for each facet
  rect_ymax = rep(1,2),
  select_std = c("top30", "top40")  # Faceting variable matching point_data
)

rect_data_2 = data.frame(
  rect_xmin = rep(-1,2),
  rect_xmax = rep(1,2),
  rect_ymin = rep(-1,2),  # Different ymin and ymax for each facet
  rect_ymax = rep(0,2),
  select_std = c("top30", "top40")  # Faceting variable matching point_data
)

#"#4c72b0ff" "#dd8452ff" "#
dat_text_dom = data.frame(
  label = c("Dominance", " "),
  select_std = c("top30", "top40")
)
dat_text_dom$select_std = factor(dat_text_dom$select_std,
                                 levels = c('top30',
                                            'top40'))

Fig.S2_compare_domin_original_pie = ggplot()+
  geom_rect(data = rect_data_1,aes(xmin = rect_xmin, xmax = rect_xmax,
                                   ymin = rect_ymin, ymax = rect_ymax),
            fill = '#ccb974ff',
            alpha = 0.4)+
  geom_rect(data = rect_data_2,aes(xmin = rect_xmin, xmax = rect_xmax,
                                   ymin = rect_ymin, ymax = rect_ymax),
            fill = '#55a868ff',
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
                  pie_scale = 130,
                  cols= c('j_exclusion_i_probs', 'i_exclusion_j_probs',
                          'coexists_probs' ,'priority_probs'))+
  coord_fixed(xlim =c(-0.05, 0.5), ylim= c(-0.15, 0.4))+
  geom_point(aes(x = nd, y = lgfd, col = stage_ij_domin), size = 6,
             inter_all_fordomin_outcomes)+ 
  geom_path(aes(nd, lgfd), arrow = arrow(length = unit(0.3, "cm"), ends = 'first'), size =1.5,
            data = inter_all_fordomin_outcomes)+
  geom_text(data = dat_text_dom,
            mapping = aes(x = 0.04, y = 0.4, label = label), fontface = 'bold'
  ) + 
  labs(x = 'Niche difference (ND)', 
       y = '  ') +
  scale_colour_manual(values = c(colors_4d[2],colors_4d[1]),
                      labels = c("estab.domin_native" = "Successful",
                                 "estab.nodomin_native" = "Failed"),
                      name = ' ')+ 
  scale_fill_manual(values = c('#55a868ff', '#ccb974ff', 'grey55','grey85'),
                    labels = c("i_exclusion_j_probs" = "Exotics Win",
                               "j_exclusion_i_probs" = "Natives Win",
                               "coexists_probs" = "Coexistence",
                               "priority_probs" = "Priority Effect"),
                    name = 'Outcome')+
  guides(
    fill = guide_legend(order = 1),  # Put shape legend first
    color = guide_legend(order = 2)   # Put color legend second
  ) +
  facet_wrap(~select_std,
             labeller=as_labeller(c("top30"="Top 15 exotics vs. Top 15 natives",
                                    "top40"="Top 20 exotics vs. Top 20 natives")),
             nrow = 2, ncol = 1)+
  theme_for_coe_plot()+
  theme(legend.position = 'none',
        legend.direction = "vertical", legend.box = "horizontal"
        #legend.box.background = element_rect(colour = 'black',
        #                                   linewidth = 2),
        #legend.title = element_text(margin = margin(t = 140)),  # create more white space above the legend titles
        #legend.box.margin = margin(b = -80)
  ) 
Fig.S2_compare_domin_original_pie


(Fig.S2_compare_domin_original = ggplot(NULL) +
    geom_pointrange(data = inter_all_fordomin_msd,
                    mapping = aes(x = nd_mean,y = lgfd_mean,
                                  ymin = lgfd_mean-lgfd_sd,
                                  ymax = lgfd_mean+lgfd_sd,
                                  color = stage,
                                  shape = stage)
                    ,fatten = 8,linewidth = 0.8
    )+
    geom_pointrange(data = inter_all_fordomin_msd,
                    mapping = aes(x = nd_mean,y = lgfd_mean,
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
    scale_x_continuous(limits=c(-0.015, 0.75)) +
    scale_y_continuous(limits=c(-0.35, 0.35)) +
    geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.3) +
    labs(x = 'Niche difference (ND)', 
         y = '  ') +
    stat_function(fun=spss.f, colour="black", linewidth = 0.3)+
    stat_function(fun=spss.ff, colour="black", linewidth = 0.3)+
    geom_text(
      data  = dat_text_dom,
      mapping = aes(x = 0.14, y = 0.35, label = label),
      fontface = 'bold'
    ) + 
    guides(colour=guide_legend(title=NULL,
                               override.aes = list(linewidth = 0.3),
                               family = 'Arial'),
           shape=guide_legend(title=NULL,
                              family = 'Arial')) +
    facet_wrap(~select_std,
               labeller=as_labeller(c("top30"="Top 15 exotics vs. Top 15 natives",
                                      "top40"="Top 20 exotics vs. Top 20 natives")),
               nrow = 3, ncol = 1)+
    theme_for_coe_plot() + 
    theme(legend.position = "none")
)


###### Fig. S2 #########
library(cowplot)

Fig.S2 = ggdraw() +
  draw_plot(Fig.S2_compare_domin_original_pie, 0.5, 0, 0.47, .98) +
  draw_plot(Fig.S2_compare_estab_original_pie, 0.06, 0, 0.47, .98) +
  draw_plot_label(c("(a)", "(b)"),
                  c(0.098, 0.536),
                  c(0.97, 0.97)
                  ,size = 12
  )
Fig.S2

emf('results/figures_ages1_35_top50/Fig.S2.emf',
    width = 20, height = 17, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S2
dev.off() #turn off device and finalize file
