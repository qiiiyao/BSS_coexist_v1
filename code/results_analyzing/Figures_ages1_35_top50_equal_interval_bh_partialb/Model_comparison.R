### Clear work place ###
rm(list = ls())

setwd("D:/R projects/BSS")

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
library(patchwork)
library(ggbreak)
windowsFonts(Arial=windowsFont("Arial"))
library(scatterpie)
#loadfonts(device = "win")
source('code/function/plot_func.R')

scientific_format = function(x) {
  sprintf("%.0e", x)  # Adjust "%.1e" for desired precision
}


gre = "#55a868ff"
ora = "#dd8452ff"
yel = "#ccb974ff"
blu = "#4c72b0ff"
colors_4d = c(blu,ora,gre,yel)
colors_2d = c(ora, blu)

load("results/fit_results/BSS_exclude_trees_raw/plot_ages1_35_top50_equal_interval_model_comparison/df_loo_all.Rdata")

colnames(df_loo_all)

df_loo_all$model = factor(df_loo_all$model,
                          levels = c('bh_partialb',
                                     'bh_allb',
                                     'bh',
                                     'null',
                                     'ricker',
                                     'ricker_logn'))
df_loo_all_summary = df_loo_all %>% group_by(model) %>% 
  summarise_at(vars('looic', 'se_looic', 'elpd_diff', 'se_diff'), c(mean))
colnames(df_loo_all_summary)[2:ncol(df_loo_all_summary)] = c(
  'looic_mean', 'se_looic_mean', 'elpd_diff_mean', 'se_diff_mean'
)

##### Draw plot

(mod_comprison_looic = ggplot()+
    geom_boxplot(data = df_loo_all,
                 mapping = aes(x = model,
                               y = looic
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    geom_violin(data = df_loo_all,
                 mapping = aes(x = model,
                               y = looic
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    geom_point(data = df_loo_all_summary,
                 mapping = aes(x = model,
                               y = looic_mean
                               #,color = model
                 ),
                 color = colors_2d[1], shape = 17, 
               size = 4
    )+
    labs(x = 'Model', y = 'LOOIC'
         , fontface = 'bold'
    )+
    scale_x_discrete(labels = c('Partial_EBH', 'EBH', 'BH', 
                                'Exp', 'Ricker', 'Log_Ricker')) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(50, 170))+
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular() 
)

(mod_comprison_selooic = ggplot()+
    geom_boxplot(data = df_loo_all,
                 mapping = aes(x = model,
                               y = se_looic
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    geom_violin(data = df_loo_all,
                mapping = aes(x = model,
                              y = se_looic
                              #,color = model
                ),
                #width = 1,
                alpha = 0#, linewidth = 2
    )+
    geom_point(data = df_loo_all_summary,
               mapping = aes(x = model,
                             y = se_looic_mean
                             #,color = model
               ),
               color = colors_2d[1], shape = 17, 
               size = 4
    )+
    labs(x = 'Model', y = 'ELPD'
         , fontface = 'bold'
    ) +
    #scale_x_discrete(labels = c('Partial_EBH', 'EBH', 'BH', 
    #                    'Exp', 'Ricker', 'Log_Ricker')) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(0, 100))+
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular() 
)


(mod_comprison_elpd_diff_1 = ggplot()+
    geom_boxplot(data = df_loo_all,
                 mapping = aes(x = model,
                               y = elpd_diff
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    geom_point(data = df_loo_all_summary,
               mapping = aes(x = model,
                             y = elpd_diff_mean
                             #,color = model
               ),
               color = colors_2d[1], shape = 17, 
               size = 4
    )+
    labs(x = NULL, y = NULL
         , fontface = 'bold'
    )+
    scale_x_discrete(labels = NULL) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(-25, 0))+
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular() 
)


(mod_comprison_elpd_diff_2 = ggplot()+
    geom_boxplot(data = df_loo_all,
                 mapping = aes(x = model,
                               y = elpd_diff
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    geom_point(data = df_loo_all_summary,
               mapping = aes(x = model,
                             y = elpd_diff_mean
                             #,color = model
               ),
               color = colors_2d[1], shape = 17, 
               size = 4
    )+
    labs(x = NULL, y = NULL
         , fontface = 'bold'
    )+
    scale_x_discrete(labels = NULL) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(-2.5e+06, -1.3e+05))+
    scale_y_continuous(labels = scientific_format) +
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular() 
)


(mod_comprison_elpd_diff_3 = ggplot()+
    geom_boxplot(data = df_loo_all,
                 mapping = aes(x = model,
                               y = elpd_diff
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    geom_point(data = df_loo_all_summary,
               mapping = aes(x = model,
                             y = elpd_diff_mean
                             #,color = model
               ),
               color = colors_2d[1], shape = 17, 
               size = 4
    )+
    labs(x = 'Model', y = NULL
         , fontface = 'bold'
    )+
    scale_x_discrete(labels = c('Partial_EBH', 'EBH', 'BH', 
                                'Exp', 'Ricker', 'Log_Ricker')) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(-5e+08, -2e+08))+
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular() + 
    theme(axis.text.x = element_text(angle = 60,color = '#000000', hjust = 1))
)

y_axis_title_elpd_diff = ggplot() + 
  labs(y = expression(Delta ~ "ELPD")) +
  theme_void() +
  theme(
    axis.title.y = element_text(angle = 90, size = 14, hjust = 0.5)
  )

# Combine plots with the specified layout
mod_comprison_elpd_diff = (y_axis_title_elpd_diff | (mod_comprison_elpd_diff_1 /
                                             mod_comprison_elpd_diff_2 /
                                             mod_comprison_elpd_diff_3) + 
  plot_layout(heights = c(0.5, 0.25, 0.25))) +
  plot_layout(widths = c(0.1, 2))

(mod_comprison_se_elpd_diff_3 = ggplot()+
    geom_boxplot(data = df_loo_all,
                 mapping = aes(x = model,
                               y = se_diff
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    #geom_violin(data = df_loo_all,
             #   mapping = aes(x = model,
             #                 y = se_diff
                              #,color = model
             #   ),
                #width = 1,
            #    alpha = 0#, linewidth = 2
   # )+
    geom_point(data = df_loo_all_summary,
               mapping = aes(x = model,
                             y = se_diff_mean
                             #,color = model
               ),
               color = colors_2d[2], shape = 17, 
               size = 4
    )+
    labs(x = 'Model', y = NULL
         , fontface = 'bold'
    ) +
    scale_x_discrete(labels = c('Partial_EBH', 'EBH', 'BH', 
                                'Exp', 'Ricker', 'Log_Ricker')) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(0, 60))+
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular() + 
    theme(axis.text.x = element_text(angle = 60,color = '#000000', hjust = 1))
)

summary((df_loo_all %>% filter(model == 'ricker_logn'))$se_diff)

(mod_comprison_se_elpd_diff_2 = ggplot()+
    geom_boxplot(data = df_loo_all,
                 mapping = aes(x = model,
                               y = se_diff
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    #geom_violin(data = df_loo_all,
    #   mapping = aes(x = model,
    #                 y = se_diff
    #,color = model
    #   ),
    #width = 1,
    #    alpha = 0#, linewidth = 2
    # )+
    geom_point(data = df_loo_all_summary,
               mapping = aes(x = model,
                             y = se_diff_mean
                             #,color = model
               ),
               color = colors_2d[2], shape = 17, 
               size = 4
    )+
    labs(x = NULL, y = NULL
         , fontface = 'bold'
    ) +
    scale_x_discrete(labels = NULL) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(3.4e+04, 5e+05))+
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular() 
)

(mod_comprison_se_elpd_diff_1 = ggplot()+
    geom_boxplot(data = df_loo_all,
                 mapping = aes(x = model,
                               y = se_diff
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    #geom_violin(data = df_loo_all,
    #   mapping = aes(x = model,
    #                 y = se_diff
    #,color = model
    #   ),
    #width = 1,
    #    alpha = 0#, linewidth = 2
    # )+
    geom_point(data = df_loo_all_summary,
               mapping = aes(x = model,
                             y = se_diff_mean
                             #,color = model
               ),
               color = colors_2d[2], shape = 17, 
               size = 4
    )+
    labs(x = NULL, y = NULL
         , fontface = 'bold'
    ) +
    scale_x_discrete(labels = NULL) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(1e+8, 2e+8))+
    scale_y_continuous(labels = scientific_format) +
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular() 
)

y_axis_title_se_elpd_diff = ggplot() + 
  labs(y = expression('SE ' ~ Delta ~ "ELPD")) +
  theme_void() +
  theme(
    axis.title.y = element_text(angle = 90, size = 14, hjust = 0.5)
  )

# Combine plots with the specified layout
mod_comprison_se_elpd_diff = (y_axis_title_se_elpd_diff | (mod_comprison_se_elpd_diff_1 /
                                                       mod_comprison_se_elpd_diff_2 /
                                                       mod_comprison_se_elpd_diff_3) + 
                             plot_layout(heights = c(0.5, 0.25, 0.25))) +
  plot_layout(widths = c(0.1, 2))





#Fig.S1
# Merge and export eight plots
library(devEMF)
Fig_model_comparion = ggarrange(mod_comprison_elpd_diff, mod_comprison_se_elpd_diff,
                              nrow = 1, ncol = 2,
                              labels = c('a',
                                         'b'),
                              
                              font.label = list(size = 18)
                              #,vjust = 1.8
)
emf('results/figures_ages1_35_top50_equal_interval_bh_partialb/Fig_model_comparion.emf',
    width = 20*0.8, height = 25*0.8, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig_model_comparion
dev.off() #turn off device and finalize file



