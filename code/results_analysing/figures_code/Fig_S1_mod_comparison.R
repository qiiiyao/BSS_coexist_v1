### Clear work place ###
rm(list = ls())

setwd("~/BSS_coexist_v1")

library(lme4)
library(lmerTest)
library(FSA)
library(multcompView)
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
pur = "#c44e52ff"
colors_4d = c(blu,ora,gre,yel)
colors_2d = c(ora, blu)
colors_5d = c(blu,ora,gre,yel,pur)

load("results/fit_results/plot_ages1_35_top50_equal_interval_model_comparison/df_loo_all.Rdata")

colnames(df_loo_all)

### filter the data from ricker_log model since it has too many outliers
df_loo_all_c = df_loo_all %>% filter(model != 'ricker_logn')

df_loo_all_c_l = split(df_loo_all_c, df_loo_all_c$model)
hist(df_loo_all_c_l$bh_partialb$elpd_diff)
hist(df_loo_all_c_l$bh_allb$elpd_diff)
hist(df_loo_all_c_l$bh$elpd_diff)
hist(df_loo_all_c_l$null$elpd_diff)
hist(df_loo_all_c_l$ricker$elpd_diff)

df_loo_all_c$model = factor(df_loo_all_c$model,
                          levels = c('bh_partialb',
                                     'bh_allb',
                                     'bh',
                                     'null',
                                     'ricker'))
df_loo_all_c_summary = df_loo_all_c %>% group_by(model) %>% 
  summarise_at(vars('looic', 'se_looic', 'elpd_diff', 'se_diff'), c(mean))

colnames(df_loo_all_c_summary)[2:ncol(df_loo_all_c_summary)] = c(
  'looic_mean', 'se_looic_mean', 'elpd_diff_mean', 'se_diff_mean'
)


y_axis_title_elpd_diff = ggplot() + 
  labs(y = expression(Delta ~ "ELPD")) +
  theme_void() +
  theme(
    axis.title.y = element_text(angle = 90, size = 14, hjust = 0.5)
  )

mean(df_loo_all_c %>% filter(model == 'null') %>% pull('se_diff'))

#### multi-comparison Tukey HSD for elpd_diff
# Fit ANOVA
anova_res = aov(elpd_diff ~ model, data = df_loo_all_c)
# Tukey HSD
elpd_diff_hsd_res = data.frame(Comparison = rownames(TukeyHSD(anova_res)$model),
                               TukeyHSD(anova_res)$model)

elpd_diff_comparisons = elpd_diff_hsd_res$Comparison
elpd_diff_p_adj = elpd_diff_hsd_res$p.adj

elpd_diff_groups_split = strsplit(elpd_diff_comparisons, "-")
mat = matrix(1, nrow=length(unique(df_loo_all_c$model)),
             ncol=length(unique(df_loo_all_c$model)))
rownames(mat) = colnames(mat) = levels(df_loo_all_c$model)

for (i in seq_along(elpd_diff_groups_split)) {
  g1 = elpd_diff_groups_split[[i]][1]
  g2 = elpd_diff_groups_split[[i]][2]
  mat[g1, g2] = elpd_diff_p_adj[i]
  mat[g2, g1] = elpd_diff_p_adj[i]
}

elpd_diff_letters = multcompLetters(mat, threshold = 0.05)
elpd_diff_letters = data.frame(model = names(elpd_diff_letters$Letters), 
                               signif_elpd_diff = elpd_diff_letters$Letters)
df_loo_all_c_summary = df_loo_all_c_summary %>% left_join(elpd_diff_letters,
                                                          by = 'model')

#### multi-comparison Tukey HSD for se_elpd_diff
# Fit ANOVA
anova_res = aov(se_diff ~ model, data = df_loo_all_c)
# Tukey HSD
se_diff_hsd_res = data.frame(Comparison = rownames(TukeyHSD(anova_res)$model),
                               TukeyHSD(anova_res)$model)

se_diff_comparisons = se_diff_hsd_res$Comparison
se_diff_p_adj = se_diff_hsd_res$p.adj

se_diff_groups_split = strsplit(se_diff_comparisons, "-")
mat = matrix(1, nrow=length(unique(df_loo_all_c$model)),
             ncol=length(unique(df_loo_all_c$model)))
rownames(mat) = colnames(mat) = levels(df_loo_all_c$model)

for (i in seq_along(se_diff_groups_split)) {
  g1 = se_diff_groups_split[[i]][1]
  g2 = se_diff_groups_split[[i]][2]
  mat[g1, g2] = se_diff_p_adj[i]
  mat[g2, g1] = se_diff_p_adj[i]
}

se_diff_letters = multcompLetters(mat, threshold = 0.05, reversed = T)
se_diff_letters = data.frame(model = names(se_diff_letters$Letters), 
                               signif_se_diff = se_diff_letters$Letters)
df_loo_all_c_summary = df_loo_all_c_summary %>% left_join(se_diff_letters,
                                                          by = 'model')




#### multi-comparison for elpd_diff
kruskal.test(elpd_diff ~ model, data = df_loo_all_c)
elpd_diff_dunn_res = dunnTest(elpd_diff ~ model,
                              data = df_loo_all_c,
                              method = "bh")

elpd_diff_pvals = elpd_diff_dunn_res$res
elpd_diff_comparisons = elpd_diff_pvals$Comparison
elpd_diff_p_adj = elpd_diff_pvals$P.adj

elpd_diff_groups_split = strsplit(elpd_diff_comparisons, " - ")
mat = matrix(1, nrow=length(unique(df_loo_all_c$model)),
             ncol=length(unique(df_loo_all_c$model)))
rownames(mat) = colnames(mat) = levels(df_loo_all_c$model)

for (i in seq_along(elpd_diff_groups_split)) {
  g1 = elpd_diff_groups_split[[i]][1]
  g2 = elpd_diff_groups_split[[i]][2]
  mat[g1, g2] = elpd_diff_p_adj[i]
  mat[g2, g1] = elpd_diff_p_adj[i]
}

elpd_diff_letters = multcompLetters(mat, threshold = 0.05)
elpd_diff_letters = data.frame(model = names(elpd_diff_letters$Letters), 
                               signif_elpd_diff = elpd_diff_letters$Letters)


#### multi-comparison for se_elpd_diff
kruskal.test(se_diff ~ model, data = df_loo_all_c)
se_diff_dunn_res = dunnTest(se_diff ~ model,
                            data = df_loo_all_c,
                            method = "s")

se_diff_pvals = se_diff_dunn_res$res
se_diff_comparisons = se_diff_pvals$Comparison
se_diff_p_adj = se_diff_pvals$P.adj

se_diff_groups_split = strsplit(se_diff_comparisons, " - ")
mat = matrix(1, nrow=length(unique(df_loo_all_c$model)),
             ncol=length(unique(df_loo_all_c$model)))
rownames(mat) = colnames(mat) = levels(df_loo_all_c$model)

for (i in seq_along(se_diff_groups_split)) {
  g1 = se_diff_groups_split[[i]][1]
  g2 = se_diff_groups_split[[i]][2]
  mat[g1, g2] = se_diff_p_adj[i]
  mat[g2, g1] = se_diff_p_adj[i]
}
se_diff_letters = multcompLetters(mat, threshold = 0.05)
se_diff_letters$Letters
se_diff_letters = data.frame(model = names(se_diff_letters$Letters), 
                             signif_se_diff = se_diff_letters$Letters)



#### multi-comparison for looic
kruskal.test(looic ~ model, data = df_loo_all_c)
looic_dunn_res = dunnTest(looic ~ model,
                            data = df_loo_all_c,
                            method = "bonferroni")

looic_pvals = looic_dunn_res$res
looic_comparisons = looic_pvals$Comparison
looic_p_adj = looic_pvals$P.adj

looic_groups_split = strsplit(looic_comparisons, " - ")
mat = matrix(1, nrow=length(unique(df_loo_all_c$model)),
             ncol=length(unique(df_loo_all_c$model)))
rownames(mat) = colnames(mat) = levels(df_loo_all_c$model)

for (i in seq_along(looic_groups_split)) {
  g1 = looic_groups_split[[i]][1]
  g2 = looic_groups_split[[i]][2]
  mat[g1, g2] = looic_p_adj[i]
  mat[g2, g1] = looic_p_adj[i]
}

looic_letters = multcompLetters(mat, threshold = 0.05)
looic_letters$Letters


#### multi-comparison for se_looic
kruskal.test(se_looic ~ model, data = df_loo_all_c)
se_looic_dunn_res = dunnTest(se_looic ~ model,
                          data = df_loo_all_c,
                          method = "bonferroni")

se_looic_pvals = se_looic_dunn_res$res
se_looic_comparisons = se_looic_pvals$Comparison
se_looic_p_adj = se_looic_pvals$P.adj

se_looic_groups_split = strsplit(se_looic_comparisons, " - ")
mat = matrix(1, nrow=length(unique(df_loo_all_c$model)),
             ncol=length(unique(df_loo_all_c$model)))
rownames(mat) = colnames(mat) = levels(df_loo_all_c$model)

for (i in seq_along(se_looic_groups_split)) {
  g1 = se_looic_groups_split[[i]][1]
  g2 = se_looic_groups_split[[i]][2]
  mat[g1, g2] = se_looic_p_adj[i]
  mat[g2, g1] = se_looic_p_adj[i]
}

se_looic_letters = multcompLetters(mat, threshold = 0.05)
se_looic_letters$Letters


##### Draw plot

(mod_comprison_looic = ggplot()+
    geom_boxplot(data = df_loo_all_c,
                 mapping = aes(x = model,
                               y = looic
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    geom_violin(data = df_loo_all_c,
                mapping = aes(x = model,
                              y = looic
                              #,color = model
                ),
                #width = 1,
                alpha = 0#, linewidth = 2
    )+
    geom_point(data = df_loo_all_c_summary,
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
    geom_boxplot(data = df_loo_all_c,
                 mapping = aes(x = model,
                               y = se_looic
                               #,color = model
                 ),
                 #width = 1,
                 alpha = 0#, linewidth = 2
    )+
    geom_violin(data = df_loo_all_c,
                mapping = aes(x = model,
                              y = se_looic
                              #,color = model
                ),
                #width = 1,
                alpha = 0#, linewidth = 2
    )+
    geom_point(data = df_loo_all_c_summary,
               mapping = aes(x = model,
                             y = se_looic_mean
                             #,color = model
               ),
               color = colors_2d[1], shape = 17, 
               size = 4
    )+
    labs(x = 'Model', y = 'se_looic'
         , fontface = 'bold'
    ) +
    #scale_x_discrete(labels = c('Partial_EBH', 'EBH', 'BH', 
    #                    'Exp', 'Ricker', 'Log_Ricker')) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(0, 100))+
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular() 
)


ggthemr::ggthemr(palette = "fresh", layout = "clean")
(mod_comprison_elpd_diff = ggplot()+
    #geom_jitter(data = df_loo_all_c,
                #mapping = aes(x = model,
                          #    y = elpd_diff,
                          #    color = model
                #),#
                #width = 0.2,
                #color = 'grey',
                #fill = 'grey',
                #alpha = 0.1)+
    geom_boxplot(data = df_loo_all_c,
                 mapping = aes(x = model,
                               y = elpd_diff,
                               color = model
                 ),
                 #width = 1,
                 alpha = 0, linewidth = 0.5
    )+
    geom_text(data = df_loo_all_c_summary,
              aes(x = model,
                  y = max(df_loo_all_c$elpd_diff) + 3,
                  label = signif_elpd_diff,
                  color = model),
              size = 5) +
    geom_point(data = df_loo_all_c_summary,
               mapping = aes(x = model,
                             y = elpd_diff_mean
                             ,color = model
               ),
               #color = 'grey',
               shape = 17, 
               size = 4
    )+
    
    labs(x = 'Model', y = expression(Delta ~ "ELPD")
         , fontface = 'bold'
    )+
    scale_x_discrete(labels = c('Partial_EBH', 'EBH', 'BH', 
                                         'Exp', 'Ricker', 'Log_Ricker')) +
    scale_color_manual(values = colors_5d) + 
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(-25, max(df_loo_all_c$elpd_diff) + 3))+
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular_2() + 
    theme(axis.text.x = element_text(angle = 60,color = '#000000', hjust = 1))
)

ggthemr::ggthemr(palette = "fresh", layout = "clean")
(mod_comprison_se_elpd_diff = ggplot()+
    #geom_jitter(data = df_loo_all_c,
               #        mapping = aes(x = model,
               #                      y = se_diff,
                     #                color = model
                   # ),
                #width = 0.2,
                    #color = 'grey',
                    #fill = 'grey',
                  #  alpha = 0.1,)+
    geom_text(data = df_loo_all_c_summary,
              aes(x = model,
                  y = max(df_loo_all_c$se_diff) + 3,
                  label = signif_se_diff,
                  color = model),
              size = 5) +
    geom_boxplot(data = df_loo_all_c,
                 mapping = aes(x = model,
                               y = se_diff,
                               color = model
                 ),
                 #width = 1,
                 alpha = 0, linewidth = 0.5
    )+
    #geom_violin(data = df_loo_all_c,
      # mapping = aes(x = model,
      #               y = se_diff
      # ),
   # width = 0.5,
     #   alpha = 0#, linewidth = 2
    # )+
    geom_point(data = df_loo_all_c_summary,
               mapping = aes(x = model,
                             y = se_diff_mean
                             ,color = model
               ),
               #color = 'grey',
               shape = 17, 
               size = 4
    )+
    labs(x = 'Model', y = expression('SE ' ~ Delta ~ "ELPD")
         , fontface = 'bold'
    ) +
    scale_color_manual(values = colors_5d) + 
    scale_x_discrete(labels = c('Partial_EBH', 'EBH', 'BH', 
                                'Exp', 'Ricker', 'Log_Ricker')) +
    coord_cartesian(#xlim=c(-0.06, 0.1),
      ylim=c(0, max(df_loo_all_c$se_diff) + 3))+
    theme(axis.title.y = element_text(face = "italic"))+
    theme_regular_2() +
    theme(axis.text.x = element_text(angle = 60,color = '#000000', hjust = 1))
)



#Fig.S1
# Merge and export eight plots
library(devEMF)
library(cowplot)

# Grid settings
gap_x = 1 / 100   # ≈ 0.033, same as original horizontal gap
gap_y = 0.001     # vertical gap between rows

# Compute new plot size to fill canvas minus gap
plot_width  = (1 - gap_x) / 2 
plot_height = (1 - gap_y-0.03)

# Vertical positions
top_row_y    = gap_y

# Horizontal positions
x1 = 0
x2 = plot_width + gap_x

# Create 2x2 plot with same gaps and full coverage
Fig_model_comparion = ggdraw() +
  # Top row
  draw_plot(mod_comprison_elpd_diff, x = x1, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(mod_comprison_se_elpd_diff, x = x2, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b"),
    x = c(x1 + 0.01, x2 + 0.01),
    y = c(top_row_y + plot_height + 0.03,  # similar to original 0.995
          top_row_y + plot_height + 0.03),
    hjust = 0, vjust = 1.1, size = 18, color = 'black'
  )
Fig_model_comparion

emf('results/figures/Fig_model_comparion.emf',
    width = 30*0.65, height = 15*0.8, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig_model_comparion
dev.off() #turn off device and finalize file



