############### Fast Start ####################
rm(list = ls())
setwd("~/BSS_coexist_v1")

### load packages
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(ggpubr)
library(devEMF)
library(ggeffects)
library(stringr)
library(viridis)
library(showtext)
library(MuMIn)
library(extrafont)
windowsFonts(Arial=windowsFont("Arial"))
library(scatterpie)
source('code/function/plot_func.R')

gre = "#55a868ff"
ora = "#dd8452ff"
yel = "#ccb974ff"
blu = "#4c72b0ff"
colors_4d = c(blu,ora,gre,yel)
colors_2d = c(ora, blu)

load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_allb_data/dat_suc_sp.rdata")
dat_suc_sp_bh_allb = dat_suc_sp
load("code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/dat_suc_sp.rdata")
dat_suc_sp_bh_partialb = dat_suc_sp
rm(dat_suc_sp)

dat_suc_sp_bh_allb$f_p_sp = paste(dat_suc_sp_bh_allb$f_p,
                                  dat_suc_sp_bh_allb$species,
                                  sep = '_')

dat_suc_sp_bh_partialb$f_p_sp = paste(dat_suc_sp_bh_partialb$f_p,
                                      dat_suc_sp_bh_partialb$species,
                                      sep = '_')

dat_suc_sp_comparison = dat_suc_sp_bh_partialb %>% 
  left_join(dat_suc_sp_bh_allb[,c('f_p_sp', "mnd.a" ,"mlgfd.a")],
            by = 'f_p_sp') %>% 
  filter(!is.na(mnd.a.y))

plot(dat_suc_sp_comparison$mnd.a.x,
     dat_suc_sp_comparison$mnd.a.y)



### the corrletion between estimated ND values from model_bh_partialb and from model_bh_allb
nd_model_cor = lmer(mnd.a.y ~ mnd.a.x + (1|f_p),
                    data = dat_suc_sp_comparison,
                    REML = T)
r.squaredGLMM(nd_model_cor)

nd_model_cor = lm(mnd.a.y ~ mnd.a.x,
                  data = dat_suc_sp_comparison)
summary(nd_model_cor)

#### Get the model predictive value
v = seq(min(dat_suc_sp_comparison$mnd.a.x),
        max(dat_suc_sp_comparison$mnd.a.x),
        length.out = 100)

nd_model_cor_pre = ggpredict(nd_model_cor,
                             type = 'fixed', ci_level = 0.99,
                             terms = 'mnd.a.x[v]')

range_nd_x = range(dat_suc_sp_comparison$mnd.a.x)
range_nd_y = range(dat_suc_sp_comparison$mnd.a.y)

#ggthemr::ggthemr(palette = "fresh", layout = "clean")
fig_nd_model_cor = ggplot()+
  geom_point(
    data = dat_suc_sp_comparison,
    mapping = aes(x = mnd.a.x, y = mnd.a.y),
    shape = 16,            # solid circle
    #color = "grey30",      # softer dark grey
    color = colors_4d[1],
    size = 1.5,            # smaller, subtle points
    alpha = 0.5
  )+
  geom_ribbon(data = nd_model_cor_pre,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = colors_4d[1], alpha = 0.5) + # CI ribbon 
  geom_line(data = nd_model_cor_pre,
            mapping = aes(x = x, y = predicted),
            linewidth = 1,
            color = "grey30")+
  labs(x = 'Partial exponential Beverton-Holt model',
       y = 'Exponential Beverton-Holt model')+
  annotate("text", x=range_nd_x[1]+(range_nd_x[2]-range_nd_x[1])*0.15,
           y=-0.4+(0.8--0.4)*0.99,
           label = expression(R^2 ~ "= 0.39"), 
           size = 5
  )+
  annotate("text", x=range_nd_x[1]+(range_nd_x[2]-range_nd_x[1])*0.15,
           y=-0.4+(0.8--0.4)*0.92,
           label = expression(italic(P) < 0.001), 
           size = 5
  )+
  ggtitle("Niche difference") + 
  scale_x_continuous(limits = c(range_nd_x[1], 0.82)) +
  scale_y_continuous(breaks = seq(-0.4, 0.8, 0.4), limits = c(range_nd_y[1], 0.82)) +
  theme_regular_2() +
  theme(plot.title = element_text(face = 'plain', hjust = 0.5, vjust = 0.8,
                                  size = 15))
fig_nd_model_cor

## the corrletion between estimated lgfd values from model_bh_partialb and from model_bh_allb
lgfd_model_cor = lmer(mlgfd.a.y ~ mlgfd.a.x + (1|f_p),
                      data = dat_suc_sp_comparison,
                      REML = T)
cor.test(dat_suc_sp_comparison$mlgfd.a.x,
         dat_suc_sp_comparison$mlgfd.a.y)
summary(lgfd_model_cor)

lgfd_model_cor = lm(mlgfd.a.y ~ mlgfd.a.x,
                    data = dat_suc_sp_comparison)
summary(lgfd_model_cor)

#### Get the model predictive value
v = seq(min(dat_suc_sp_comparison$mlgfd.a.x),
        max(dat_suc_sp_comparison$mlgfd.a.x),
        length.out = 100)
lgfd_model_cor_pre = ggpredict(lgfd_model_cor, ci_level = 0.99,
                               type = 'fixed', terms = 'mlgfd.a.x[v]')

range_lgfd_x = range(dat_suc_sp_comparison$mlgfd.a.x)
range_lgfd_y = range(dat_suc_sp_comparison$mlgfd.a.y)

#ggthemr::ggthemr(palette = "fresh", layout = "clean")
fig_lgfd_model_cor = ggplot()+
  geom_point(
    data = dat_suc_sp_comparison,
    mapping = aes(x = mlgfd.a.x, y = mlgfd.a.y),
    shape = 16,            # solid circle
    #color = "grey30",      # softer dark grey
    color = colors_4d[2],
    size = 1.5,            # smaller, subtle points
    alpha = 0.5
  )+
  geom_ribbon(data = lgfd_model_cor_pre,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = colors_4d[2], alpha = 0.2) + # CI ribbon 
  geom_line(data = lgfd_model_cor_pre,
            mapping = aes(x = x, y = predicted),
            linewidth = 1,
            color = "grey30")+
  annotate("text", x=range_lgfd_x[1]+(range_lgfd_x[2]-range_lgfd_x[1])*0.15,
           y=-2.1+(2.1--2.1)*0.99,
           #label = expression("Relative fitness difference"["RFD"[ab]])
           label = expression(R^2 ~ "= 0.72")
            , size = 5
  )+
  annotate("text", x=range_lgfd_x[1]+(range_lgfd_x[2]-range_lgfd_x[1])*0.15,
           y=-2.1+(2.1--2.1)*0.92,
           label = expression(italic(P) < 0.001), 
           size = 5
  )+
  labs(x = 'Partial exponential Beverton-Holt model', y = ' ')+
  ggtitle("Relative fitness difference") + 
  scale_x_continuous(limits = c(range_lgfd_x[1], 1.8), 
                     breaks = seq(-1.6, 1.6, 0.8)) +
  scale_y_continuous(limits = c(-2.1, 2.1), 
                     breaks = seq(-2, 2, 1)) +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.1) ) +
  theme_regular_2() +
  theme(plot.title = element_text(face = 'plain', hjust = 0.5, vjust = 0.8,
                                  size = 15))
fig_lgfd_model_cor

#### Merge and export Fig.S1 ####
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
Fig._nfd_model_cor = ggdraw() +
  # Top row
  draw_plot(fig_nd_model_cor, x = x1, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(fig_lgfd_model_cor, x = x2, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot_label(
    label = c("a", "b"),
    x = c(x1 + 0.01, x2 + 0.01),
    y = c(top_row_y + plot_height,  # similar to original 0.995
          top_row_y + plot_height),
    hjust = 0, vjust = 1.1, size = 18, color = 'black'
  )
Fig._nfd_model_cor

emf('results/figures/Fig.S2_mod_cor.emf',
    width = 30*0.75, height = 15*0.8, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig._nfd_model_cor
dev.off() #turn off device and finalize file
