#### Fast start for analyzing mnd ~ mpd, mnd ~ mconti_func_d, mfd ~ mpd, mfd ~ mconti_func_d ####
rm(list = ls())
setwd("~/BSS_coexist_v1")

require(INLA)
require(phyr)
require(inlabru)
require(tibble)
require(lme4)
require(nlme)
require(lmerTest)
require(minpack.lm)

## Functions for plot 
source('code/function/plot_func.R')

gre = "#55a868ff"
ora = "#dd8452ff"
yel = "#ccb974ff"
blu = "#4c72b0ff"
colors_4d = c(blu,ora,gre,yel)
colors_2d = c(ora, blu)

load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/dat_suc_sp.rdata')

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




#### Relationships among abundance weighted mean differences ####
##### predictive curve for mnd.a ~ mpd.a_all #####
lincombs.data.mnd.a.mpd.a_all = data.frame(mpd.a_all=seq(min(dat_suc_sp$mpd.a_all),
                                                         max(dat_suc_sp$mpd.a_all),
                                                         length=100),
                                           mconti_func_d.a_all=mean(dat_suc_sp$mconti_func_d.a_all))

lincombs.matrix.mnd.a.mpd.a_all=model.matrix(~mpd.a_all,
                                             data=lincombs.data.mnd.a.mpd.a_all)
lincombs.matrix.mnd.a.mpd.a_all=as.data.frame(lincombs.matrix.mnd.a.mpd.a_all)
lincombs.mnd.a.mpd.a_all=inla.make.lincombs(lincombs.matrix.mnd.a.mpd.a_all)

inla.model_lincombs.mnd.a.mpd.a_all = pglmm(mnd.a ~ mpd.a_all +#(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.mnd.a.mpd.a_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.mnd.a.mpd.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.a.mpd.a_all$predicted.value=inla.model_lincombs.mnd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.a.mpd.a_all$lower=inla.model_lincombs.mnd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.a.mpd.a_all$upper=inla.model_lincombs.mnd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.a.mpd.a_all

save(lincombs.data.mnd.a.mpd.a_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.a.mpd.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.a.mpd.a_all.rdata')
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
(Fig_mnd.a_mpd.a_all = ggplot(data=lincombs.data.mnd.a.mpd.a_all, aes(x=mpd.a_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mpd.a_all, y=mnd.a),
               shape = 16,            # solid circle
               #color = "grey30",      # softer dark grey
               color = colors_4d[1],
               size = 1.5,            # smaller, subtle points
               alpha = 0.5,
               position=position_jitter(height=0.2))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1], alpha=0.2)+
    geom_line(linetype = 1, linewidth = 1, color = "grey30")+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(limits = c(-0.83, max(dat_suc_sp$mnd.a)+0.1),
                       breaks = seq(-0.8, 0.8, 0.4))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x = ' ',
         #y = expression('ND'[ab])
         y = 'Niche difference (ND)')+
    annotate(geom="text",
             x=min(dat_suc_sp$mpd.a_all)+(max(dat_suc_sp$mpd.a_all)-min(dat_suc_sp$mpd.a_all))*0.01,
             y=c(-0.58,-0.78),
             label=c("italic(β)['MPD'[ab]] == '-0.00002'",
                     "'95% CI' == '[-0.00004, -0.00001]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_2() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)


##### predictive curve for mnd.a ~ mconti_func_d.a_all #####
lincombs.data.mnd.a.mconti_func_d.a_all = data.frame(mconti_func_d.a_all=seq(min(dat_suc_sp$mconti_func_d.a_all),
                                                                 max(dat_suc_sp$mconti_func_d.a_all),
                                                                 length=100),
                                               mpd.a_all=mean(dat_suc_sp$mpd.a_all))

lincombs.matrix.mnd.a.mconti_func_d.a_all=model.matrix(~mconti_func_d.a_all,
                                                 data=lincombs.data.mnd.a.mconti_func_d.a_all)
lincombs.matrix.mnd.a.mconti_func_d.a_all=as.data.frame(lincombs.matrix.mnd.a.mconti_func_d.a_all)
lincombs.mnd.a.mconti_func_d.a_all=inla.make.lincombs(lincombs.matrix.mnd.a.mconti_func_d.a_all)

inla.model_lincombs.mnd.a.mconti_func_d.a_all = pglmm(mnd.a ~ mconti_func_d.a_all+#(1|species) + 
                                                  (1|f_p) + (1|field), data = dat_suc_sp,
                                                family = "gaussian", cov_ranef = list(species = tree),
                                                bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                            config = TRUE),
                                                                     quantiles=c(0.025,0.5,0.975),
                                                                     lincomb=lincombs.mnd.a.mconti_func_d.a_all,
                                                                     control.predictor=list(compute=T)),
                                                bayes = T)

inla.model_lincombs.mnd.a.mconti_func_d.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.a.mconti_func_d.a_all$predicted.value=inla.model_lincombs.mnd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.a.mconti_func_d.a_all$lower=inla.model_lincombs.mnd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.a.mconti_func_d.a_all$upper=inla.model_lincombs.mnd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.a.mconti_func_d.a_all

save(lincombs.data.mnd.a.mconti_func_d.a_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.a.mconti_func_d.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.a.mconti_func_d.a_all.rdata')
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
(Fig_mnd.a_mconti_func_d.a_all = ggplot(data=lincombs.data.mnd.a.mconti_func_d.a_all, aes(x=mconti_func_d.a_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d.a_all, y=mnd.a),
               shape = 16,            # solid circle
               #color = "grey30",      # softer dark grey
               color = colors_4d[1],
               size = 1.5,            # smaller, subtle points
               alpha = 0.5,
               position=position_jitter(height=0.2))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill= colors_4d[1],alpha=0.2)+
    geom_line(linetype = 1, linewidth = 1, color = "grey30")+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(limits = c(-0.83, max(dat_suc_sp$mnd.a)+0.1),
                       breaks = seq(-0.8, 0.8, 0.4))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x= ' ',
         y= ' ')+
    annotate(geom="text",
             x=min(dat_suc_sp$mconti_func_d.a_all)+
               (max(dat_suc_sp$mconti_func_d.a_all)-min(dat_suc_sp$mconti_func_d.a_all))*0.01,
             y=c(-0.58,-0.78),
             label=c("italic(β)['MFD'[ab]] == '-0.00223'",
                     "'95% CI' == '[-0.00429, -0.00017]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_2() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)
#mconti_func_d.a_all -0.00223   -0.00429   -0.00017


##### predictive curve for mlgfd.a ~ mpd.a_all #####
lincombs.data.mlgfd.a.mpd.a_all = data.frame(mpd.a_all=seq(min(dat_suc_sp$mpd.a_all),
                                                           max(dat_suc_sp$mpd.a_all),
                                                           length=100),
                                             mconti_func_d.a_all=mean(dat_suc_sp$mconti_func_d.a_all))

lincombs.matrix.mlgfd.a.mpd.a_all=model.matrix(~mpd.a_all,
                                               data=lincombs.data.mlgfd.a.mpd.a_all)
lincombs.matrix.mlgfd.a.mpd.a_all=as.data.frame(lincombs.matrix.mlgfd.a.mpd.a_all)
lincombs.mlgfd.a.mpd.a_all=inla.make.lincombs(lincombs.matrix.mlgfd.a.mpd.a_all)

inla.model_lincombs.mlgfd.a.mpd.a_all = pglmm(mlgfd.a ~ mpd.a_all+#(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "gaussian", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.mlgfd.a.mpd.a_all,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)


inla.model_lincombs.mlgfd.a.mpd.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.a.mpd.a_all$predicted.value=inla.model_lincombs.mlgfd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.a.mpd.a_all$lower=inla.model_lincombs.mlgfd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.a.mpd.a_all$upper=inla.model_lincombs.mlgfd.a.mpd.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.a.mpd.a_all

save(lincombs.data.mlgfd.a.mpd.a_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.a.mpd.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.a.mpd.a_all.rdata')
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
(Fig_mlgfd.a_mpd.a_all = ggplot(data=lincombs.data.mlgfd.a.mpd.a_all, aes(x=mpd.a_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mpd.a_all, y=mlgfd.a),
               shape = 16,            # solid circle
               #color = "grey30",      # softer dark grey
               color = colors_4d[2],
               size = 1.5,            # smaller, subtle points
               alpha = 0.5,
               position=position_jitter(height=0.2))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(linetype = 1, linewidth = 1, color = "grey30")+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(limits = c(min(dat_suc_sp$mlgfd.a), 1.7),
                       breaks = seq(-1.6, 1.6, 0.8))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(#x= expression('MPD'[ab]),
         #y=expression('RFD'[ab])
         x = 'Phylogenetic difference (PD)',
         y = 'Relative fitness difference (RFD)'
         )+
    annotate(geom="text",
             x=min(dat_suc_sp$mpd.a_all)+
               (max(dat_suc_sp$mpd.a_all)-min(dat_suc_sp$mpd.a_all))*0.01,
             y=c(-1.2,-1.6),
             label=c("italic(β)['MPD'[ab]] == '-0.00024'",
                     "'95% CI' == '[-0.00031, -0.00017]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_2() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)
# mpd.a_all   -0.00024   -0.00031   -0.00017


##### predictive curve for mlgfd.a ~ mconti_func_d.a_all #####
lincombs.data.mlgfd.a.mconti_func_d.a_all = data.frame(mconti_func_d.a_all=seq(
                                                                   min(dat_suc_sp$mconti_func_d.a_all),
                                                                   max(dat_suc_sp$mconti_func_d.a_all),
                                                                   length=100),
                                                 mpd.a_all=mean(dat_suc_sp$mpd.a_all))

lincombs.matrix.mlgfd.a.mconti_func_d.a_all=model.matrix(~mconti_func_d.a_all,
                                                   data=lincombs.data.mlgfd.a.mconti_func_d.a_all)
lincombs.matrix.mlgfd.a.mconti_func_d.a_all=as.data.frame(lincombs.matrix.mlgfd.a.mconti_func_d.a_all)
lincombs.mlgfd.a.mconti_func_d.a_all=inla.make.lincombs(lincombs.matrix.mlgfd.a.mconti_func_d.a_all)

inla.model_lincombs.mlgfd.a.mconti_func_d.a_all = pglmm(mlgfd.a ~ mconti_func_d.a_all+#(1|species) + 
                                                    (1|f_p) + (1|field), data = dat_suc_sp,
                                                  family = "gaussian", cov_ranef = list(species = tree),
                                                  bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                              config = TRUE),
                                                                       quantiles=c(0.025,0.5,0.975),
                                                                       lincomb=lincombs.mlgfd.a.mconti_func_d.a_all,
                                                                       control.predictor=list(compute=T)),
                                                  bayes = T)

inla.model_lincombs.mlgfd.a.mconti_func_d.a_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.a.mconti_func_d.a_all$predicted.value=inla.model_lincombs.mlgfd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.a.mconti_func_d.a_all$lower=inla.model_lincombs.mlgfd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.a.mconti_func_d.a_all$upper=inla.model_lincombs.mlgfd.a.mconti_func_d.a_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.a.mconti_func_d.a_all

save(lincombs.data.mlgfd.a.mconti_func_d.a_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.a.mconti_func_d.a_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.a.mconti_func_d.a_all.rdata')
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
(Fig_mlgfd.a_mconti_func_d.a_all = ggplot(data=lincombs.data.mlgfd.a.mconti_func_d.a_all, aes(x=mconti_func_d.a_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d.a_all, y=mlgfd.a),
               shape = 16,            # solid circle
               #color = "grey30",      # softer dark grey
               color = colors_4d[2],
               size = 1.5,            # smaller, subtle points
               alpha = 0.5,
               position=position_jitter(height=0.2))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(linetype = 1, linewidth = 1, color = "grey30")+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(limits = c(min(dat_suc_sp$mlgfd.a), 1.7),
                       breaks = seq(-1.6, 1.6, 0.8))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(#x= expression('MFD'[ab]),
         x = 'Functional difference (FD)',
         y= ' ')+
    annotate(geom="text",
             x=min(dat_suc_sp$mconti_func_d.a_all)+
               (max(dat_suc_sp$mconti_func_d.a_all)-min(dat_suc_sp$mconti_func_d.a_all))*0.01,
             y=c(-1.2,-1.6),
             label=c("italic(β)['MFD'[ab]] == '-0.02736'",
                     "'95% CI' == '[-0.03487, -0.01985]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_2() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)
#mconti_func_d.a_all -0.02736   -0.03487   -0.01985


##### merge plots for Fig. S10 #####
library(cowplot)

# Grid settings
gap_x = 1 / 100   # ≈ 0.033, same as original horizontal gap
gap_y = 0.001     # vertical gap between rows

# Compute new plot size to fill canvas minus gap
plot_width  = (1 - gap_x) / 2 
plot_height = (1 - gap_y-0.03) / 2

# Vertical positions
top_row_y    = gap_y + plot_height
bottom_row_y = 0

# Horizontal positions
x1 = 0
x2 = plot_width + gap_x

# Create 2x2 plot with same gaps and full coverage
Fig.pd_fd_nd_rfd.a = ggdraw() +
  # Top row
  draw_plot(Fig_mnd.a_mpd.a_all, x = x1, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(Fig_mnd.a_mconti_func_d.a_all, x = x2, y = top_row_y,
            width = plot_width, height = plot_height) +
  # Bottom row
  draw_plot(Fig_mlgfd.a_mpd.a_all, x = x1, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(Fig_mlgfd.a_mconti_func_d.a_all, x = x2, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  # Labels – maintain same relative offset from plots
  draw_plot_label(
    label = c("a", "b", "c", "d"),
    x = c(x1 + 0.01, x2 + 0.01, x1 + 0.01, x2 + 0.01),
    y = c(top_row_y + plot_height + 0.03,  # similar to original 0.995
          top_row_y + plot_height + 0.03,
          bottom_row_y + plot_height + 0.03,
          bottom_row_y + plot_height + 0.03),
    hjust = 0, vjust = 1.1, size = 18, color = 'black'
  )
Fig.pd_fd_nd_rfd.a

emf('results/figures/Fig.S10.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.pd_fd_nd_rfd.a
dev.off() #turn off device and finalize file


#### Relationships among mean differences ####
##### predictive curve for mnd ~ mpd_all #####
lincombs.data.mnd.mpd_all = data.frame(mpd_all=seq(min(dat_suc_sp$mpd_all),
                                                         max(dat_suc_sp$mpd_all),
                                                         length=100),
                                           mconti_func_d_all=mean(dat_suc_sp$mconti_func_d_all))

lincombs.matrix.mnd.mpd_all=model.matrix(~mpd_all,
                                             data=lincombs.data.mnd.mpd_all)
lincombs.matrix.mnd.mpd_all=as.data.frame(lincombs.matrix.mnd.mpd_all)
lincombs.mnd.mpd_all=inla.make.lincombs(lincombs.matrix.mnd.mpd_all)

inla.model_lincombs.mnd.mpd_all = pglmm(mnd ~ mpd_all +#(1|species) + 
                                              (1|f_p) + (1|field), data = dat_suc_sp,
                                            family = "gaussian", cov_ranef = list(species = tree),
                                            bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                        config = TRUE),
                                                                 quantiles=c(0.025,0.5,0.975),
                                                                 lincomb=lincombs.mnd.mpd_all,
                                                                 control.predictor=list(compute=T)),
                                            bayes = T)

inla.model_lincombs.mnd.mpd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(7)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.mpd_all$predicted.value=inla.model_lincombs.mnd.mpd_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.mpd_all$lower=inla.model_lincombs.mnd.mpd_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.mpd_all$upper=inla.model_lincombs.mnd.mpd_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.mpd_all

save(lincombs.data.mnd.mpd_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.mpd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.mpd_all.rdata')
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
(Fig_mnd_mpd_all = ggplot(data=lincombs.data.mnd.mpd_all, aes(x=mpd_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=mnd),
               shape = 16,            # solid circle
               #color = "grey30",      # softer dark grey
               color = colors_4d[1],
               size = 1.5,            # smaller, subtle points
               alpha = 0.5,
               position=position_jitter(height=0.2))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[1], alpha=0.2)+
    geom_line(linetype = 1, linewidth = 1, color = "grey30")+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(limits = c(-0.42, 0.42),
                       breaks = seq(-0.4, 0.4, 0.2))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x = ' ',
         y = 'Niche difference (ND)')+
    annotate(geom="text",
             x=min(dat_suc_sp$mpd_all)+(max(dat_suc_sp$mpd_all)-min(dat_suc_sp$mpd_all))*0.01,
             y=c(-0.32,-0.4),
             label=c("italic(β)['MPD'] == '-0.00002'",
                     "'95% CI' == '[-0.00003, > -0.00001]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_2() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)
#mpd_all     -0.00002   -0.00003    0.00000

##### predictive curve for mnd ~ mconti_func_d_all #####
lincombs.data.mnd.mconti_func_d_all = data.frame(mconti_func_d_all=seq(min(dat_suc_sp$mconti_func_d_all),
                                                                             max(dat_suc_sp$mconti_func_d_all),
                                                                             length=100),
                                                     mpd_all=mean(dat_suc_sp$mpd_all))

lincombs.matrix.mnd.mconti_func_d_all=model.matrix(~mconti_func_d_all,
                                                       data=lincombs.data.mnd.mconti_func_d_all)
lincombs.matrix.mnd.mconti_func_d_all=as.data.frame(lincombs.matrix.mnd.mconti_func_d_all)
lincombs.mnd.mconti_func_d_all=inla.make.lincombs(lincombs.matrix.mnd.mconti_func_d_all)

inla.model_lincombs.mnd.mconti_func_d_all = pglmm(mnd ~ mconti_func_d_all+#(1|species) + 
                                                        (1|f_p) + (1|field), data = dat_suc_sp,
                                                      family = "gaussian", cov_ranef = list(species = tree),
                                                      bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                  config = TRUE),
                                                                           quantiles=c(0.025,0.5,0.975),
                                                                           lincomb=lincombs.mnd.mconti_func_d_all,
                                                                           control.predictor=list(compute=T)),
                                                      bayes = T)

inla.model_lincombs.mnd.mconti_func_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mnd.mconti_func_d_all$predicted.value=inla.model_lincombs.mnd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mnd.mconti_func_d_all$lower=inla.model_lincombs.mnd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mnd.mconti_func_d_all$upper=inla.model_lincombs.mnd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mnd.mconti_func_d_all

save(lincombs.data.mnd.mconti_func_d_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.mconti_func_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mnd.mconti_func_d_all.rdata')
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
(Fig_mnd_mconti_func_d_all = ggplot(data=lincombs.data.mnd.mconti_func_d_all,
                                    aes(x=mconti_func_d_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d_all, y=mnd),
               shape = 16,            # solid circle
               #color = "grey30",      # softer dark grey
               color = colors_4d[1],
               size = 1.5,            # smaller, subtle points
               alpha = 0.5,
               position=position_jitter(height=0.2))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill= colors_4d[1],alpha=0.2)+
    geom_line(linetype = 1, linewidth = 1, color = "grey30")+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(limits = c(-0.42, 0.42),
                       breaks = seq(-0.4, 0.4, 0.2))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(x= ' ',
         y= ' ')+
    annotate(geom="text",
             x=min(dat_suc_sp$mconti_func_d_all)+
               (max(dat_suc_sp$mconti_func_d_all)-min(dat_suc_sp$mconti_func_d_all))*0.01,
             y=c(-0.32,-0.4),
             label=c("italic(β)['MFD'] == '0.00270'",
                     "'95% CI' == '[0.00108, 0.00432]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_2() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)
#mconti_func_d_all  0.00270    0.00108    0.00432

##### predictive curve for mlgfd ~ mpd_all #####
lincombs.data.mlgfd.mpd_all = data.frame(mpd_all=seq(min(dat_suc_sp$mpd_all),
                                                           max(dat_suc_sp$mpd_all),
                                                           length=100),
                                             mconti_func_d_all=mean(dat_suc_sp$mconti_func_d_all))

lincombs.matrix.mlgfd.mpd_all=model.matrix(~mpd_all,
                                               data=lincombs.data.mlgfd.mpd_all)
lincombs.matrix.mlgfd.mpd_all=as.data.frame(lincombs.matrix.mlgfd.mpd_all)
lincombs.mlgfd.mpd_all=inla.make.lincombs(lincombs.matrix.mlgfd.mpd_all)

inla.model_lincombs.mlgfd.mpd_all = pglmm(mlgfd ~ mpd_all+#(1|species) + 
                                                (1|f_p) + (1|field), data = dat_suc_sp,
                                              family = "gaussian", cov_ranef = list(species = tree),
                                              bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                          config = TRUE),
                                                                   quantiles=c(0.025,0.5,0.975),
                                                                   lincomb=lincombs.mlgfd.mpd_all,
                                                                   control.predictor=list(compute=T)),
                                              bayes = T)


inla.model_lincombs.mlgfd.mpd_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.mpd_all$predicted.value=inla.model_lincombs.mlgfd.mpd_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.mpd_all$lower=inla.model_lincombs.mlgfd.mpd_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.mpd_all$upper=inla.model_lincombs.mlgfd.mpd_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.mpd_all

save(lincombs.data.mlgfd.mpd_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.mpd_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.mpd_all.rdata')
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
range(dat_suc_sp$mlgfd)
(Fig_mlgfd_mpd_all = ggplot(data=lincombs.data.mlgfd.mpd_all, aes(x=mpd_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=mlgfd),
               shape = 16,            # solid circle
               #color = "grey30",      # softer dark grey
               color = colors_4d[2],
               size = 1.5,            # smaller, subtle points
               alpha = 0.5,
               position=position_jitter(height=0.2))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(linetype = 1, linewidth = 1, color = "grey30")+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(limits = c(-1.05, 1.51), breaks = seq(-1.0, 1.5, 0.5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(#x= expression('MPD'),
         #y=expression('RFD')
         x = 'Phylogenetic difference (PD)',
         y = 'Relative fitness difference (RFD)')+
    annotate(geom="text",
             x=min(dat_suc_sp$mpd_all)+
               (max(dat_suc_sp$mpd_all)-min(dat_suc_sp$mpd_all))*0.01,
             y=c(-0.76,-1),
             label=c("italic(β)['MPD'] == '-0.00026'",
                     "'95% CI' == '[-0.00036, -0.00017]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_2() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)
# mpd_all  -0.00026   -0.00036   -0.00017




##### predictive curve for mlgfd ~ mconti_func_d_all #####
lincombs.data.mlgfd.mconti_func_d_all = data.frame(mconti_func_d_all=seq(
  min(dat_suc_sp$mconti_func_d_all),
  max(dat_suc_sp$mconti_func_d_all),
  length=100),
  mpd_all=mean(dat_suc_sp$mpd_all))

lincombs.matrix.mlgfd.mconti_func_d_all=model.matrix(~mconti_func_d_all,
                                                         data=lincombs.data.mlgfd.mconti_func_d_all)
lincombs.matrix.mlgfd.mconti_func_d_all=as.data.frame(lincombs.matrix.mlgfd.mconti_func_d_all)
lincombs.mlgfd.mconti_func_d_all=inla.make.lincombs(lincombs.matrix.mlgfd.mconti_func_d_all)

inla.model_lincombs.mlgfd.mconti_func_d_all = pglmm(mlgfd ~ mconti_func_d_all+#(1|species) + 
                                                          (1|f_p) + (1|field), data = dat_suc_sp,
                                                        family = "gaussian", cov_ranef = list(species = tree),
                                                        bayes_options = list(control.compute = list(dic=T,waic=T,cpo=T,
                                                                                                    config = TRUE),
                                                                             quantiles=c(0.025,0.5,0.975),
                                                                             lincomb=lincombs.mlgfd.mconti_func_d_all,
                                                                             control.predictor=list(compute=T)),
                                                        bayes = T)

inla.model_lincombs.mlgfd.mconti_func_d_all$inla.model$summary.fixed[c(1,3,5)]%>%round(5)##Extracting effects and confidence intervals for prediction curves from raw data

lincombs.data.mlgfd.mconti_func_d_all$predicted.value=inla.model_lincombs.mlgfd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,2]
lincombs.data.mlgfd.mconti_func_d_all$lower=inla.model_lincombs.mlgfd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,4]
lincombs.data.mlgfd.mconti_func_d_all$upper=inla.model_lincombs.mlgfd.mconti_func_d_all$inla.model$summary.lincomb.derived[1:100,6]
#lincombs.data.mlgfd.mconti_func_d_all

save(lincombs.data.mlgfd.mconti_func_d_all, 
     file = 'code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.mconti_func_d_all.rdata')

### Draw plot
load('code/results_analyzing/analysing_ages1_35_top50_equal_interval_bh_partialb_data/lincombs.data.mlgfd.mconti_func_d_all.rdata')
#ggthemr::ggthemr(palette = "fresh", layout = "clean")
(Fig_mlgfd_mconti_func_d_all = ggplot(data=lincombs.data.mlgfd.mconti_func_d_all, aes(x=mconti_func_d_all, y=predicted.value))+
    geom_point(data=dat_suc_sp, aes(x=mconti_func_d_all, y=mlgfd),
               shape = 16,            # solid circle
               #color = "grey30",      # softer dark grey
               color = colors_4d[2],
               size = 1.5,            # smaller, subtle points
               alpha = 0.5,
               position=position_jitter(height=0.2))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=colors_4d[2],alpha=0.2)+
    geom_line(linetype = 3, linewidth = 1, color = "grey30")+
    scale_color_gradient(low = turbo(4)[3],
                         high = turbo(4)[4])+
    scale_y_continuous(limits = c(-1.05, 1.51), breaks = seq(-1.0, 1.5, 0.5))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    labs(#x= expression('MFD'),
         x = 'Functional difference (FD)',
         y= ' ')+
    annotate(geom="text",
             x=min(dat_suc_sp$mconti_func_d_all)+
               (max(dat_suc_sp$mconti_func_d_all)-min(dat_suc_sp$mconti_func_d_all))*0.01,
             y=c(-0.78,-1),
             label=c("italic(β)['MFD'] == '0.00050'",
                     "'95% CI' == '[-0.00833, 0.00933]'"),
             parse=T,size=5, hjust = 0)+
    #ggtitle('Establishment')+
    theme_regular_2() 
  #+ theme(plot.title = element_text(hjust = 0.5, face = 'plain', family = 'Arial'))
)
# mconti_func_d_all  0.00050   -0.00833    0.00933

##### merge plots for Fig. S11 #####
library(cowplot)

# Grid settings
gap_x = 1 / 100   # ≈ 0.033, same as original horizontal gap
gap_y = 0.001     # vertical gap between rows

# Compute new plot size to fill canvas minus gap
plot_width  = (1 - gap_x) / 2 
plot_height = (1 - gap_y-0.03) / 2

# Vertical positions
top_row_y    = gap_y + plot_height
bottom_row_y = 0

# Horizontal positions
x1 = 0
x2 = plot_width + gap_x

# Create 2x2 plot with same gaps and full coverage
Fig.pd_fd_nd_rfd = ggdraw() +
  # Top row
  draw_plot(Fig_mnd_mpd_all, x = x1, y = top_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(Fig_mnd_mconti_func_d_all, x = x2, y = top_row_y,
            width = plot_width, height = plot_height) +
  # Bottom row
  draw_plot(Fig_mlgfd_mpd_all, x = x1, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  draw_plot(Fig_mlgfd_mconti_func_d_all, x = x2, y = bottom_row_y,
            width = plot_width, height = plot_height) +
  # Labels – maintain same relative offset from plots
  draw_plot_label(
    label = c("a", "b", "c", "d"),
    x = c(x1 + 0.01, x2 + 0.01, x1 + 0.01, x2 + 0.01),
    y = c(top_row_y + plot_height + 0.03,  # similar to original 0.995
          top_row_y + plot_height + 0.03,
          bottom_row_y + plot_height + 0.03,
          bottom_row_y + plot_height + 0.03),
    hjust = 0, vjust = 1.1, size = 18, color = 'black'
  )
Fig.pd_fd_nd_rfd

emf('results/figures/Fig.S11.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.pd_fd_nd_rfd
dev.off() #turn off device and finalize file

