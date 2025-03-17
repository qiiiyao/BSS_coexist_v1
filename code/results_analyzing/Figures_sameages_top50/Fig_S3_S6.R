############ Fast start for analyzing invasion success probability ~ mnd+mfd+mpd+mfunc_d for all species #########
rm(list = ls())
setwd("D:/R projects/BSS")

library(phyr)
library(tibble)
library(lme4)
require(ape)
library(scales)
library(ggthemes)

load('code/results_analyzing/analysing_sameages_top50_data/dat_suc_sp.rdata')
numcols = grep("^m",names(dat_suc_sp))
dat_suc_sps = dat_suc_sp
dat_suc_sps[,numcols] = scale(dat_suc_sps[,numcols])

## Functions for plot 
theme_regular = function(x){
  ggplot2::theme_test() + 
    theme(text = element_text(family = 'Arial'),
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


#### estab ####
pc_prior = list(prec=list("pc.prec", param=c(0.1,0.01)))
dat_suc_sps$species_1 = as.factor(dat_suc_sps$species)
estab_sp_names = unique(dat_suc_sps$species)

tree = read.tree('data/original data/phylo_tree332.txt')
estab_tree_fit = keep.tip(tree, estab_sp_names)
estab_vcv_tree = ape::vcv(estab_tree_fit, model = "Brownian", corr = FALSE)
estab_vcv_tree_sparse = inla.as.sparse(solve(estab_vcv_tree))


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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment~mnd
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mnd.rdata")

(estab.mnd.partial.logistic=ggplot(data=lincombs.data.estab.mnd,
                                   aes(x=mnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[3],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[3],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnd, y=estab),
               color = stata_pal("s2color")(4)[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x=' ',
         y="Establishment probability")+
    annotate(geom="text",x=c((min(dat_suc_sp$mnd) + max(dat_suc_sp$mnd))/2,
                             (min(dat_suc_sp$mnd) + max(dat_suc_sp$mnd))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MND'] == -2.88",
                     "'95%CI' == '[-3.71, -2.05]'"),parse=T,size=3.5)+
    theme_regular()
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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mlgfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mlgfd
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mlgfd.rdata")
(estab.mlgfd.partial.logistic=ggplot(data=lincombs.data.estab.mlgfd,aes(x=mlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[2],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[2],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mlgfd, y=estab),
               color = stata_pal("s2color")(4)[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=c(c((min(dat_suc_sp$mlgfd) + max(dat_suc_sp$mlgfd))/2,
                               (min(dat_suc_sp$mlgfd) + max(dat_suc_sp$mlgfd))/2)),
             y=c(0.80,0.70),
             label=c("italic(β)['MRFD'] == 2.42",
                     "'95%CI' == '[1.54, 3.30]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mpd_all, y=estab),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",x=c(c((min(dat_suc_sp$mpd_all) + max(dat_suc_sp$mpd_all))/2,
                               (min(dat_suc_sp$mpd_all) + max(dat_suc_sp$mpd_all))/2)),
             y=c(0.80,0.70),
             label=c("italic(β)['MPD'] == '-0.02'",
                     "'95%CI' == '[-0.03, 0.00]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mfunc_d_all, y=estab),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=c((min(dat_suc_sp$mfunc_d_all) + max(dat_suc_sp$mfunc_d_all))/2,
                               (min(dat_suc_sp$mfunc_d_all) + max(dat_suc_sp$mfunc_d_all))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MFD'] == 8.45",
                     "'95%CI' == '[3.93, 12.94]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
inla.model_lincombs.estab.mnnd$inla.model$summary.fixed[c(1,3,5)]%>%round(2)##Extracting effects and confidence intervals for prediction curves from raw data

#inla.model1_lincombs.estab.mnnd$summary.lincomb,This is logit scale, can not be used directly, 
# you need to use the posterior distribution to the original scale, 
# and then calculate the mean and confidence. The following method:
lincombs.data.estab.mnnd$predicted.value=unlist(lapply(lincombs.posterior.estab.mnnd,
                                                       function(x)inla.emarginal(fun=plogis,x)))
lincombs.data.estab.mnnd$lower=unlist(lapply(lincombs.posterior.estab.mnnd,function(x)inla.qmarginal(0.025,inla.tmarginal(fun=plogis,x))))
lincombs.data.estab.mnnd$upper=unlist(lapply(lincombs.posterior.estab.mnnd,function(x)inla.qmarginal(0.975,inla.tmarginal(fun=plogis,x))))

save(lincombs.data.estab.mnnd, 
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mnnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment~mnnd
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mnnd.rdata")

(estab.mnnd.partial.logistic=ggplot(data=lincombs.data.estab.mnnd,
                                    aes(x=mnnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[3],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[3],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnnd, y=estab),
               color = stata_pal("s2color")(4)[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x=' ',
         y="Establishment probability")+
    annotate(geom="text",x=c((min(dat_suc_sp$mnnd) + max(dat_suc_sp$mnnd))/2,
                             (min(dat_suc_sp$mnnd) + max(dat_suc_sp$mnnd))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MNND'] == -0.80",
                     "'95%CI' == '[-1.08, -1.53]'"),parse=T,size=3.5) + 
    theme_regular()
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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mnlgfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for establishment ~ mnlgfd
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.estab.mnlgfd.rdata")
(estab.mnlgfd.partial.logistic=ggplot(data=lincombs.data.estab.mnlgfd,aes(x=mnlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[2],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[2],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnlgfd, y=estab),
               color = stata_pal("s2color")(4)[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=c((min(dat_suc_sp$mnlgfd) + max(dat_suc_sp$mnlgfd))/2,
                             (min(dat_suc_sp$mnlgfd) + max(dat_suc_sp$mnlgfd))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MNRFD'] == 1.84",
                     "'95%CI' == '[1.34, 2.35]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
(estab.mntd_all.partial.logistic=ggplot(data=lincombs.data.estab.mntd_all,
                                        aes(x=mntd_all, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[4],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[4],size=1,
              linetype = 3)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mntd_all, y=estab),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='Establishment probability')+
    annotate(geom="text",x=c((min(dat_suc_sp$mntd_all) + max(dat_suc_sp$mntd_all))/2,
                             (min(dat_suc_sp$mntd_all) + max(dat_suc_sp$mntd_all))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MNTD'] == '-0.00'",
                     "'95%CI' == '[-0.00, 0.00]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_suc_sp, aes(x=mnfunc_d_all, y=estab),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x=' ', y='  ')+
    annotate(geom="text",x=c((min(dat_suc_sp$mnfunc_d_all) + max(dat_suc_sp$mnfunc_d_all))/2,
                             (min(dat_suc_sp$mnfunc_d_all) + max(dat_suc_sp$mnfunc_d_all))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MNFD'] == 0.97",
                     "'95%CI' == '[-1.24, 3.19]'"),
             parse=T,size=3.5) + 
    theme_regular()
)



#### domin ####
dat_dom_sp = dat_suc_sp %>% filter(stage %in% c('establish', 'dominant'))
dat_dom_sps = dat_suc_sps %>% filter(stage %in% c('establish', 'dominant'))

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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance~mnd
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mnd.rdata")
(domin.mnd.partial.logistic=ggplot(data=lincombs.data.domin.mnd,
                                   aes(x=mnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[3],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[3],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnd, y=domin),
               color = stata_pal("s2color")(4)[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='MND',
         y="Dominance probability")+
    annotate(geom="text",x=c((min(dat_dom_sp$mnd) + max(dat_dom_sp$mnd))/2,
                             (min(dat_dom_sp$mnd) + max(dat_dom_sp$mnd))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MND'] == -2.42",
                     "'95%CI' == '[-3.97, -0.87]'"),parse=T,size=3.5) + 
    theme_regular()
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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mlgfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mlgfd
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mlgfd.rdata")
(domin.mlgfd.partial.logistic=ggplot(data=lincombs.data.domin.mlgfd,aes(x=mlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[2],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[2],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mlgfd, y=domin),
               color = stata_pal("s2color")(4)[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='MRFD', y='  ')+
    annotate(geom="text",x=c((min(dat_dom_sp$mlgfd) + max(dat_dom_sp$mlgfd))/2,
                             (min(dat_dom_sp$mlgfd) + max(dat_dom_sp$mlgfd))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MRFD'] == 3.36",
                     "'95%CI' == '[1.68, 5.04]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mpd_all, y=domin),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='MPD', y='Dominance probability')+
    annotate(geom="text",x=c((min(dat_dom_sp$mpd_all) + max(dat_dom_sp$mpd_all))/2,
                             (min(dat_dom_sp$mpd_all) + max(dat_dom_sp$mpd_all))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MPD'] == '-0.02'",
                     "'95%CI' == '[-0.03, 0.00]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mfunc_d_all, y=domin),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='MFD', y='  ')+
    annotate(geom="text",x=c((min(dat_dom_sp$mfunc_d_all) + max(dat_dom_sp$mfunc_d_all))/2,
                             (min(dat_dom_sp$mfunc_d_all) + max(dat_dom_sp$mfunc_d_all))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MFD'] == 7.40",
                     "'95%CI' == '[-0.45, 15.07]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mnnd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance~mnnd
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mnnd.rdata")
(domin.mnnd.partial.logistic=ggplot(data=lincombs.data.domin.mnnd,
                                    aes(x=mnnd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[3],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[3],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnnd, y=domin),
               color = stata_pal("s2color")(4)[3],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    
    labs(x='MNND',
         y="Dominance probability")+
    annotate(geom="text",x=c((min(dat_dom_sp$mnnd) + max(dat_dom_sp$mnnd))/2,
                             (min(dat_dom_sp$mnnd) + max(dat_dom_sp$mnnd))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MNND'] == -0.91",
                     "'95%CI' == '[-1.55, -0.27]'"),parse=T,size=3.5) + 
    theme_regular()
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
     file = 'code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mnlgfd.rdata')

# Remarks: plogis() turns the posterior to the original scale, 
# inla.emarginal(fun,x) calculates the expected value after the turn,
# i.e., the mean, packed as a function of. inla.tmarginal(fun,x) turns the posterior,
# and then inla.qmarginal() calculates confidence intervals for the turned posterior.
# The above inla.emarginal() is accomplishing these two steps at once.

# Draw logistic curve for dominance ~ mnlgfd
load("code/results_analyzing/analysing_sameages_top50_data/lincombs.data.domin.mnlgfd.rdata")
(domin.mnlgfd.partial.logistic=ggplot(data=lincombs.data.domin.mnlgfd,aes(x=mnlgfd, y=predicted.value))+
    geom_ribbon(aes(ymin=lower,ymax=upper),
                fill=stata_pal("s2color")(4)[2],alpha=0.2)+
    geom_line(color=stata_pal("s2color")(4)[2],size=1)+
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnlgfd, y=domin),
               color = stata_pal("s2color")(4)[2],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='MNRFD', y='  ')+
    annotate(geom="text",x=c((min(dat_dom_sp$mnlgfd) + max(dat_dom_sp$mnlgfd))/2,
                             (min(dat_dom_sp$mnlgfd) + max(dat_dom_sp$mnlgfd))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MNRFD'] == 3.41",
                     "'95%CI' == '[2.32, 4.49]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mntd_all, y=domin),
               color = stata_pal("s2color")(4)[4],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='MNTD', y='Dominance probability')+
    annotate(geom="text",x=c((min(dat_dom_sp$mntd_all) + max(dat_dom_sp$mntd_all))/2,
                             (min(dat_dom_sp$mntd_all) + max(dat_dom_sp$mntd_all))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MNTD'] == '-0.00'",
                     "'95%CI' == '[-0.00, 0.00]'"),
             parse=T,size=3.5) + 
    theme_regular()
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
    scale_color_gradient(low = turbo(4)[2],
                         high = turbo(4)[1])+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    #theme(plot.margin=unit(c(0.4,0,0.4,0.4),units="lines"))+
    geom_point(data=dat_dom_sp, aes(x=mnfunc_d_all, y=domin),
               color = stata_pal("s2color")(4)[1],
               shape=1,
               alpha = 0.2,
               position=position_jitter(height=0.02))+
    labs(x='MNFD', y='  ')+
    annotate(geom="text",x=c((min(dat_dom_sp$mnfunc_d_all) + max(dat_dom_sp$mnfunc_d_all))/2,
                             (min(dat_dom_sp$mnfunc_d_all) + max(dat_dom_sp$mnfunc_d_all))/2),
             y=c(0.80,0.70),
             label=c("italic(β)['MNFD'] == 4.17",
                     "'95%CI' == '[0.22, 8.11]'"),
             parse=T,size=3.5) + 
    theme_regular()
)


#### Merge plots for Fig. S3 & S4 ####
library(cowplot)

Fig.S3 = ggdraw() +
  draw_plot(estab.mnd.partial.logistic, 0.02, 0.5, .48, .46) +
  draw_plot(estab.mlgfd.partial.logistic, 0.5, 0.5, .48, .46) +
  draw_plot(domin.mnd.partial.logistic, 0.02, 0.08, .48, .46) +
  draw_plot(domin.mlgfd.partial.logistic, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12
  )
emf('results/figures_sameages_top50/Fig.S3.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S3
dev.off() #turn off device and finalize file


Fig.S4 = ggdraw() +
  draw_plot(estab.mpd_all.partial.logistic, 0.02, 0.5, .48, .46) +
  draw_plot(estab.mfunc_d_all.partial.logistic, 0.5, 0.5, .48, .46) +
  draw_plot(domin.mpd_all.partial.logistic, 0.02, 0.08, .48, .46) +
  draw_plot(domin.mfunc_d_all.partial.logistic, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12
  )
emf('results/figures_sameages_top50/Fig.S4.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S4
dev.off() #turn off device and finalize file

#### Merge plots for Fig. S5 & S6 ####
library(cowplot)

Fig.S5 = ggdraw() +
  draw_plot(estab.mnnd.partial.logistic, 0.02, 0.5, .48, .46) +
  draw_plot(estab.mnlgfd.partial.logistic, 0.5, 0.5, .48, .46) +
  draw_plot(domin.mnnd.partial.logistic, 0.02, 0.08, .48, .46) +
  draw_plot(domin.mnlgfd.partial.logistic, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12
  )
emf('results/figures_sameages_top50/Fig.S5.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S5
dev.off() #turn off device and finalize file


Fig.S6 = ggdraw() +
  draw_plot(estab.mntd_all.partial.logistic, 0.02, 0.5, .48, .46) +
  draw_plot(estab.mnfunc_d_all.partial.logistic, 0.5, 0.5, .48, .46) +
  draw_plot(domin.mntd_all.partial.logistic, 0.02, 0.08, .48, .46) +
  draw_plot(domin.mnfunc_d_all.partial.logistic, 0.5, 0.08, .48, .46) +
  draw_plot_label(c("(a)", "(b)", "(c)", "(d)"),
                  c(0.015, 0.495, 0.015, 0.495),
                  c(0.97, 0.97, 0.55, 0.55)
                  ,size = 12)

emf('results/figures_sameages_top50/Fig.S6.emf',
    width = 20, height = 20, coordDPI = 1200, 
    emfPlusFontToPath=TRUE, # setting emfPlusFontToPath=TRUE to 
    # ensure text looks correct on the viewing system
    units = 'cm')
Fig.S6
dev.off() #turn off device and finalize file
