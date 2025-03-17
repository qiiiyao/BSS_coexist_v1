#### Test INLA
load('code/data preparation/fit_stage_fp_alltime.rdata')
library(INLA)
pcprior = list(prec=list("pc.prec", param=c(0.1,0.01)))


x = inva_stage_fp_all_time[[1]]
data = x[[5]]
sp_name = colnames(data$X)
x = data[[1]]

  data_fit = data.frame(G = x$G, time = x$time,
                        x$X)
  data_fit_l = as.list(data_fit)
  
  formula = eval(parse(text=paste(colnames(data_fit)[1],
                                  paste(paste(paste(colnames(data_fit)[3:ncol(data_fit)],
                                                    collapse = '+'),
                                              sep = ''),
                                        '+f(time, model= "ar1")',
                                        sep = ''),
                                  sep = '~')))
  
  formula_1 = eval(parse(text=paste(colnames(data_fit)[1],
                                    paste(paste(paste(colnames(data_fit)[3:ncol(data_fit)],
                                                      collapse = '+'),
                                                sep = ''),
                                          sep = ''),
                                    sep = '~')))
  
  fit = inla(formula,
             #family = 'gaussian',
             data = data_fit_l,
             control.compute = list(dic=T,cpo=T,
                                    waic = T))
  
  fit_1 = inla(formula_1,
               #family = 'gaussian',
               data = data_fit_l,
               control.compute = list(dic=T,cpo=T,
                                      waic = T)
  )
  dd = summary(fit)
  dd_1 = summary(fit_1)
  dd$fixed
  
  fit_all = list()
for (i in 1:464) {
  x = inva_stage_fp_all_time[[i]]
  data = x[[5]]
  x$re_cover_ab_f
  fit_l = lapply(data, function(y){
    y = data[[10]]
    data_fit = data.frame(G = y$G, time = y$time,
                          y$X)
    data_fit_l = as.list(data_fit)
    
    formula = eval(parse(text=paste(colnames(data_fit)[1],
                                    paste(paste(paste(colnames(data_fit)[3:ncol(data_fit)],
                                                      collapse = '+'),
                                                sep = ''),
                                          '+f(time, model= "ar1")',
                                          sep = ''),
                                    sep = '~')))
    
    fit = inla(formula,
               #family = 'gaussian',
               data = data_fit_l,
               control.compute = list(dic=T,cpo=T,
                                      waic = T),
               verbose=TRUE)
    
    dd = summary(fit)
    return(dd)
  })
  fit_all[[i]] = fit_l
}



lCtr <- lmeControl(maxIter = 500, msMaxIter = 500,
                   tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

formula = eval(parse(text=paste(colnames(data_fit)[1],
                                paste(paste(paste(colnames(data_fit)[3:ncol(data_fit)],
                                                  collapse = '+'),
                                            sep = '')),
                                sep = '~')))
library(nlme)
fit_1 <- gls(formula,
                  data= data_fit, 
                  control=lCtr,
                  correlation=corAR1(form=~time), method='REML',
                  na.action=na.omit)
summary(fit_1)
fit_lme = lm(eval(parse(text=paste(colnames(data_fit)[1],
                                   paste(paste(paste(colnames(data_fit)[3:ncol(data_fit)],
                                                     collapse = '+'),
                                               sep = '')),
                                   sep = '~'))),
             data = data_fit)
summary(fit_lme)