library('cmdstanr')
library("coda")
library("deSolve")
library("BayesianTools")
library("parallel")
library("posterior")
library("dplyr")
library('brms')

source('F:\\PG\\loo\\twotwo_ran\\long\\data preparation\\data preparation for abcd_sp_rare_fg.R') # get all abcd sp obsevration data
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
setwd("F:/PG/loo/twotwo_ran/recent_2.0_t2.5")
### load data

### load work dictionary
file1 = file.path("Ricker_random_time_loo.stan")
mod1 = cmdstan_model(file1)
file2 = file.path("BH_random_time_loo.stan")
mod2 = cmdstan_model(file2)
file3 = file.path("BH_b_random_time_loo.stan")
mod3 = cmdstan_model(file3)
file4 = file.path("exp_r.stan")
mod4 = cmdstan_model(file4)
file5 = file.path("BH_log(n+1)_random_time_loo.stan")
mod5 = cmdstan_model(file5)
file6 = file.path("BH_1+b_random_time_loo.stan")
mod6 = cmdstan_model(file6)
file7 = file.path("Ricker_log(n+1)random_time_loo.stan")
mod7 = cmdstan_model(file7)
file8 = file.path("BH_n^a_random_time.stan")
mod8 = cmdstan_model(file8)



chains = 4 # 4 chains in parallel
sampling = 2500
warmup = 1000
thin = 1
badloo_mod1 <- vector()
badloo_mod2 <- vector()
badloo_mod3 <- vector()
badloo_mod4 <- vector()
badloo_mod5 <- vector()
badloo_mod6 <- vector()
badloo_mod7 <- vector()
badloo_mod8 <- vector()
all_abcd_sp[[76]] <- all_abcd_sp[[76]][-c(1,2),]
all_abcd_sp[[77]] <- all_abcd_sp[[77]][-c(1,2),]
all_abcd_sp[[78]] <- all_abcd_sp[[78]][-c(1,2),]
# 27 29 63 64
for (z in c(65:79)) {
  # z = 1
  coexist_trans_sp = all_abcd_sp[[z]]
  coexist_trans_sp[coexist_trans_sp == 0] = 0.01 # change 0 to 1e-2
  coexist_trans_sp = arrange(coexist_trans_sp,
                             coexist_trans_sp$year,
                             coexist_trans_sp$rep)
  
  year = unique(coexist_trans_sp$year)
  rep = unique(coexist_trans_sp$rep)
  
  time = length(unique(coexist_trans_sp$year))
  nrep = length(unique(coexist_trans_sp$rep))
  nplot = length(unique(coexist_trans_sp$plot))
  
    sp_bio = colnames(coexist_trans_sp[6:ncol(coexist_trans_sp)])
    library('dplyr')
    
    
    ### since 1991
    
    coexist_trans_sp = arrange(coexist_trans_sp, coexist_trans_sp$year)
    
    coexist_trans_sp_since = coexist_trans_sp %>%
      filter(rep %in% c(1:6))
    coexist_trans_sp_since_1 = coexist_trans_sp_since[(nrep+1):(nrep*(time-1)),] %>%
      slice(rep(1:n(), each = 2))
    
    coexist_trans_sp_since_2 = coexist_trans_sp_since[c(1:nrep,((nrep*(time-1))+1):(nrep*time)),]
    
    coexist_trans_sp_since_3 = rbind(coexist_trans_sp_since_1,
                                     coexist_trans_sp_since_2)
    
    coexist_trans_sp_since_3 = arrange(coexist_trans_sp_since_3,
                                       coexist_trans_sp_since_3$rep,
                                       coexist_trans_sp_since_3$year)
    
    coexist_trans_sp_since_4 = cbind(plot = coexist_trans_sp_since_3[,1],
                                     year2 = rep(c(1,2), nrep*(time-1)),
                                     coexist_trans_sp_since_3[,
                                                              2:ncol(coexist_trans_sp_since_3)])
    
    row.names(coexist_trans_sp_since_4) = c(1:nrow(coexist_trans_sp_since_4))
    initime_since = filter(coexist_trans_sp_since_4, year2 == '1')
    nextime_since = filter(coexist_trans_sp_since_4, year2 == '2')
    initime_since_1 = initime_since[,-c(1:6)]
    dndt_since = nextime_since[,7:ncol(nextime_since)]/initime_since[,7:ncol(nextime_since)]
    colnames(dndt_since) = paste('dndt', colnames(dndt_since), sep = '.')
    mod_dat_since = cbind(dndt_since,initime_since_1)
    
    mod_dat_trans_sp_since = cbind(time = rep(c(1:(time-1)),nrep),
                                   rep = c(rep(1,(time-1)), rep(2,(time-1)), rep(3,(time-1)),
                                           rep(4,(time-1)),
                                           rep(5,(time-1)), rep(6,(time-1))), mod_dat_since)
    
    mod_dat_trans_sp = mod_dat_trans_sp_since
    
    mod_dat_trans_sp$ftime = factor(mod_dat_trans_sp$time)
    
    Growthrates_trans_sp = colnames(mod_dat_trans_sp)[3:(3+length(sp_bio)-1)]
    Bio_trans_sp = colnames(mod_dat_trans_sp)[(3+length(sp_bio)):(3+2*length(sp_bio)-1)]
    
    ### hierarchical bayesian model
    
    sdscal_trans_sp = vector()
    Bio_trans_sp
    for (j in 1:length(sp_bio)) {
      gr = Growthrates_trans_sp[j]
      sdscal_trans_sp[j] = sd(residuals(lm(paste(gr,
                                                 paste(Bio_trans_sp, collapse = "+"),
                                                 sep = '~'),
                                           data = mod_dat_trans_sp)))
    }
    
   
    list_trans1 = list()
    list_trans2 = list()
    list_trans3 = list()
    list_trans4 = list()
    list_trans5 = list()
    list_trans6 = list()
    list_trans7 = list()
    list_trans8 = list()
    ### fit and extract
    
    for (k in 1:length(sp_bio)) {
      # k = 1
      data = list(N = length(mod_dat_trans_sp[,1]),
                  P = length(sp_bio),
                  X = as.matrix(mod_dat_trans_sp[,(2+length(sp_bio)+1):(2+2*length(sp_bio))]),
                  ftime = length(unique(mod_dat_trans_sp$ftime)),
                  time = as.numeric(mod_dat_trans_sp$ftime),
                  #sdscal = sdscal_trans_sp[k],
                  G = mod_dat_trans_sp[,2+k])
      
      list_trans1[[k]] = mod1$sample(data = data, seed = 123,
                                     chains = chains,
                                     #  init = 50,
                                     parallel_chains = chains,
                                     iter_warmup = warmup,
                                     iter_sampling = sampling,
                                     refresh = 10, # print update every 10 iters
                                     max_treedepth = 10,
                                     adapt_delta = 0.8)
      
      loofit1 <- list_trans1[[k]]$loo()
      if(max(loofit1[["diagnostics"]][["pareto_k"]]) < 0.7){
        estimates <- as.matrix(loofit1[["estimates"]])
        write.csv(loofit1[["estimates"]], paste(paste("./loo/loo_", unique(coexist_trans_sp$plot),'_',sp_bio[k],"_mod1.csv", sep="")))
      } else{badloo_mod1 <- c(badloo_mod1,paste(paste(unique(coexist_trans_sp$plot),'_',sp_bio[k], sep=""))) }
      write.csv(list_trans1[[k]] $summary(), paste("./rds/r_n_eff_", unique(coexist_trans_sp$plot),'_',sp_bio[k] ,"_mod1.csv", sep=""))
      list_trans2[[k]] = mod2$sample(data = data, seed = 123,
                                     chains = chains,
                                     #  init = 50,
                                     parallel_chains = chains,
                                     iter_warmup = warmup,
                                     iter_sampling = sampling,
                                     refresh = 10, # print update every 10 iters
                                     max_treedepth = 10,
                                     adapt_delta = 0.8)
      
      loofit2 <- list_trans2[[k]]$loo()
      if(max(loofit2[["diagnostics"]][["pareto_k"]]) < 0.7){
        estimates <- as.matrix(loofit2[["estimates"]])
        write.csv(loofit2[["estimates"]], paste(paste("./loo/loo_", unique(coexist_trans_sp$plot),'_',sp_bio[k],"_mod2.csv", sep="")))
      } else{badloo_mod2 <- c(badloo_mod2,paste(paste(unique(coexist_trans_sp$plot),'_',sp_bio[k], sep=""))) }
      write.csv(list_trans2[[k]] $summary(), paste("./rds/r_n_eff_", unique(coexist_trans_sp$plot),'_',sp_bio[k] ,"_mod2.csv", sep=""))
      list_trans3[[k]] = mod3$sample(data = data, seed = 123,
                                     chains = chains,
                                     #  init = 50,
                                     parallel_chains = chains,
                                     iter_warmup = warmup,
                                     iter_sampling = sampling,
                                     refresh = 10, # print update every 10 iters
                                     max_treedepth = 10,
                                     adapt_delta = 0.8)
      
      loofit3 <- list_trans3[[k]]$loo()
      if(max(loofit3[["diagnostics"]][["pareto_k"]]) < 0.7){
        estimates <- as.matrix(loofit3[["estimates"]])
        write.csv(loofit3[["estimates"]], paste(paste("./loo/loo_", unique(coexist_trans_sp$plot),'_',sp_bio[k],"_mod3.csv", sep="")))
      } else{badloo_mod3 <- c(badloo_mod3,paste(paste(unique(coexist_trans_sp$plot),'_',sp_bio[k], sep=""))) }
      write.csv(list_trans3[[k]] $summary(), paste("./rds/r_n_eff_", unique(coexist_trans_sp$plot),'_',sp_bio[k] ,"_mod3.csv", sep=""))
      
      list_trans4[[k]] = mod4$sample(data = data, seed = 123,
                                     chains = chains,
                                     #  init = 50,
                                     parallel_chains = chains,
                                     iter_warmup = warmup,
                                     iter_sampling = sampling,
                                     refresh = 10, # print update every 10 iters
                                     max_treedepth = 10,
                                     adapt_delta = 0.8)
      
      loofit4 <- list_trans4[[k]]$loo()
      if(max(loofit4[["diagnostics"]][["pareto_k"]]) < 0.7){
        estimates <- as.matrix(loofit4[["estimates"]])
        write.csv(loofit4[["estimates"]], paste(paste("./loo/loo_", unique(coexist_trans_sp$plot),'_',sp_bio[k],"_mod4.csv", sep="")))
      } else{badloo_mod4 <- c(badloo_mod4,paste(paste(unique(coexist_trans_sp$plot),'_',sp_bio[k], sep=""))) }
      write.csv(list_trans4[[k]] $summary(), paste("./rds/r_n_eff_", unique(coexist_trans_sp$plot),'_',sp_bio[k] ,"_mod4.csv", sep=""))
      list_trans5[[k]] = mod5$sample(data = data, seed = 123,
                                     chains = chains,
                                     #  init = 50,
                                     parallel_chains = chains,
                                     iter_warmup = warmup,
                                     iter_sampling = sampling,
                                     refresh = 10, # print update every 10 iters
                                     max_treedepth = 10,
                                     adapt_delta = 0.8)
      
      loofit5 <- list_trans5[[k]]$loo()
      if(max(loofit5[["diagnostics"]][["pareto_k"]]) < 0.7){
        estimates <- as.matrix(loofit5[["estimates"]])
        write.csv(loofit5[["estimates"]], paste(paste("./loo/loo_", unique(coexist_trans_sp$plot),'_',sp_bio[k],"_mod5.csv", sep="")))
      } else{badloo_mod5 <- c(badloo_mod5,paste(paste(unique(coexist_trans_sp$plot),'_',sp_bio[k], sep=""))) }
      write.csv(list_trans5[[k]] $summary(), paste("./rds/r_n_eff_", unique(coexist_trans_sp$plot),'_',sp_bio[k] ,"_mod5.csv", sep=""))
      
      list_trans6[[k]] = mod6$sample(data = data, seed = 123,
                         chains = chains,
                         #  init = 50,
                         parallel_chains = chains,
                         iter_warmup = warmup,
                         iter_sampling = sampling,
                         refresh = 10, # print update every 10 iters
                         max_treedepth = 10,
                         adapt_delta = 0.8)
      
      loofit6 <- list_trans6[[k]]$loo()
      if(max(loofit6[["diagnostics"]][["pareto_k"]]) < 0.7){
        estimates <- as.matrix(loofit6[["estimates"]])
        write.csv(loofit6[["estimates"]], paste(paste("./loo/loo_", unique(coexist_trans_sp$plot),'_',sp_bio[k],"_mod6.csv", sep="")))
      } else{badloo_mod6 <- c(badloo_mod6,paste(paste(unique(coexist_trans_sp$plot),'_',sp_bio[k], sep=""))) }
      write.csv(list_trans6[[k]] $summary(), paste("./rds/r_n_eff_", unique(coexist_trans_sp$plot),'_',sp_bio[k] ,"_mod6.csv", sep=""))
      
      list_trans7[[k]] = mod7$sample(data = data, seed = 123,
                                     chains = chains,
                                     #  init = 50,
                                     parallel_chains = chains,
                                     iter_warmup = warmup,
                                     iter_sampling = sampling,
                                     refresh = 10, # print update every 10 iters
                                     max_treedepth = 10,
                                     adapt_delta = 0.8)
      
      loofit7 <- list_trans7[[k]]$loo()
      if(max(loofit7[["diagnostics"]][["pareto_k"]]) < 0.7){
        estimates <- as.matrix(loofit7[["estimates"]])
        write.csv(loofit7[["estimates"]], paste(paste("./loo/loo_", unique(coexist_trans_sp$plot),'_',sp_bio[k],"_mod7.csv", sep="")))
      } else{badloo_mod7 <- c(badloo_mod7,paste(paste(unique(coexist_trans_sp$plot),'_',sp_bio[k], sep=""))) }
      write.csv(list_trans7[[k]] $summary(), paste("./rds/r_n_eff_", unique(coexist_trans_sp$plot),'_',sp_bio[k] ,"_mod7.csv", sep=""))
      
      list_trans8[[k]] = mod8$sample(data = data, seed = 123,
                                     chains = chains,
                                     #  init = 50,
                                     parallel_chains = chains,
                                     iter_warmup = warmup,
                                     iter_sampling = sampling,
                                     refresh = 10, # print update every 10 iters
                                     max_treedepth = 10,
                                     adapt_delta = 0.8)
      
      loofit8 <- list_trans8[[k]]$loo()
      if(max(loofit8[["diagnostics"]][["pareto_k"]]) < 0.7){
        estimates <- as.matrix(loofit8[["estimates"]])
        write.csv(loofit8[["estimates"]], paste(paste("./loo/loo_", unique(coexist_trans_sp$plot),'_',sp_bio[k],"_mod8.csv", sep="")))
      } else{badloo_mod8 <- c(badloo_mod8,paste(paste(unique(coexist_trans_sp$plot),'_',sp_bio[k], sep=""))) }
      write.csv(list_trans8[[k]] $summary(), paste("./rds/r_n_eff_", unique(coexist_trans_sp$plot),'_',sp_bio[k] ,"_mod8.csv", sep=""))
      
    }
}



### load packages
library('cmdstanr')
library("coda")
library("deSolve")
library("BayesianTools")
library("parallel")
library("posterior")
library("dplyr")
library(brms)

check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

### load data

write.csv(badloo_mod1,'badloo_mod1.csv')
write.csv(badloo_mod2,'badloo_mod2.csv')
write.csv(badloo_mod3,'badloo_mod3.csv')
write.csv(badloo_mod4,'badloo_mod4.csv')
write.csv(badloo_mod5,'badloo_mod5.csv')
write.csv(badloo_mod6,'badloo_mod6.csv')
write.csv(badloo_mod7,'badloo_mod7.csv')
write.csv(badloo_mod8,'badloo_mod8.csv')

###### merge loo ####
library(dplyr)
library(stringr)
library(data.table)
summary_r = list.files('./loo', full.names = T)[grep(("loo.*mod8.csv"),
                                                     list.files('./loo', full.names = T))]

trans_summary_all_mod8 <- data.frame()
for (i in 1:length(summary_r)) {
  # i = 12
  trans_summary = read.csv(summary_r[i], header = T)
  trans_summary$mod <- 'mod8'
  plot_sp <- sub(".*/loo_(.*?)_mod8\\.csv", "\\1", summary_r[i])
  trans_summary$plot_sp <- plot_sp
  trans_summary_all_mod8 = rbind(trans_summary_all_mod8,trans_summary)
}

loo_all <- rbind(trans_summary_all_mod1, trans_summary_all_mod2,trans_summary_all_mod3,
                 trans_summary_all_mod4,trans_summary_all_mod5,trans_summary_all_mod6,
  trans_summary_all_mod7,trans_summary_all_mod8)
write.csv(loo_all,'loo_all.csv')
#### plot_loo####
trans_summary_all <- read.csv('loo_all.csv')
# 计算某一列中每个值出现的次数
value_counts <- table(trans_summary_all$plot_sp)

# 找出出现次数为8的值
values_with_nine_occurrences <- names(value_counts[value_counts == 24])

# 根据出现次数为3*3的值来筛选行
subset_df <- trans_summary_all[trans_summary_all$plot_sp %in% values_with_nine_occurrences, ]
subset_df <- subset(subset_df,subset_df$X == 'elpd_loo')
# 按照'plot'列分割数据框
split_df <- split(subset_df, subset_df$plot_sp)
colnames(subset_df)
# 对每个小数据框进行操作
modified_dfs <- lapply(split_df, function(sub_df) {
  sub_df <- arrange(sub_df, desc(Estimate) )
  sub_df$adjust_esti <- sub_df$Estimate + abs(max(sub_df$Estimate)) 
  return(sub_df)
})

# 重新拼接成一个数据框
new_df <- do.call(rbind, modified_dfs)
unique(new_df$plot_sp)
colnames(new_df)
library(ggsignif)
library(ggplot2)

ran_figure_ab <- ggplot(new_df, aes(x = mod, y = adjust_esti))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 3) +
  scale_x_discrete(labels = c('ricker',"BH", "BH(1+na)^b",'null', "BH(log)",'BH1+(na)^b','ricker(log n+1)','ricker_n^a')) 
ggsave("ran_figure_ab.png", ran_figure_ab, width = 7, height = 8, units = "in")
ran_figure <- ggplot(new_df, aes(x = mod, y = Estimate))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 3) +
  scale_x_discrete(labels = c('ricker',"BH", "BH(1+na)^b",'null', "BH(log)",'BH1+(na)^b','ricker(log n+1)','ricker_n^a')) 
ggsave("ran_figure.png", ran_figure, width = 7, height = 8, units = "in")
###### check rhat ####
library(dplyr)
library(stringr)
library(data.table)
badhat <- vector()
summary_r = list.files('./rds',
                       full.names = T)[grep(("r_n_eff"),
                                            list.files('./rds',
                                                       full.names = T))]

for (i in 1:length(summary_r)) {
  #i = 12
  trans_summary = read.csv(summary_r[i], header = T)
  trans_summary_c = trans_summary %>% 
    filter(grepl("^a", variable) |
             grepl('r\\[', variable)|
             grepl("^b", variable)) 
  if(max(trans_summary_c$rhat) > 1.05){
    badhat <- c(badhat,summary_r[i])  }
}
write.csv(badhat,'badhat.csv')
