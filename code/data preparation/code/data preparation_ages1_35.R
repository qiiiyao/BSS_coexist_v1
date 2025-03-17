###spatio-temporal data
###data sorting
rm(list = ls())
library(dplyr)
library(data.table)

##### Data loading
sp_cover = read.csv('data/original data/BSS_community_332.csv',
                    header = T)
field = read.csv('data/original data/FIELDS.csv',
                 header = T)
trait = read.csv('data/original data/traits332.csv',
                 header = T)

sp_cover$f_p = paste(sp_cover$Field, sp_cover$Plot, sep = '_')
sp_cover = sp_cover %>% relocate(f_p, .after = Plot)
trait$species_id = paste('sp', trait$species_id, sep = '_')
sp_real = vector()
sp = colnames(sp_cover)[8:ncol(sp_cover)]

for (i in 1:length(sp)) {
  sp_real[i] = trait[trait$species_id == sp[i],]$Species
}

colnames(sp_cover)[8:ncol(sp_cover)] = sp_real # change the species name of data

### Find the time point that trees weren't appeared
source('code/data preparation/code/function/get_consec_rep.R')
sp_cover_tree = sp_cover %>% select(all_of(c(colnames(sp_cover)[1:7],
                                             intersect(sp_real,
                                                       trait[(trait$Growth=='Tree'),]$Species))))
sp_cover_tree_l = split(sp_cover_tree, sp_cover_tree$f_p)
sp_cover_tree_rowsums = sapply(sp_cover_tree_l, function(x){
  y = cbind(x[,1:7],
            tree = rowSums(x[,8:ncol(sp_cover_tree)]))
  y_1 = y %>% filter(tree < 30)
  standard_seq = unique(y_1$Age)
  ddd = standard_seq %in% y$Age
  sep_t_l = max(get_longest_consec_rep(ddd))
  sep_t_l = y$Age[sep_t_l]
  return(sep_t_l)})
mean(sp_cover_tree_rowsums)
### mean age 27, at this age tree cannot occupy more than 30% cover; therefore,
### we can use fitted results based on fake age 1-25 data to verify that our results
### can be achieved using no-tree data.

### Exclude trees for fitting
sp_cover = sp_cover %>% select(all_of(c(colnames(sp_cover)[1:7],
                                        intersect(sp_real, trait[!(trait$Growth=='Tree'),]$Species))))
sort(unique(sp_cover$Absolute_year))
sp_cover_f1 = sp_cover %>% filter(Age %in% c(1:51))
sp_cover_fp = split(sp_cover, sp_cover$f_p)
sp_cover_year = split(sp_cover, sp_cover$Absolute_year)
sp_cover_age = split(sp_cover, sp_cover$Age)
ages = sapply(sp_cover_age, function(x){y = unique(x$f_p)})

##### Follow Li et al., (2015), condense the two year
##### (half of field survey) to
##### one year and Follow Yin et al., (2022),
##### our study restricted to 3-51 age
sp_cover_f1 = sp_cover %>% filter(Absolute_year %in% c(1958:2017))
sp_cover_f1$fake_year = sp_cover_f1$Absolute_year
sp_cover_f1$fake_age = sp_cover_f1$Age
full_f_p = unique(sp_cover_f1$f_p)
for (i in seq(1959, 2017, 2)){
  #i = 1991
  f_p_year_2 = unique(sp_cover_f1[sp_cover_f1$Absolute_year == i-1,]$f_p)
  diff_f_p_year_2 = setdiff(full_f_p, f_p_year_2)
  row_diff = nrow(sp_cover_f1[sp_cover_f1$Absolute_year == i & 
                                sp_cover_f1$f_p %in% diff_f_p_year_2,
  ])
  if (row_diff == 0) {print(paste('cannot make it up', 'for',
                                  i-1, sep = ' '))}
  else {
    sp_cover_f1[sp_cover_f1$Absolute_year == i & sp_cover_f1$f_p %in% 
                  diff_f_p_year_2,]$fake_year = i-1
    age = unique(sp_cover_f1[sp_cover_f1$Absolute_year == i & sp_cover_f1$f_p %in% 
                               diff_f_p_year_2,]$Age)
    sp_cover_f1[sp_cover_f1$Absolute_year == i & sp_cover_f1$f_p %in%
                  diff_f_p_year_2,]$fake_age = rep(age-1,
                                                   each = row_diff/length(age-1))
  }
}

sp_cover_f1 = sp_cover_f1 %>% relocate(fake_year, .before = Age)
sp_cover_f1 = sp_cover_f1 %>% relocate(fake_age, .before = Field)
sp_cover_f1_y = split(sp_cover_f1, sp_cover_f1$fake_age)

sp_cover_f1_ages1_35 = sp_cover_f1 %>% filter(fake_age %in% c(1:35))

sp_cover_f1_ages1_35_fp = split(sp_cover_f1_ages1_35,
                                sp_cover_f1_ages1_35$f_p)
re_cover_ab_f1_ages1_35_fp = lapply(sp_cover_f1_ages1_35_fp, 
                                    function(w){
                                      re_cover_ab = cbind(w[,1:9],
                                                          w[,10:(ncol(w))][colSums(w[,
                                                                                     10:(ncol(w))]) != 0])
                                    })
save(re_cover_ab_f1_ages1_35_fp,
     file = 'code/data preparation/transformed_data/re_cover_ab_f1_ages1_35_fp.rdata')


sp_cover_exclude_tree_f1_allages = sp_cover_exclude_tree_f1 %>% filter(fake_age %in% seq(1, 51, 1))


###### Sum the species for field level and all level
sp_cover_f2_field = sp_cover_f2 %>% group_by(Absolute_year,             
                                             Year,                   
                                             Age,  Field) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)

sp_cover_f2_field = split(sp_cover_f2_field, sp_cover_f2_field$Field)

sp_cover_f2_all = sp_cover_f2 %>% group_by(Age) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)


######## Filter top species for fit under plot/field and all field level
### top 50: top25 exotics and top25 natives
freq_sps = apply(sp_cover_f1_ages1_35[,c(10:ncol(sp_cover_f1_ages1_35))], 2, function(x){length(x[x>0])})
freq_sps_exotic = freq_sps[(trait %>% filter(Origin == 'Exotic'))$Species]
freq_sps_exotic = freq_sps_exotic[which(!is.na(freq_sps_exotic))]
freq_sps_native = freq_sps[(trait %>% filter(Origin == 'Native'))$Species] 
freq_sps_native = freq_sps_native[which(!is.na(freq_sps_native))]
freq_sps_exotic_top_25 = sort(freq_sps_exotic, decreasing = T)[1:25]
freq_sps_native_top_25 = sort(freq_sps_native, decreasing = T)[1:25]
freq_sps_top_50 = c(freq_sps_exotic_top_25, freq_sps_native_top_25)

sp_cover_f1_ages1_35_top50 = cbind(sp_cover_f1_ages1_35[,c(1:9)],
                          sp_cover_f1_ages1_35 %>% select(names(freq_sps_top_50)))
sp_cover_f1_ages1_35_top50_fp = split(sp_cover_f1_ages1_35_top50, sp_cover_f1_ages1_35_top50$f_p)

### top25 exotics and all natives
freq_sps_native_all = freq_sps[(trait %>% filter(Origin == 'Native'))$Species] 
freq_sps_native_all = freq_sps_native_all[which(!is.na(freq_sps_native_all))]
freq_sps_top_25exotic_allnative = c(freq_sps_exotic_top_25, freq_sps_native_all)

sp_cover_ages1_35_top25_exotic = cbind(sp_cover_f1_ages1_35[,c(1:9)],
                                       sp_cover_f1_ages1_35 %>% select(names(freq_sps_top_25exotic_allnative)))
sp_cover_ages1_35_top25_exotic_fp = split(sp_cover_ages1_35_top25_exotic,
                                          sp_cover_ages1_35_top25_exotic$f_p)
save(sp_cover_ages1_35_top25_exotic_fp,
     file = 'code/data preparation/transformed data/sp_cover_ages1_35_top25_exotic_fp.rdata')


### top 40: top20 exotics and top20 natives
freq_sps = apply(sp_cover_f1_ages1_35[,c(10:ncol(sp_cover_f1_ages1_35))], 2, function(x){length(x[x>0])})
freq_sps_exotic = freq_sps[(trait %>% filter(Origin == 'Exotic'))$Species]
freq_sps_exotic = freq_sps_exotic[which(!is.na(freq_sps_exotic))]
freq_sps_native = freq_sps[(trait %>% filter(Origin == 'Native'))$Species] 
freq_sps_native = freq_sps_native[which(!is.na(freq_sps_native))]
freq_sps_exotic_top_20 = sort(freq_sps_exotic, decreasing = T)[1:20]
freq_sps_native_top_20 = sort(freq_sps_native, decreasing = T)[1:20]
freq_sps_top_40 = c(freq_sps_exotic_top_20, freq_sps_native_top_20)

sp_cover_f1_ages1_35_top40 = cbind(sp_cover_f1_ages1_35[,c(1:9)],
                                   sp_cover_f1_ages1_35 %>% select(names(freq_sps_top_40)))
sp_cover_f1_ages1_35_top40_fp = split(sp_cover_f1_ages1_35_top40, sp_cover_f1_ages1_35_top40$f_p)

### top20 exotics and all natives
freq_sps_native_all = freq_sps[(trait %>% filter(Origin == 'Native'))$Species] 
freq_sps_native_all = freq_sps_native_all[which(!is.na(freq_sps_native_all))]
freq_sps_top_20exotic_allnative = c(freq_sps_exotic_top_20, freq_sps_native_all)

sp_cover_ages1_35_top20_exotic = cbind(sp_cover_f1_ages1_35[,c(1:9)],
                                       sp_cover_f1_ages1_35 %>% select(names(freq_sps_top_20exotic_allnative)))
sp_cover_ages1_35_top20_exotic_fp = split(sp_cover_ages1_35_top20_exotic,
                                          sp_cover_ages1_35_top20_exotic$f_p)
save(sp_cover_ages1_35_top20_exotic_fp,
     file = 'code/data preparation/transformed data/sp_cover_ages1_35_top20_exotic_fp.rdata')


### top 30: top15 exotics and top15 natives
freq_sps = apply(sp_cover_f1_ages1_35[,c(10:ncol(sp_cover_f1_ages1_35))], 2, function(x){length(x[x>0])})
freq_sps_exotic = freq_sps[(trait %>% filter(Origin == 'Exotic'))$Species]
freq_sps_exotic = freq_sps_exotic[which(!is.na(freq_sps_exotic))]
freq_sps_native = freq_sps[(trait %>% filter(Origin == 'Native'))$Species] 
freq_sps_native = freq_sps_native[which(!is.na(freq_sps_native))]
freq_sps_exotic_top_15 = sort(freq_sps_exotic, decreasing = T)[1:15]
freq_sps_native_top_15 = sort(freq_sps_native, decreasing = T)[1:15]
freq_sps_top_30 = c(freq_sps_exotic_top_15, freq_sps_native_top_15)

sp_cover_f1_ages1_35_top30 = cbind(sp_cover_f1_ages1_35[,c(1:9)],
                                   sp_cover_f1_ages1_35 %>% select(names(freq_sps_top_30)))
sp_cover_f1_ages1_35_top30_fp = split(sp_cover_f1_ages1_35_top30, sp_cover_f1_ages1_35_top30$f_p)

### top15 exotics and all natives
freq_sps_native_all = freq_sps[(trait %>% filter(Origin == 'Native'))$Species] 
freq_sps_native_all = freq_sps_native_all[which(!is.na(freq_sps_native_all))]
freq_sps_top_15exotic_allnative = c(freq_sps_exotic_top_15, freq_sps_native_all)

sp_cover_ages1_35_top15_exotic = cbind(sp_cover_f1_ages1_35[,c(1:9)],
                                       sp_cover_f1_ages1_35 %>% select(names(freq_sps_top_15exotic_allnative)))
sp_cover_ages1_35_top15_exotic_fp = split(sp_cover_ages1_35_top15_exotic,
                                          sp_cover_ages1_35_top15_exotic$f_p)
save(sp_cover_ages1_35_top15_exotic_fp,
     file = 'code/data preparation/transformed data/sp_cover_ages1_35_top15_exotic_fp.rdata')


###### Preparing fitting ######
source('code/data preparation/code/function/get_consec_rep.R') 
source('code/data preparation/code/function/data collation_function.R')
source('code/data preparation/code/function/fit_prepa_function_filtered.R')
source('code/data preparation/code/function/fit_prepa_function_top50.R')
#source('code/data preparation/code/function/fit_prepa_function_original_field.R')

## top 50 sps for early succession (ages 1-35)
abundance_raw = sapply(re_cover_ab_f1_ages1_35_fp,
                       function(x){sum(x[,c(10:ncol(x))])})
richness_raw = sapply(re_cover_ab_f1_ages1_35_fp,
                      function(x){ncol(x[,c(10:ncol(x))])})

fit_fp_top50_ages1_35 = lapply(sp_cover_f1_ages1_35_top50_fp, fit_prepa_top)
fit_fp_top50_ages1_35_prop = lapply(sp_cover_f1_ages1_35_top50_fp, fit_prepa_top_prop)
intro_sp_top50_early_suc = sapply(fit_fp_top50_ages1_35,
                                  function(x){length(strsplit(x[["stage_sp_name"]][2], ', ')[[1]])})
length(intro_sp_top50_early_suc[which(intro_sp_top50_early_suc == 0)])/length(intro_sp_top50_early_suc)
domin_sp_top50_early_suc = sapply(fit_fp_top50_ages1_35,
                                  function(x){length(strsplit(x[["stage_sp_name"]][4], ', ')[[1]])})
length(intro_sp_top50_early_suc[which(intro_sp_top50_early_suc == 0)])/length(intro_sp_top50_early_suc)
mean(sapply(fit_fp_top50_ages1_35,
            function(x){ncol(x[["re_cover_ab"]])-9})) ## mean 37 species
abundance_pro = mean(sapply(fit_fp_top50_ages1_35,
                            function(x){sum(x[["re_cover_ab"]][,c(10:ncol(x[["re_cover_ab"]]))])})/abundance_raw)

richness_pro = mean(sapply(fit_fp_top50_ages1_35,
                           function(x){ncol(x[["re_cover_ab"]][,c(10:ncol(x[["re_cover_ab"]]))])})/richness_raw)
mean(sapply(fit_fp_top50_ages1_35,
            function(x){nrow(x[["re_cover_ab"]])})) 

save(fit_fp_top50_ages1_35,
     file = "code/data preparation/transformed data/fit_fp_top50_ages1_35.RData")
save(fit_fp_top50_ages1_35_prop,
     file = "code/data preparation/transformed data/fit_fp_top50_ages1_35_prop.RData")

## top 40 sps for early succession (ages 1-35)
fit_fp_top40_ages1_35 = lapply(sp_cover_f1_ages1_35_top40_fp, fit_prepa_top)
intro_sp_top40_early_suc = sapply(fit_fp_top40_ages1_35,
                                  function(x){length(strsplit(x[["stage_sp_name"]][2], ', ')[[1]])})
length(intro_sp_top40_early_suc[which(intro_sp_top40_early_suc == 0)])/length(intro_sp_top40_early_suc)
domin_sp_top40_early_suc = sapply(fit_fp_top40_ages1_35,
                                  function(x){length(strsplit(x[["stage_sp_name"]][4], ', ')[[1]])})
length(intro_sp_top40_early_suc[which(intro_sp_top40_early_suc == 0)])/length(intro_sp_top40_early_suc)
mean(sapply(fit_fp_top40_ages1_35,
            function(x){ncol(x[["re_cover_ab"]])-9})) ## mean 31.33958 species
save(fit_fp_top40_ages1_35,
     file = "code/data preparation/transformed data/fit_fp_top40_ages1_35.RData")


## top 30 sps for early succession (ages 1-35)
fit_fp_top30_ages1_35 = lapply(sp_cover_f1_ages1_35_top30_fp, fit_prepa_top)
intro_sp_top30_early_suc = sapply(fit_fp_top30_ages1_35,
                                  function(x){length(strsplit(x[["stage_sp_name"]][2], ', ')[[1]])})
length(intro_sp_top30_early_suc[which(intro_sp_top30_early_suc == 0)])/length(intro_sp_top30_early_suc)
domin_sp_top30_early_suc = sapply(fit_fp_top30_ages1_35,
                                  function(x){length(strsplit(x[["stage_sp_name"]][4], ', ')[[1]])})
length(intro_sp_top30_early_suc[which(intro_sp_top30_early_suc == 0)])/length(intro_sp_top30_early_suc)
mean(sapply(fit_fp_top30_ages1_35,
            function(x){ncol(x[["re_cover_ab"]])-9})) ## mean 24.675 species
save(fit_fp_top30_ages1_35,
     file = "code/data preparation/transformed data/fit_fp_top30_ages1_35.RData")


## top 50 sps merged into native/intro/estab and domin in early succession (ages 1-35)
source('code/data preparation/code/function/fit_prepa_function_neighbor_group.R')
fit_fp_same_ages_mergetop50_early_suc = lapply(sp_cover_f2_top50_early_suc_fp,
                                               fit_prepa_neighbor_group)
save(fit_fp_same_ages_mergetop50_early_suc,
     file = "code/data preparation/transformed data/fit_fp_same_ages_mergetop50_early_suc.RData")

#### Using filtered time-series data to judge different invasion stage
fake_age_all = sort(unique(sp_cover_f2$fake_age))

fit_stage_fp_t_all_filter_1 = list()
for (i in (1:14)) {
  #i = 14
  seq_t = c(i:(i+12))
  if (i == 1) {
    sp_cover_f2_f3_fp_f = lapply(sp_cover_f2_f3_fp, function(x){
      if (nrow(x) < 26) {
        x = NULL} else {
          x = x
        } 
    })
    sp_cover_f2_f3_fp_f <- sp_cover_f2_f3_fp_f[!sapply(sp_cover_f2_f3_fp_f,is.null)]
    fit_stage_fp_t_all_filter_1[[i]] = lapply(sp_cover_f2_f3_fp_f,
                                              fit_prepa,
                                              t = seq_t,
                                              seq_t_all = fake_age_all)
  } else {
    fit_stage_fp_t_all_filter_1[[i]] = lapply(sp_cover_f2_f3_fp,
                                              fit_prepa,
                                              t = seq_t,
                                              seq_t_all = fake_age_all)
  } 
  
}
#f_p_other = setdiff(names(fit_stage_fp_t_all_filter_1[[2]]),
#                    names(fit_stage_fp_t_all_filter_1[[1]]))
#fit_stage_fp_t_other = lapply(fit_stage_fp_t_all_filter_1,
#                               function(x){x[f_p_other]})

#save(fit_stage_fp_t_other,
#    file = './code/data preparation//transformed data/fit_stage_fp_t_other.rdata')
save(fit_stage_fp_t_all_filter_1,
     file = './code/data preparation//transformed data/fit_stage_fp_t_all_filter_1.rdata')

#### Using real time-series data to judge different invasion stage
fit_stage_fp_t_all_filter_2 = list()
for (i in (1:14)) {
  seq_t = c(i:(i+12))
  if (i == 1) {
    sp_cover_f2_f3_fp_f = lapply(sp_cover_f2_f3_fp, function(x){
      if (nrow(x) < 26) {
        x = NULL} else {
          x = x
        } 
    })
    sp_cover_f2_f3_fp_f <- sp_cover_f2_f3_fp_f[!sapply(sp_cover_f2_f3_fp_f,is.null)]
    fit_stage_fp_t_all_filter_2[[i]] = lapply(sp_cover_f2_f3_fp_f,
                                              fit_prepa_2,
                                              sp_cover_f2_origin = sp_cover_f2_fp,
                                              t = seq_t,
                                              seq_t_all = fake_age_all)
  } else {
    fit_stage_fp_t_all_filter_2[[i]] = lapply(sp_cover_f2_f3_fp,
                                              fit_prepa_2,
                                              sp_cover_f2_origin = sp_cover_f2_fp,
                                              t = seq_t,
                                              seq_t_all = fake_age_all)
  } 
}
save(fit_stage_fp_t_all_filter_2,
     file = './code/data preparation//transformed data/fit_stage_fp_t_all_filter_2.rdata')


#### Field level 
### filter some plots cannot cover all time span
fit_stage_field_all_time_strict = lapply(sp_cover_f2_f2_field, fit_prepa_field_strict,
                                         seq_t_all = fake_age_all)
save(fit_stage_field_all_time_strict,
     file = "code/data preparation//transformed data/fit_stage_field_all_time_strict.RData")

### If established species appeared in the top ***% most abundant species,
### defined it as dominant species
fit_stage_field_all_time_modest = lapply(sp_cover_f2_f2_field, fit_prepa_field_modest,
                                         seq_t_all = fake_age_all,
                                         domin_prop = 0.02)
save(fit_stage_field_all_time_modest,
     file = "code/data preparation/transformed data/fit_stage_field_all_time_modest.RData")

######## Time windows

fit_stage_field_t_all_filter_1 = list()
for (i in (1:14)) {
  seq_t = c(i:(i+12))
  if (i == 1) {
    sp_cover_f2_f3_field_f = lapply(sp_cover_f2_f3_field, function(x){
      if (nrow(x) < 26) {
        x = NULL} else {
          x = x
        } 
    })
    sp_cover_f2_f3_field_f <- sp_cover_f2_f3_field_f[!sapply(sp_cover_f2_f3_field_f,is.null)]
    fit_stage_field_t_all_filter_1[[i]] = lapply(sp_cover_f2_f3_field_f,
                                                 fit_prepa_field_modest,
                                                 t = seq_t,
                                                 seq_t_all = fake_age_all,
                                                 domin_prop = 0.02)
  } else {
    fit_stage_field_t_all_filter_1[[i]] = lapply(sp_cover_f2_f3_field,
                                                 fit_prepa_field_modest,
                                                 t = seq_t,
                                                 seq_t_all = fake_age_all,
                                                 domin_prop = 0.02)
  } 
}

save(fit_stage_field_t_all_filter_1,
     file = './code/data preparation/transformed data/fit_stage_field_t_all_filter_1.rdata')


