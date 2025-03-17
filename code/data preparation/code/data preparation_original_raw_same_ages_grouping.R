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

### Exclude trees for fitting
sp_cover = sp_cover %>% select(all_of(c(colnames(sp_cover)[1:7],intersect(sp_real, trait[!(trait$Growth=='Tree'),]$Species))))
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
sp_cover_f2 = sp_cover_f1 %>% filter(fake_age %in% c(1:51))
sp_cover_f2_fp = split(sp_cover_f2, sp_cover_f2$f_p)


######## Filter abundant species for fit under plot/field and all field level
source('code/data preparation/code/function/get_consec_rep.R')
source('code/data preparation/code/function/data collation_function.R')
source('code/data preparation/code/function/fit_prepa_function_filtered.R')
source('code/data preparation/code/function/fit_prepa_function_top50.R')
#source('code/data preparation/code/function/fit_prepa_function_original_field.R')

freq_sps = apply(sp_cover_f2[,c(10:ncol(sp_cover_f2))], 2, function(x){length(x[x>0])})
freq_sps_exotic = freq_sps[(trait %>% filter(Origin == 'Exotic'))$Species]
freq_sps_exotic = freq_sps_exotic[which(!is.na(freq_sps_exotic))]
freq_sps_native = freq_sps[(trait %>% filter(Origin == 'Native'))$Species] 
freq_sps_native = freq_sps_native[which(!is.na(freq_sps_native))]
freq_sps_exotic_top_25 = sort(freq_sps_exotic, decreasing = T)[1:25]
freq_sps_native_top_25 = sort(freq_sps_native, decreasing = T)[1:25]
freq_sps_top_50 = c(freq_sps_exotic_top_25, freq_sps_native_top_25)

sp_cover_f2_top50 = cbind(sp_cover_f2[,c(1:9)],
                          sp_cover_f2 %>% select(names(freq_sps_top_50)))

sp_cover_f2_top50_fp = split(sp_cover_f2_top50, sp_cover_f2_top50$f_p)

#### Plot level 
fit_fp_same_ages_top50 = lapply(sp_cover_f2_top50_fp, fit_prepa_top)
intro_sp_top50 = sapply(fit_fp_same_ages_top50,
                      function(x){length(strsplit(x[["stage_sp_name"]][2], ', ')[[1]])})
length(intro_sp_top50[which(intro_sp_top50 == 0)])/length(intro_sp_top50)
domin_sp_top50 = sapply(fit_fp_same_ages_top50,
                        function(x){length(strsplit(x[["stage_sp_name"]][4], ', ')[[1]])})
length(intro_sp_top50[which(intro_sp_top50 == 0)])/length(intro_sp_top50)
mean(sapply(fit_fp_same_ages_top50,
            function(x){ncol(x[["re_cover_ab"]])-9})) ## mean 40 species
save(fit_fp_same_ages_top50,
     file = "code/data preparation/transformed data/fit_fp_same_ages_top50.RData")

fit_fp_same_ages_0.3_raw = lapply(sp_cover_f2_fp, fit_prepa,
                                  nati_per = 0.3,intro_per = 0.3,
                                  estab_per = 0.3, domin_per = 0.3)
intro_sp_0.3 = sapply(fit_fp_same_ages_0.3_raw,
                      function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
length(intro_sp_0.3[which(intro_sp_0.3 == 0)])/length(intro_sp_0.3)
#save(fit_fp_same_ages_0.3_raw,
 #    file = "code/data preparation/transformed data/fit_fp_same_ages_0.3_raw.RData")

fit_fp_same_ages_0.25_raw = lapply(sp_cover_f2_fp, fit_prepa, nati_per = 0.25,
                                   intro_per = 0.25, estab_per = 0.25,
                                   domin_per = 0.25)
intro_sp_0.25 = sapply(fit_fp_same_ages_0.25_raw,
                       function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
mean(sapply(fit_fp_same_ages_0.25_raw,
            function(x){ncol(x[["re_cover_ab_f"]])-9}))
mean(sapply(fit_fp_same_ages_0.25_raw,
            function(x){ncol(x[["re_cover_ab"]])-9}))
mean(sapply(fit_fp_same_ages_0.25_raw,
            function(x){nrow(x[["re_cover_ab_f"]])}))
35.3*15/(15*15+15)
length(intro_sp_0.25[which(intro_sp_0.25 == 0)])/length(intro_sp_0.25)
save(fit_fp_same_ages_0.25_raw,
     file = "code/data preparation/transformed data/fit_fp_same_ages_0.25_raw.RData")

fit_fp_same_ages_0.20_raw = lapply(sp_cover_f2_fp, fit_prepa, nati_per = 0.20,
                                   intro_per = 0.20, estab_per = 0.20,
                                   domin_per = 0.20)
intro_sp_0.20 = sapply(fit_fp_same_ages_0.20_raw,
                       function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
mean(sapply(fit_fp_same_ages_0.20_raw,
            function(x){ncol(x[["re_cover_ab_f"]])-9}))
sp = unique(unlist(sapply(fit_fp_same_ages_0.20_raw,
                          function(x){colnames(x[["re_cover_ab_f"]])[10:ncol(x[["re_cover_ab_f"]])]})))
abundance_pro = mean(sapply(fit_fp_same_ages_0.20_raw,
                            function(x){sum(x[["re_cover_ab_f"]][,c(10:ncol(x[["re_cover_ab_f"]]))])/sum(x[["re_cover_ab"]][,c(10:ncol(x[["re_cover_ab"]]))])}))

richness_pro = mean(sapply(fit_fp_same_ages_0.20_raw,
                           function(x){ncol(x[["re_cover_ab_f"]][,c(10:ncol(x[["re_cover_ab_f"]]))])/ncol(x[["re_cover_ab"]][,c(10:ncol(x[["re_cover_ab"]]))])}))
length(intro_sp_0.20[which(intro_sp_0.20 == 0)])/length(intro_sp_0.20)
save(fit_fp_same_ages_0.20_raw,
     file = "code/data preparation/transformed data/fit_fp_same_ages_0.20_raw.RData")

fit_fp_same_ages_0.15_raw = lapply(sp_cover_f2_fp, fit_prepa,
                                   nati_per = 0.15,
                                   intro_per = 0.15,
                                   estab_per = 0.15,
                                   domin_per = 0.15)
intro_sp_0.15 = sapply(fit_fp_same_ages_0.15_raw,
                       function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
length(intro_sp_0.15[which(intro_sp_0.15 == 0)])/length(intro_sp_0.15)
save(fit_fp_same_ages_0.15_raw,
     file = "code/data preparation/transformed data/fit_fp_same_ages_0.15_raw.RData")


fit_field_same_ages_0.40_raw = lapply(sp_cover_f2_field, 
                                      fit_prepa_original_field,
                                      nati_per = 0.40,
                                      intro_per = 0.40,
                                      estab_per = 0.40,
                                      domin_per = 0.40)
intro_sp_0.40 = sapply(fit_field_same_ages_0.40_raw,
                       function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
mean(sapply(fit_field_same_ages_0.40_raw,
            function(x){ncol(x[["re_cover_ab_f"]])-4}))
sp = unique(unlist(sapply(fit_field_same_ages_0.40_raw,
                          function(x){colnames(x[["re_cover_ab_f"]])[5:ncol(x[["re_cover_ab_f"]])]})))
length(intro_sp_0.40[which(intro_sp_0.40 == 0)])/length(intro_sp_0.40)
save(fit_field_same_ages_0.40_raw,
     file = "code/data preparation/transformed data/fit_field_same_ages_0.40_raw.RData")

fit_field_same_ages_0.80_raw = lapply(sp_cover_f2_field, fit_prepa_original_field,
                                         nati_per = 0.80,
                                         intro_per = 0.80, estab_per = 0.80, domin_per = 0.80)
unique(sp_cover_f2)
intro_sp_0.80 = sapply(fit_field_same_ages_0.80_raw,
                       function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
mean(sapply(fit_field_same_ages_0.80_raw,
            function(x){ncol(x[["re_cover_ab_f"]])-4}))
mean(sapply(fit_field_same_ages_0.80_raw,
            function(x){nrow(x[["re_cover_ab_f"]])}))
sp = unique(unlist(sapply(fit_field_same_ages_0.80_raw,
                          function(x){colnames(x[["re_cover_ab_f"]])[5:ncol(x[["re_cover_ab_f"]])]})))
length(intro_sp_0.80[which(intro_sp_0.80 == 0)])/length(intro_sp_0.80)
save(fit_field_same_ages_0.80_raw,
     file = "code/data preparation/transformed data/fit_field_same_ages_0.80_raw.RData")

#### Using filtered time-series data to judge different invasion stage
fake_age_all = sort(unique(sp_cover_f2_f3$fake_age))

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


