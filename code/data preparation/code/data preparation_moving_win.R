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

sp_cover_fp = split(sp_cover, sp_cover$f_p)
sp_cover_year = split(sp_cover, sp_cover$Absolute_year)
sp_cover_age = split(sp_cover, sp_cover$Age)

sapply(sp_cover_fp, function(x) {sort(unique(x$Absolute_year))})
sapply(sp_cover_year, function(x) {sort(unique(x$f_p))})

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

sp_cover_fakeage1 = sp_cover_f1 %>% filter(fake_age == 1)

sp_cover_f2 = sp_cover_f1 %>% filter(fake_age %in% seq(1, 51, 2) &
                                     f_p %in% sp_cover_fakeage1$f_p)
sp_cover_f2 = sp_cover_f2 %>% relocate(fake_year, .before = Age)
sp_cover_f2 = sp_cover_f2 %>% relocate(fake_age, .before = Field)
sp_cover_f2_y = split(sp_cover_f2, sp_cover_f2$fake_age)
sapply(sp_cover_f2_y, function(x) {sort(unique(x$f_p))})
unique(sp_cover_f2$fake_age)
sp_cover_f2_fp = split(sp_cover_f2, sp_cover_f2$f_p)

#### Not delete other plots for including age1
sp_cover_f3 = sp_cover_f1 %>% filter(fake_age %in% seq(1, 51, 2))
sp_cover_f3 = sp_cover_f3 %>% relocate(fake_year, .before = Age)
sp_cover_f3 = sp_cover_f3 %>% relocate(fake_age, .before = Field)
sp_cover_f3_y = split(sp_cover_f3, sp_cover_f3$fake_age)
sapply(sp_cover_f3_y, function(x) {sort(unique(x$f_p))})
unique(sp_cover_f3$fake_age)
sp_cover_f3_fp = split(sp_cover_f3, sp_cover_f3$f_p)


###### Sum the species for field level and all level
sp_cover_f2_field = sp_cover_f2 %>% group_by(Absolute_year,             
                                             Year, fake_year,                  
                                             Age, fake_age, Field) %>% 
                    summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
                                 sum, na.rm = T)
sp_cover_f2_field = split(sp_cover_f2_field, sp_cover_f2_field$Field)
sp_cover_f2_all = sp_cover_f2 %>% group_by(fake_age) %>% 
                  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
                  sum, na.rm = T)

###### Sum the species for field level and all level
sp_cover_f3_field = sp_cover_f3 %>% group_by(Absolute_year,             
                                             Year, fake_year,                  
                                             Age, fake_age, Field) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)
sp_cover_f3_field = split(sp_cover_f3_field, sp_cover_f3_field$Field)
sp_cover_f3_all = sp_cover_f3 %>% group_by(fake_age) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)

######## Filter abundant species for fit under plot/field and all field level
source('code/data preparation/code/function/get_consec_rep.R')
source('code/data preparation/code/function/data collation_function.R')
source('code/data preparation/code/function/fit_prepa_function.R')
source('code/data preparation/code/function/fit_prepa_function_2.R')
source('code/data preparation/code/function/fit_prepa_field_strict.R')
source('code/data preparation/code/function/fit_prepa_field_modest.R')

#### Plot level 
fit_stage_fp_alltime = lapply(sp_cover_f2_fp, fit_prepa, nati_per = 0.3,
                              intro_per = 0.15, estab_per = 0.3, domin_per = 0.3)
intro_sp = sapply(fit_stage_fp_alltime,
                  function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
save(fit_stage_fp_alltime,
     file = "code/data preparation/transformed data/fit_stage_fp_alltime.RData")

fit_stage_fp_alltime_same_0.3 = lapply(sp_cover_f2_fp, fit_prepa, nati_per = 0.3,
                                       intro_per = 0.3, estab_per = 0.3, domin_per = 0.3)
intro_sp_0.3 = sapply(fit_stage_fp_alltime_same_0.3,
                      function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
length(intro_sp_0.3[which(intro_sp_0.3 == 0)])/length(intro_sp_0.3)
save(fit_stage_fp_alltime_same_0.3,
     file = "code/data preparation/transformed data/fit_stage_fp_alltime_same_0.3.RData")

fit_stage_fp_alltime_same_0.25 = lapply(sp_cover_f2_fp, fit_prepa, nati_per = 0.25,
                                       intro_per = 0.25, estab_per = 0.25, domin_per = 0.25)
intro_sp_0.25 = sapply(fit_stage_fp_alltime_same_0.25,
                      function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
length(intro_sp_0.25[which(intro_sp_0.25 == 0)])/length(intro_sp_0.25)
save(fit_stage_fp_alltime_same_0.25,
     file = "code/data preparation/transformed data/fit_stage_fp_alltime_same_0.25.RData")

fit_stage_fp_alltime_same_0.20 = lapply(sp_cover_f2_fp, fit_prepa, nati_per = 0.20,
                                        intro_per = 0.20, estab_per = 0.20, domin_per = 0.20)
intro_sp_0.20 = sapply(fit_stage_fp_alltime_same_0.20,
                       function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
length(intro_sp_0.20[which(intro_sp_0.20 == 0)])/length(intro_sp_0.20)
save(fit_stage_fp_alltime_same_0.20,
     file = "code/data preparation/transformed data/fit_stage_fp_alltime_same_0.20.RData")

fit_stage_fp_alltime_same_0.15 = lapply(sp_cover_f2_fp, fit_prepa, nati_per = 0.15,
                                        intro_per = 0.15, estab_per = 0.15, domin_per = 0.15)
intro_sp_0.15 = sapply(fit_stage_fp_alltime_same_0.15,
                       function(x){length(strsplit(x[["stage_sp_name_f"]][3], ', ')[[1]])})
length(intro_sp_0.15[which(intro_sp_0.15 == 0)])/length(intro_sp_0.15)
save(fit_stage_fp_alltime_same_0.15,
     file = "code/data preparation/transformed data/fit_stage_fp_alltime_same_0.15.RData")


#### Using filtered time-series data to judge different invasion stage
fake_age_all = sort(unique(sp_cover_f3$fake_age))

fit_stage_fp_t_all_filter_1 = list()
for (i in (1:14)) {
  #i = 14
  seq_t = c(i:(i+12))
  if (i == 1) {
    sp_cover_f3_fp_f = lapply(sp_cover_f3_fp, function(x){
                              if (nrow(x) < 26) {
                                x = NULL} else {
                                  x = x
                                } 
    })
    sp_cover_f3_fp_f <- sp_cover_f3_fp_f[!sapply(sp_cover_f3_fp_f,is.null)]
    fit_stage_fp_t_all_filter_1[[i]] = lapply(sp_cover_f3_fp_f,
                                  fit_prepa,
                                  t = seq_t,
                                  seq_t_all = fake_age_all)
  } else {
    fit_stage_fp_t_all_filter_1[[i]] = lapply(sp_cover_f3_fp,
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
    sp_cover_f3_fp_f = lapply(sp_cover_f3_fp, function(x){
      if (nrow(x) < 26) {
        x = NULL} else {
          x = x
        } 
    })
    sp_cover_f3_fp_f <- sp_cover_f3_fp_f[!sapply(sp_cover_f3_fp_f,is.null)]
    fit_stage_fp_t_all_filter_2[[i]] = lapply(sp_cover_f3_fp_f,
                                              fit_prepa_2,
                                              sp_cover_origin = sp_cover_fp,
                                              t = seq_t,
                                              seq_t_all = fake_age_all)
  } else {
    fit_stage_fp_t_all_filter_2[[i]] = lapply(sp_cover_f3_fp,
                                              fit_prepa_2,
                                              sp_cover_origin = sp_cover_fp,
                                              t = seq_t,
                                              seq_t_all = fake_age_all)
  } 
}
save(fit_stage_fp_t_all_filter_2,
     file = './code/data preparation//transformed data/fit_stage_fp_t_all_filter_2.rdata')


#### Field level 
### filter some plots cannot cover all time span
fit_stage_field_all_time_strict = lapply(sp_cover_f2_field, fit_prepa_field_strict,
                                         seq_t_all = fake_age_all)
save(fit_stage_field_all_time_strict,
     file = "code/data preparation//transformed data/fit_stage_field_all_time_strict.RData")

### If established species appeared in the top ***% most abundant species,
### defined it as dominant species
fit_stage_field_all_time_modest = lapply(sp_cover_f2_field, fit_prepa_field_modest,
                                         seq_t_all = fake_age_all,
                                         domin_prop = 0.02)
save(fit_stage_field_all_time_modest,
     file = "code/data preparation/transformed data/fit_stage_field_all_time_modest.RData")

######## Time windows

fit_stage_field_t_all_filter_1 = list()
for (i in (1:14)) {
  seq_t = c(i:(i+12))
  if (i == 1) {
    sp_cover_f3_field_f = lapply(sp_cover_f3_field, function(x){
      if (nrow(x) < 26) {
        x = NULL} else {
          x = x
        } 
    })
    sp_cover_f3_field_f <- sp_cover_f3_field_f[!sapply(sp_cover_f3_field_f,is.null)]
    fit_stage_field_t_all_filter_1[[i]] = lapply(sp_cover_f3_field_f,
                                              fit_prepa_field_modest,
                                              t = seq_t,
                                              seq_t_all = fake_age_all,
                                              domin_prop = 0.02)
  } else {
    fit_stage_field_t_all_filter_1[[i]] = lapply(sp_cover_f3_field,
                                              fit_prepa_field_modest,
                                              t = seq_t,
                                              seq_t_all = fake_age_all,
                                              domin_prop = 0.02)
  } 
}

save(fit_stage_field_t_all_filter_1,
     file = './code/data preparation/transformed data/fit_stage_field_t_all_filter_1.rdata')


