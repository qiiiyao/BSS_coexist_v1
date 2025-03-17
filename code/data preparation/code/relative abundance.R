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
sapply(sp_cover_fp, function(x) {max(sort(unique(x$Absolute_year)))})
sapply(sp_cover_year, function(x) {sort(unique(x$f_p))})

##### Follow Li et al., (2015), condense the two year (half of field survey) to
##### one year and Follow Yin et al., (2022),, our study restricted to 3-51 age
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

sp_cover_fakeage1 = sp_cover_f1 %>% filter(fake_age == 1)
sp_cover_f2 = sp_cover_f1 %>% filter(fake_age %in% seq(1, 51, 2) &
                                       f_p %in% sp_cover_fakeage1$f_p)

sp_cover_f2 = sp_cover_f2 %>% relocate(fake_year, .before = Age)
sp_cover_f2 = sp_cover_f2 %>% relocate(fake_age, .before = Field)
save(sp_cover_f2, file = 'code/data preparation/transformed data/sp_cover_f2.rdata')
save(sp_cover_f1, file = 'code/data preparation/transformed data/sp_cover_f1.rdata')

#### Not delete other plots for including age1
sp_cover_f3 = sp_cover_f1 %>% filter(fake_age %in% seq(1, 51, 2))
sp_cover_f3 = sp_cover_f3 %>% relocate(fake_year, .before = Age)
sp_cover_f3 = sp_cover_f3 %>% relocate(fake_age, .before = Field)
sp_cover_f3_y = split(sp_cover_f3, sp_cover_f3$fake_age)
sapply(sp_cover_f3_y, function(x) {sort(unique(x$f_p))})
unique(sp_cover_f3$fake_age)
sp_cover_f3_fp = split(sp_cover_f3, sp_cover_f3$f_p)
save(sp_cover_f3, file = 'code/data preparation/transformed data/sp_cover_f3.rdata')

###### f1 Sum the species for field level and all level
sp_cover_f1_field = sp_cover_f1 %>% group_by(Absolute_year,             
                                             Year, fake_year,                  
                                             Age, fake_age, Field) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)
#sp_cover_f1_field = split(sp_cover_f1_field, sp_cover_f1_field$Field)
sp_cover_f1_all = sp_cover_f1 %>% group_by(fake_age) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)
save(sp_cover_f1_field, file = 'code/data preparation/transformed data/sp_cover_f1_field.rdata')
save(sp_cover_f1_all, file = 'code/data preparation/transformed data/sp_cover_f1_all.rdata')

###### f2 Sum the species for field level and all level
sp_cover_f2_field = sp_cover_f2 %>% group_by(Absolute_year,             
                                             Year, fake_year,                  
                                             Age, fake_age, Field) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)
#sp_cover_f2_field = split(sp_cover_f2_field, sp_cover_f2_field$Field)
sp_cover_f2_all = sp_cover_f2 %>% group_by(fake_age) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)
save(sp_cover_f2_field, file = 'code/data preparation/transformed data/sp_cover_f2_field.rdata')
save(sp_cover_f2_all, file = 'code/data preparation/transformed data/sp_cover_f2_all.rdata')

###### f3 Sum the species for field level and all level
sp_cover_f3_field = sp_cover_f3 %>% group_by(Absolute_year,             
                                             Year, fake_year,                  
                                             Age, fake_age, Field) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)
#sp_cover_f3_field = split(sp_cover_f3_field, sp_cover_f3_field$Field)
sp_cover_f3_all = sp_cover_f3 %>% group_by(fake_age) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)
save(sp_cover_f3_field, file = 'code/data preparation/transformed data/sp_cover_f3_field.rdata')
save(sp_cover_f3_all, file = 'code/data preparation/transformed data/sp_cover_f3_all.rdata')

#### Create the relative abundance for weighted mean calculating in the plot level
rm(list = ls())
library(dplyr)
library(data.table)
library(vegan)
library(doParallel)

load('code/data preparation/transformed data/sp_cover_f1.rdata')
load('code/data preparation/transformed data/sp_cover_f2.rdata')
load('code/data preparation/transformed data/sp_cover_f3.rdata')

fake_age_f2 = sort(unique(sp_cover_f2$fake_age))
sp_real = colnames(sp_cover_f2)[10:ncol(sp_cover_f2)]

#### For Age all time fit ####
### Real time data
sp_cover_f1_mean_fp = sp_cover_f1 %>%
  group_by(f_p) %>% 
  summarise_at(vars(all_of(sp_real)),
               mean, na.rm = T)

sp_racover_f1_mean_fp = cbind(f_p = sp_cover_f1_mean_fp$f_p,
                              decostand(sp_cover_f1_mean_fp[,
                                                            c(2:ncol(sp_cover_f1_mean_fp))],
                                        method = 'total',
                                        MARGIN = 1))
sp_racover_f1_mean_fp_alltime_raw = rbindlist(apply(sp_racover_f1_mean_fp, 1, function(x){
  y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
}))
colnames(sp_racover_f1_mean_fp_alltime_raw) = c('f_p', 'species', 'ra_m')
sp_racover_f1_mean_fp_alltime_raw$f_p_species = paste(sp_racover_f1_mean_fp_alltime_raw$f_p,
                                                  sp_racover_f1_mean_fp_alltime_raw$species,
                                                  sep = '_')


sp_racover_f1_mean_fp_alltime_raw = sp_racover_f1_mean_fp_alltime_raw %>% relocate(f_p_species,
                                                                               .after = species)
sp_racover_f1_mean_fp_alltime_raw$ra_m = as.numeric(sp_racover_f1_mean_fp_alltime_raw$ra_m)
save(sp_racover_f1_mean_fp_alltime_raw,
     file = 'code/data preparation/transformed data/sp_racover_f1_mean_fp_alltime_raw.rdata')

#### For Age 1-51 fit ####
### Real time data
sp_cover_f1_mean_fp = sp_cover_f1 %>%
  filter(fake_age %in% seq(1, 51, 1)) %>% 
  group_by(f_p) %>% 
  summarise_at(vars(all_of(sp_real)),
               mean, na.rm = T)

sp_racover_f1_mean_fp = cbind(f_p = sp_cover_f1_mean_fp$f_p,
                             decostand(sp_cover_f1_mean_fp[,
                                                          c(2:ncol(sp_cover_f1_mean_fp))],
                                       method = 'total',
                                       MARGIN = 1))
sp_racover_f1_mean_fp_alltime = rbindlist(apply(sp_racover_f1_mean_fp, 1, function(x){
  y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
}))
colnames(sp_racover_f1_mean_fp_alltime) = c('f_p', 'species', 'ra_m')
sp_racover_f1_mean_fp_alltime$f_p_species = paste(sp_racover_f1_mean_fp_alltime$f_p,
                                                 sp_racover_f1_mean_fp_alltime$species,
                                                 sep = '_')
### Fitted time data
sp_cover_f2_mean_fp = sp_cover_f2 %>%
  group_by(f_p) %>% 
  summarise_at(vars(all_of(sp_real)),
               mean, na.rm = T)

sp_racover_f2_mean_fp = cbind(f_p = sp_cover_f2_mean_fp$f_p,
                             decostand(sp_cover_f2_mean_fp[,
                                                          c(2:ncol(sp_cover_f2_mean_fp))],
                                       method = 'total',
                                       MARGIN = 1))
sp_racover_f2_mean_fp_alltime = rbindlist(apply(sp_racover_f2_mean_fp, 1, function(x){
  y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
}))
colnames(sp_racover_f2_mean_fp_alltime) = c('f_p', 'species', 'ra_m')
sp_racover_f2_mean_fp_alltime$f_p_species = paste(sp_racover_f2_mean_fp_alltime$f_p,
                                                 sp_racover_f2_mean_fp_alltime$species,
                                                 sep = '_')

## merge real time and fitted time
sp_racover_f1_2_mean_fp_alltime = sp_racover_f1_mean_fp_alltime %>% left_join(sp_racover_f2_mean_fp_alltime[,c(3,4)],
                                                                            by = join_by(f_p_species))
colnames(sp_racover_f1_2_mean_fp_alltime)[c(3,5)] = c('ra_m_real_t', 'ra_m_fit_t')
sp_racover_f1_2_mean_fp_alltime = sp_racover_f1_2_mean_fp_alltime %>% relocate(f_p_species,
                                                                             .after = species)
sp_racover_f1_2_mean_fp_alltime$ra_m_real_t = as.numeric(sp_racover_f1_2_mean_fp_alltime$ra_m_real_t)
sp_racover_f1_2_mean_fp_alltime$ra_m_fit_t = as.numeric(sp_racover_f1_2_mean_fp_alltime$ra_m_fit_t)
save(sp_racover_f1_2_mean_fp_alltime,
     file = 'code/data preparation/transformed data/sp_racover_f1_2_mean_fp_alltime.rdata')

#### For early succession 1-26 fit ####
### Real time data
sp_cover_f1_mean_fp = sp_cover_f1 %>%
  filter(fake_age < 26) %>% 
  group_by(f_p) %>% 
  summarise_at(vars(all_of(sp_real)),
               mean, na.rm = T)

sp_racover_f1_mean_fp = cbind(f_p = sp_cover_f1_mean_fp$f_p,
                              decostand(sp_cover_f1_mean_fp[,
                                                            c(2:ncol(sp_cover_f1_mean_fp))],
                                        method = 'total',
                                        MARGIN = 1))
sp_racover_f1_mean_fp_early_suc = rbindlist(apply(sp_racover_f1_mean_fp, 1, function(x){
  y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
}))

colnames(sp_racover_f1_mean_fp_early_suc) = c('f_p', 'species', 'ra_m')
sp_racover_f1_mean_fp_early_suc$f_p_species = paste(sp_racover_f1_mean_fp_early_suc$f_p,
                                                  sp_racover_f1_mean_fp_early_suc$species,
                                                  sep = '_')
colnames(sp_racover_f1_mean_fp_early_suc)[3] = c('ra_m_real_t')
sp_racover_f1_mean_fp_early_suc = sp_racover_f1_mean_fp_early_suc %>% relocate(f_p_species,
                                                                               .after = species)
sp_racover_f1_mean_fp_early_suc$ra_m_real_t = as.numeric(sp_racover_f1_mean_fp_early_suc$ra_m_real_t)
save(sp_racover_f1_mean_fp_early_suc,
     file = 'code/data preparation/transformed data/sp_racover_f1_mean_fp_early_suc.rdata')

#### For early succession 1-35 fit ####
### Real time data
sp_cover_f1_mean_fp = sp_cover_f1 %>%
  filter(fake_age < 36) %>% 
  group_by(f_p) %>% 
  summarise_at(vars(all_of(sp_real)),
               mean, na.rm = T)

sp_racover_f1_mean_fp = cbind(f_p = sp_cover_f1_mean_fp$f_p,
                              decostand(sp_cover_f1_mean_fp[,
                                                            c(2:ncol(sp_cover_f1_mean_fp))],
                                        method = 'total',
                                        MARGIN = 1))
sp_racover_f1_mean_fp_early_suc_1_35 = rbindlist(apply(sp_racover_f1_mean_fp, 1, function(x){
  y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
}))

colnames(sp_racover_f1_mean_fp_early_suc_1_35) = c('f_p', 'species', 'ra_m')
sp_racover_f1_mean_fp_early_suc_1_35$f_p_species = paste(sp_racover_f1_mean_fp_early_suc_1_35$f_p,
                                                    sp_racover_f1_mean_fp_early_suc_1_35$species,
                                                    sep = '_')
colnames(sp_racover_f1_mean_fp_early_suc_1_35)[3] = c('ra_m_real_t')
sp_racover_f1_mean_fp_early_suc_1_35 = sp_racover_f1_mean_fp_early_suc_1_35 %>% relocate(f_p_species,
                                                                               .after = species)
sp_racover_f1_mean_fp_early_suc_1_35$ra_m_real_t = as.numeric(sp_racover_f1_mean_fp_early_suc_1_35$ra_m_real_t)
save(sp_racover_f1_mean_fp_early_suc_1_35,
     file = 'code/data preparation/transformed data/sp_racover_f1_mean_fp_early_suc_1_35.rdata')


#### For ages 1-31 fit ####
load('code/data preparation/transformed data/sp_cover_f1.rdata')
library(dplyr)
library(vegan)
### Real time data
sp_cover_f1_mean_fp = sp_cover_f1 %>%
  filter(fake_age < 32) %>% 
  group_by(f_p) %>% 
  summarise_at(vars(all_of(colnames(sp_cover_f1)[10:ncol(sp_cover_f1)])),
               mean, na.rm = T)

sp_racover_f1_mean_fp = cbind(f_p = sp_cover_f1_mean_fp$f_p,
                              decostand(sp_cover_f1_mean_fp[,
                                                            c(2:ncol(sp_cover_f1_mean_fp))],
                                        method = 'total',
                                        MARGIN = 1))
sp_racover_f1_mean_fp_ages1_31 = rbindlist(apply(sp_racover_f1_mean_fp, 1, function(x){
  y = as.data.frame(cbind(x[1], colnames(sp_cover_f1)[10:ncol(sp_cover_f1)],
                          x[2:291]))
}))

colnames(sp_racover_f1_mean_fp_ages1_31) = c('f_p', 'species', 'ra_m')
sp_racover_f1_mean_fp_ages1_31$f_p_species = paste(sp_racover_f1_mean_fp_ages1_31$f_p,
                                                         sp_racover_f1_mean_fp_ages1_31$species,
                                                         sep = '_')
colnames(sp_racover_f1_mean_fp_ages1_31)[3] = c('ra_m_real_t')
sp_racover_f1_mean_fp_ages1_31 = sp_racover_f1_mean_fp_ages1_31 %>% relocate(f_p_species,
                                                                                         .after = species)
sp_racover_f1_mean_fp_ages1_31$ra_m_real_t = as.numeric(sp_racover_f1_mean_fp_ages1_31$ra_m_real_t)
save(sp_racover_f1_mean_fp_ages1_31,
     file = 'code/data preparation/transformed data/sp_racover_f1_mean_fp_ages1_31.rdata')


#### Calculating mean relative cover for moving windows
fake_age_f3 = sort(unique(sp_cover_f3$fake_age))
sp_racover_f1_3_mean_fp_t_l = foreach(i = 1:14, .packages = c('dplyr', 'vegan',
                                                         'data.table')) %dopar% {
                                                           #i = 1
                                seq_t = seq(fake_age_f3[i], fake_age_f3[i+12],1)
                                sp_cover_f1_mean_fp_t = sp_cover_f1 %>%
                                                             filter(fake_age %in% seq_t) %>% 
                                                             group_by(f_p) %>% 
                                                             summarise_at(vars(all_of(sp_real)),
                                                                          mean, na.rm = T)
                                                           
                                sp_racover_f1_mean_fp_t = cbind(f_p = sp_cover_f1_mean_fp_t$f_p,
                                                               decostand(sp_cover_f1_mean_fp_t[,
                                                               c(2:ncol(sp_cover_f1_mean_fp_t))],
                                                              method = 'total',
                                                              MARGIN = 1))
                                sp_racover_f1_mean_fp_t = rbindlist(apply(sp_racover_f1_mean_fp_t, 1, function(x){
                                y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
                                          }))
                                colnames(sp_racover_f1_mean_fp_t) = c('f_p', 'species', 'ra_m')
                                sp_racover_f1_mean_fp_t$f_p_species = paste(sp_racover_f1_mean_fp_t$f_p,
                                                                      sp_racover_f1_mean_fp_t$species,
                                                                      sep = '_')
                                sp_racover_f1_mean_fp_t$ra_m = as.numeric(sp_racover_f1_mean_fp_t$ra_m)
                                
                                ### fitted time
                                sp_cover_f3_mean_fp_t = sp_cover_f3 %>%
                                  filter(fake_age %in% seq_t) %>% 
                                  group_by(f_p) %>% 
                                  summarise_at(vars(all_of(sp_real)),
                                               mean, na.rm = T)
                                
                                sp_racover_f3_mean_fp_t = cbind(f_p = sp_cover_f3_mean_fp_t$f_p,
                                                               decostand(sp_cover_f3_mean_fp_t[,
                                                                                              c(2:ncol(sp_cover_f3_mean_fp_t))],
                                                                         method = 'total',
                                                                         MARGIN = 1))
                                sp_racover_f3_mean_fp_t = rbindlist(apply(sp_racover_f3_mean_fp_t, 1, function(x){
                                  y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
                                }))
                                colnames(sp_racover_f3_mean_fp_t) = c('f_p', 'species', 'ra_m')
                                sp_racover_f3_mean_fp_t$f_p_species = paste(sp_racover_f3_mean_fp_t$f_p,
                                                                           sp_racover_f3_mean_fp_t$species,
                                                                           sep = '_')
                                sp_racover_f3_mean_fp_t$ra_m = as.numeric(sp_racover_f3_mean_fp_t$ra_m)
                                
                                sp_racover_f1_3_mean_fp_t = sp_racover_f1_mean_fp_t %>%
                                             left_join(sp_racover_f3_mean_fp_t[,c(3,4)],
                                            by = join_by(f_p_species))
                                colnames(sp_racover_f1_3_mean_fp_t)[c(3,5)] = c('ra_m_real_t', 'ra_m_fit_t')
                                sp_racover_f1_3_mean_fp_t = sp_racover_f1_3_mean_fp_t %>% relocate(f_p_species,
                                                            .after = species)
                                
                                return(sp_racover_f1_3_mean_fp_t)
                                                           
                                }

save(sp_racover_f1_3_mean_fp_t_l,
     file = 'code/data preparation/transformed data/sp_racover_f1_3_mean_fp_t_l.rdata')


#### Create the relative abundance for weighted mean calculating in the field level
load('code/data preparation/transformed data/sp_cover_f1_field.rdata')
load('code/data preparation/transformed data/sp_cover_f2_field.rdata')
load('code/data preparation/transformed data/sp_cover_f3_field.rdata')

fake_age_f2 = sort(unique(sp_cover_f2$fake_age))
### Real time data
sp_cover_f1_mean_field = sp_cover_f1_field %>%
  group_by(Field) %>% 
  summarise_at(vars(all_of(sp_real)),
               mean, na.rm = T)

sp_racover_f1_mean_field = cbind(field = sp_cover_f1_mean_field$Field,
                              decostand(sp_cover_f1_mean_field[,
                                                            c(2:ncol(sp_cover_f1_mean_field))],
                                        method = 'total',
                                        MARGIN = 1))
sp_racover_f1_mean_field_alltime = rbindlist(apply(sp_racover_f1_mean_field, 1, function(x){
  y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
}))

colnames(sp_racover_f1_mean_field_alltime) = c('field', 'species', 'ra_m')
sp_racover_f1_mean_field_alltime$field_species = paste(sp_racover_f1_mean_field_alltime$field,
                                                  sp_racover_f1_mean_field_alltime$species,
                                                  sep = '_')

sp_racover_f1_mean_field_alltime = sp_racover_f1_mean_field_alltime %>% relocate(field_species,
                                                                                     .after = species)
sp_racover_f1_mean_field_alltime$ra_m = as.numeric(sp_racover_f1_mean_field_alltime$ra_m)
save(sp_racover_f1_mean_field_alltime,
     file = 'code/data preparation/transformed data/sp_racover_f1_mean_field_alltime.rdata')

### Fitted time data
sp_cover_f2_mean_field = sp_cover_f2 %>%
  group_by(Field) %>% 
  summarise_at(vars(all_of(sp_real)),
               mean, na.rm = T)

sp_racover_f2_mean_field = cbind(field = sp_cover_f2_mean_field$Field,
                              decostand(sp_cover_f2_mean_field[,
                                                            c(2:ncol(sp_cover_f2_mean_field))],
                                        method = 'total',
                                        MARGIN = 1))
sp_racover_f2_mean_field_alltime = rbindlist(apply(sp_racover_f2_mean_field, 1, function(x){
  y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
}))
colnames(sp_racover_f2_mean_field_alltime) = c('f_p', 'species', 'ra_m')
sp_racover_f2_mean_field_alltime$field_species = paste(sp_racover_f2_mean_field_alltime$field,
                                                  sp_racover_f2_mean_field_alltime$species,
                                                  sep = '_')
#### Merge real time and fitted time
sp_racover_f1_2_mean_field_alltime = sp_racover_f1_mean_field_alltime %>% 
                                     left_join(sp_racover_f2_mean_field_alltime[,c(3,4)],
                                               by = join_by(field_species))
colnames(sp_racover_f1_2_mean_field_alltime)[c(3,5)] = c('ra_m_real_t', 'ra_m_fit_t')
sp_racover_f1_2_mean_field_alltime = sp_racover_f1_2_mean_field_alltime %>% relocate(field_species,
                                                                               .after = species)
sp_racover_f1_2_mean_field_alltime$ra_m_real_t = as.numeric(sp_racover_f1_2_mean_field_alltime$ra_m_real_t)
sp_racover_f1_2_mean_field_alltime$ra_m_fit_t = as.numeric(sp_racover_f1_2_mean_field_alltime$ra_m_fit_t)
save(sp_racover_f1_2_mean_field_alltime,
     file = 'code/data preparation/transformed data/sp_racover_f1_2_mean_field_alltime.rdata')

#### Calculating mean relative cover for moving windows
fake_age_f3 = sort(unique(sp_cover_f3_field$fake_age))
sp_racover_f1_3_mean_field_t_l = foreach(i = 1:14, .packages = c('dplyr', 'vegan',
                                                              'data.table')) %dopar% {
                                                                #i = 1
                                    seq_t = seq(fake_age_f3[i], fake_age_f3[i+12],1)
                                    sp_cover_f1_mean_field_t = sp_cover_f1_field %>%
                                    filter(fake_age %in% seq_t) %>% 
                                    group_by(Field) %>% 
                                    summarise_at(vars(all_of(sp_real)),
                                    mean, na.rm = T)
                                    
                                   sp_racover_f1_mean_field_t = cbind(field = sp_cover_f1_mean_field_t$Field,
                                                                      decostand(sp_cover_f1_mean_field_t[,
                                                                                c(2:ncol(sp_cover_f1_mean_field_t))],
                                                                                method = 'total',
                                                                                MARGIN = 1))
                                   sp_racover_f1_mean_field_t = rbindlist(apply(sp_racover_f1_mean_field_t,
                                                                                1,
                                                                                function(x){
                                         y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
                                        }))
                                   colnames(sp_racover_f1_mean_field_t) = c('field', 'species', 'ra_m')
                                   sp_racover_f1_mean_field_t$field_species = paste(sp_racover_f1_mean_field_t$field,
                                                                                  sep = '_')
                                   sp_racover_f1_mean_field_t$ra_m = as.numeric(sp_racover_f1_mean_field_t$ra_m)
                                                                
                                    ### fitted time
                                    sp_cover_f3_mean_field_t = sp_cover_f3_field %>%
                                    filter(fake_age %in% seq_t) %>% 
                                    group_by(Field) %>% 
                                    summarise_at(vars(all_of(sp_real)),
                                    mean, na.rm = T)
                                    
                                    sp_racover_f3_mean_field_t = cbind(field = sp_cover_f3_mean_field_t$Field,
                                                                  decostand(sp_cover_f3_mean_field_t[,
                                                                  c(2:ncol(sp_cover_f3_mean_field_t))],
                                                                  method = 'total',
                                                                  MARGIN = 1))
                                    sp_racover_f3_mean_field_t = rbindlist(apply(sp_racover_f3_mean_field_t,
                                                                                 1,
                                                                                 function(x){
                                                      y = as.data.frame(cbind(x[1], sp_real, x[2:291]))
                                                      }))
                                    colnames(sp_racover_f3_mean_field_t) = c('field', 'species', 'ra_m')
                                    sp_racover_f3_mean_field_t$field_species = paste(sp_racover_f3_mean_field_t$field,
                                                                               sp_racover_f3_mean_field_t$species,
                                                                               sep = '_')
                                    sp_racover_f3_mean_field_t$ra_m = as.numeric(sp_racover_f3_mean_field_t$ra_m)
                                    
                                    ## Merge real time and fitted time data                            
                                    sp_racover_f1_3_mean_field_t = sp_racover_f1_mean_field_t %>%
                                                                  left_join(sp_racover_f3_mean_field_t[,c(3,4)],
                                                                            by = join_by(field_species))
                                    colnames(sp_racover_f1_3_mean_field_t)[c(3,5)] = c('ra_m_real_t', 'ra_m_fit_t')
                                    sp_racover_f1_3_mean_field_t = sp_racover_f1_3_mean_field_t %>% 
                                                                   relocate(field_species,
                                                                            .after = species)
                                                                
                                     return(sp_racover_f1_3_mean_field_t)
                                                                
                                                              }

save(sp_racover_f1_3_mean_field_t_l,
     file = 'code/data preparation/transformed data/sp_racover_f1_3_mean_field_t_l.rdata')
