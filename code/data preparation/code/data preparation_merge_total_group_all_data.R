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
sp_cover_f2 = sp_cover_f1 %>% filter(fake_age %in% c(1:51))
sp_cover_f2_fp = split(sp_cover_f2, sp_cover_f2$f_p)

sp_cover_f2_ages_1_35 = sp_cover_f1 %>% filter(fake_age < 36)
sp_cover_f2_ages_1_35_fp = split(sp_cover_f2_ages_1_35, sp_cover_f2_ages_1_35$f_p)

###### Merge species into 4 stages
sp_cover_f2_field = sp_cover_f2 %>% group_by(Absolute_year,             
                                             Year,                   
                                             Age,  Field) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)

sp_cover_f2_field = split(sp_cover_f2_field, sp_cover_f2_field$Field)

sp_cover_f2_all = sp_cover_f2 %>% group_by(Age) %>% 
  summarise_at(vars(Abutilon_theophrasti:Vitis_aestivalis),
               sum, na.rm = T)


######## Filter abundant species for fit under plot/field and all field level
source('code/data preparation/code/function/get_consec_rep.R') 
source('code/data preparation/code/function/data collation_function.R')
source('code/data preparation/code/function/fit_prepa_function_filtered.R')
source('code/data preparation/code/function/fit_prepa_function_top50.R')
#source('code/data preparation/code/function/fit_prepa_function_original_field.R')

### top 50
freq_sps = apply(sp_cover_f2_ages_1_35[,c(10:ncol(sp_cover_f2_ages_1_35))], 2, function(x){length(x[x>0])})
freq_sps_exotic = freq_sps[(trait %>% filter(Origin == 'Exotic'))$Species]
freq_sps_exotic = freq_sps_exotic[which(!is.na(freq_sps_exotic))]
freq_sps_native = freq_sps[(trait %>% filter(Origin == 'Native'))$Species] 
freq_sps_native = freq_sps_native[which(!is.na(freq_sps_native))]
freq_sps_exotic_top_25 = sort(freq_sps_exotic, decreasing = T)[1:25]
freq_sps_native_top_25 = sort(freq_sps_native, decreasing = T)[1:25]
freq_sps_top_50 = c(freq_sps_exotic_top_25, freq_sps_native_top_25)

sp_cover_f2_top50_ages1_35 = cbind(sp_cover_f2_ages_1_35[,c(1:9)],
                          sp_cover_f2_ages_1_35 %>% select(names(freq_sps_top_50)))
sp_cover_f2_top50_ages1_35_fp = split(sp_cover_f2_top50_ages1_35, sp_cover_f2_top50_ages1_35$f_p)
fake_age_all = sort(unique(sp_cover_f2_top50_ages1_35$fake_age))


sp_cover_f2_top50_ages1_35_merge_field = sp_cover_f2_top50_ages1_35 %>% group_by(
  #fake_year, 
  fake_age,
  Field) %>% 
  summarise_at(colnames(sp_cover_f2_top50_ages1_35)[10:ncol(sp_cover_f2_top50_ages1_35)],
               sum, na.rm = T)

sp_cover_f2_top50_ages1_35_merge_all = sp_cover_f2_top50_ages1_35 %>% 
  group_by(#fake_year, 
    fake_age) %>% 
  summarise_at(colnames(sp_cover_f2_top50_ages1_35)[10:ncol(sp_cover_f2_top50_ages1_35)],
               sum, na.rm = T)


fit_all_top50_ages1_35 = fit_prepa_top_merge(sp_cover_f2_top50_ages1_35_merge_all,
                                                     posi_sp_start = 2)

save(fit_all_top50_ages1_35, 
     file = "code/data preparation/transformed data/fit_all_top50_ages1_35.RData")
