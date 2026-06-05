#### Load libraries
rm(list = ls())
library(ape)
library(dplyr)
library(tidyr)
library(doParallel)
library(parallel)

setwd("D:/R projects/BSS_coexist")

### Calculate the original functional distance
load("results/fit_results/plot_ages1_35_top50_equal_interval_model_comparison/bh_partialb/inter_all_c_alltime.rdata")
tree = read.tree('data/original data/phylo_tree332.txt')
#tree_li_2015 = read.table('D:/????????????/DNH_BSS/0312EL_Darwin/20150312/invader1year_native20year/MPDab_20years.txt', 
#                         header = T)
plot(tree)
tree_dis = cophenetic(tree)

pd_dat = as.data.frame(tree_dis) %>% 
  mutate(species_i = rownames(tree_dis),
         .before = 'Erigeron_strigosus') %>%
  pivot_longer(!species_i, names_to = "species_j", values_to = "Phylo_dis")

pd_dat$sp_pair = paste(pd_dat$species_i,
                       pd_dat$species_j,
                       sep = '_')

### Load the community data 
##### Data loading
sp_cover = read.csv('data/original data/BSS_community_332.csv',
                    header = T)
field = read.csv('data/original data/FIELDS.csv',
                 header = T)
trait_1 = read.csv('data/original data/traits332.csv',
                   header = T)
trait = trait_1[-nrow(trait_1),]
colnames(trait_1)


#trait_supp = read.csv('data/original data/traits332_supp.csv',
#                      header = T)
trait = arrange(trait, trait$species_id)

sp_cover$f_p = paste(sp_cover$Field, sp_cover$Plot, sep = '_')
sp_cover = sp_cover %>% relocate(f_p, .after = Plot)
trait$species_id = paste('sp', trait$species_id, sep = '_')
sp_real = vector()
sp = colnames(sp_cover)[8:ncol(sp_cover)]

for (i in 1:length(sp)) {
  sp_real[i] = trait[trait$species_id == sp[i],]$Species
}

colnames(sp_cover)[8:ncol(sp_cover)] = sp_real # change the species name of data

trait.conti_mat = as.matrix(trait_1[,12:ncol(trait)])
length(which(is.na(trait.conti_mat)))/length(trait.conti_mat)

trait_noloss = trait_1 %>% filter(!is.na(Height..m.) &  !is.na(LDMC..g.g.) &
                                    !is.na(SLA..cm2.g.) & !is.na(Seedmass..g.1000.seeds.))



select_sps = sort(unique(inter_all_c_alltime$species_i))
setdiff(select_sps, trait_noloss$Species) 
# all selected species have real trait values, do not need do imputation!
setdiff(trait_noloss$Species, select_sps)
trait_selected = trait_1 %>% filter(Species %in% select_sps) %>% arrange(Species)
trait_noloss = trait_noloss %>% arrange(Species)

#### Functional distance for BSS
trait.conti = as.matrix(trait_noloss[,12:ncol(trait)])
trait.cate = as.matrix(trait_noloss[,7:11])

hist(trait.conti[,1])
hist(trait.conti[,2])
hist(trait.conti[,3])
hist(trait.conti[,4])

trait.conti.2 = trait.conti
trait.conti.2[,1] = log(trait.conti.2[,1])
trait.conti.2[,3] = log(trait.conti.2[,3])
trait.conti.2[,4] = log(trait.conti.2[,4])
hist(trait.conti.2[,1])
hist(trait.conti.2[,2])
hist(trait.conti.2[,3])
hist(trait.conti.2[,4])
colnames(trait.conti.2)[c(1,3,4)] = c('log(Height..m.)',
                                      'log(SLA..cm2.g.)',
                                      'log(Seedmass..g.1000.seeds.)')
str(trait.conti.2)
trait.conti.2 = apply(trait.conti.2, 2, function(x){as.numeric(x)})
trait.conti.2.scaled = apply(trait.conti.2, 2, scale)

#### Check the covary trends among different traits
#library(PerformanceAnalytics)
#chart.Correlation(trait.conti.2.scaled, histogram=TRUE, pch=19)
##  the max correlation effect is 0.68, which indicates there is little reduncy among
## different traits

fd.conti.mat = as.matrix(dist(trait.conti.2.scaled,
                             method = 'euclidean'))
colnames(fd.conti.mat) = trait_noloss$Species
rownames(fd.conti.mat) = trait_noloss$Species

save(fd.conti.mat, file = 'code/data preparation/Phylo_Func/fd.conti.mat.rdata')

fd_conti_dat = as.data.frame(fd.conti.mat) %>% 
  mutate(species_i = rownames(fd.conti.mat),
         .before = 'Acalypha_rhomboidea') %>%
  pivot_longer(!species_i,
               names_to = "species_j",
               values_to = "dist")

fd_conti_dat$sp_pair = paste(fd_conti_dat$species_i,
                            fd_conti_dat$species_j,
                            sep = '_')


# Calculate Gower distance for all traits
library(FD)
library(tidyr)
trait_mat = trait_noloss[,7:ncol(trait_noloss)]
trait_mat$Growth = as.factor(trait_mat$Growth)
trait_mat$Span = as.factor(trait_mat$Span)
trait_mat$Pollination = as.factor(trait_mat$Pollination)
trait_mat$Dispersal = as.factor(trait_mat$Dispersal)
trait_mat$Clonality = as.factor(trait_mat$Clonality)
trait_mat$Height..m. = as.numeric(trait_mat$Height..m.)
trait_mat$LDMC..g.g. = as.numeric(trait_mat$LDMC..g.g.)
trait_mat$SLA..cm2.g. = as.numeric(trait_mat$SLA..cm2.g.)
trait_mat$Seedmass..g.1000.seeds. = as.numeric(trait_mat$Seedmass..g.1000.seeds.)
rownames(trait_mat) = trait_noloss$Species

fd_gower_mat = as.matrix(gowdis(trait_mat))
save(fd_gower_mat, file = 'code/data preparation/Phylo_Func/fd_gower_mat.rdata')

fd_gower_dat = as.data.frame(fd_gower_mat) %>% 
  mutate(species_i = rownames(fd_gower_mat),
         .before = 'Acalypha_rhomboidea') %>%
  pivot_longer(!species_i,
               names_to = "species_j",
               values_to = "Gower_dis")

fd_gower_dat$sp_pair = paste(fd_gower_dat$species_i,
                             fd_gower_dat$species_j,
                             sep = '_')



# Calculate single traits distance
library(FD)
library(tidyr)
trait_mat = trait_noloss[,7:ncol(trait_noloss)]
trait_mat$Growth = as.factor(trait_mat$Growth)
trait_mat$Span = as.factor(trait_mat$Span)
trait_mat$Pollination = as.factor(trait_mat$Pollination)
trait_mat$Dispersal = as.factor(trait_mat$Dispersal)
trait_mat$Clonality = as.factor(trait_mat$Clonality)
trait_mat$Height..m. = as.numeric(trait_mat$Height..m.)
trait_mat$LDMC..g.g. = as.numeric(trait_mat$LDMC..g.g.)
trait_mat$SLA..cm2.g. = as.numeric(trait_mat$SLA..cm2.g.)
trait_mat$Seedmass..g.1000.seeds. = as.numeric(trait_mat$Seedmass..g.1000.seeds.)
rownames(trait_mat) = trait_noloss$Species

fd_mat_growth = as.matrix(gowdis(as.data.frame(trait_mat$Growth)))
rownames(fd_mat_growth) = trait_noloss$Species
colnames(fd_mat_growth) = trait_noloss$Species
fd_mat_span = as.matrix(gowdis(as.data.frame(trait_mat$Span)))
rownames(fd_mat_span) = trait_noloss$Species
colnames(fd_mat_span) = trait_noloss$Species
fd_mat_pollination = as.matrix(gowdis(as.data.frame(trait_mat$Pollination)))
rownames(fd_mat_pollination) = trait_noloss$Species
colnames(fd_mat_pollination) = trait_noloss$Species
fd_mat_dispersal = as.matrix(gowdis(as.data.frame(trait_mat$Dispersal)))
rownames(fd_mat_dispersal) = trait_noloss$Species
colnames(fd_mat_dispersal) = trait_noloss$Species
fd_mat_clonality = as.matrix(gowdis(as.data.frame(trait_mat$Clonality)))
rownames(fd_mat_clonality) = trait_noloss$Species
colnames(fd_mat_clonality) = trait_noloss$Species
fd_mat_height = as.matrix(gowdis(as.data.frame(trait_mat$Height..m.)))
rownames(fd_mat_height) = trait_noloss$Species
colnames(fd_mat_height) = trait_noloss$Species
fd_mat_ldmc = as.matrix(gowdis(as.data.frame(trait_mat$LDMC..g.g.)))
rownames(fd_mat_ldmc) = trait_noloss$Species
colnames(fd_mat_ldmc) = trait_noloss$Species
fd_mat_sla = as.matrix(gowdis(as.data.frame(trait_mat$SLA..cm2.g.)))
rownames(fd_mat_sla) = trait_noloss$Species
colnames(fd_mat_sla) = trait_noloss$Species
fd_mat_seedmass = as.matrix(gowdis(as.data.frame(trait_mat$Seedmass..g.1000.seeds.)))
rownames(fd_mat_seedmass) = trait_noloss$Species
colnames(fd_mat_seedmass) = trait_noloss$Species

fd_mat_list = list(fd_mat_growth, fd_mat_span, fd_mat_pollination,
                   fd_mat_dispersal, fd_mat_clonality, fd_mat_height,
                   fd_mat_ldmc, fd_mat_sla, fd_mat_seedmass
)

fd_dat_list = lapply(fd_mat_list, function(x) {
  as.data.frame(x) %>% 
    mutate(species_i = rownames(x),
           .before = rownames(x)[1]) %>%
    pivot_longer(!species_i, names_to = "species_j", values_to = "Gower_dis")
})

fd_dat_single = cbind(fd_dat_list[[1]][,c(1:2)],
                      sapply(fd_dat_list, function(x){
                        y = x$Gower_dis
                      }))
colnames(trait)
colnames(fd_dat_single)[3:ncol(fd_dat_single)] = c("growth","span",                  
                                                   "pollination","dispersal",
                                                   "clonality","height",
                                                   "ldmc","sla","seedmass")

fd_dat_single$sp_pair = paste(fd_dat_single$species_i,
                              fd_dat_single$species_j,
                              sep = '_')
fd_dat_single = fd_dat_single %>% relocate('sp_pair', .after = 'species_j')

fd_dat_noloss = fd_dat_single %>% left_join(fd_conti_dat[,c('dist', 'sp_pair')],
                                         by = 'sp_pair') %>% 
  left_join(fd_gower_dat[,c('Gower_dis', 'sp_pair')],
            by = 'sp_pair') %>% 
  rename(Multi_conti_traits = dist) %>% 
  rename(Multi_traits = Gower_dis)

#### Save the pd fd
save(pd_dat, file = 'code/data preparation/Phylo_Func/pd.rdata')
save(fd_dat_noloss, file = 'code/data preparation/Phylo_Func/fd_dat_noloss.rdata')


