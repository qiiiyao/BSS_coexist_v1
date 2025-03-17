#### Load libraries
rm(list = ls())
library(ape)
library(dplyr)
library(tidyr)
library(doParallel)
library(parallel)

### Calculate the original functional distance
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


#### Functional distance for BSS
trait.conti = as.matrix(trait[,12:ncol(trait)])
trait.cate = as.matrix(trait[,7:11])

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
trait.conti.2.scaled = apply(trait.conti.2, 2, 
                             scale)

#### Check the covary trends among different traits
library(PerformanceAnalytics)
chart.Correlation(trait.conti.2.scaled, histogram=TRUE, pch=19)
##  the max correlation effect is 0.68, which indicates there is little reduncy among
## different traits

### Impute the miss values
library(Rphylopars)

pre = phylopars(trait_data = cbind(species = trait$Species,
                                   as_tibble(trait.conti.2.scaled)),
                tree = tree,
                pheno_error = TRUE,
                phylo_correlated = TRUE,
                pheno_correlated = TRUE)
trait.conti.2_fake = cbind(species = row.names(pre$anc_recon)[1:nrow(trait.conti.2.scaled)]
  ,as_tibble(pre$anc_recon[1:nrow(trait.conti.2.scaled),])) %>% 
  arrange(species)

pc.dist.mat = as.matrix(dist(trait.conti.2_fake[,2:ncol(trait.conti.2_fake)],
                             method = 'euclidean'))
colnames(pc.dist.mat) = trait.conti.2_fake$species
row.names(pc.dist.mat) = trait.conti.2_fake$species

save(pc.dist.mat, file = 'code/Phylo_Func/pc.dist.mat.rdata')

pc_dist_dat = as.data.frame(pc.dist.mat) %>% 
  mutate(species_i = rownames(pc.dist.mat),
         .before = 'Abutilon_theophrasti') %>%
  pivot_longer(!species_i,
               names_to = "species_j",
               values_to = "dist")

pc_dist_dat$sp_pair = paste(pc_dist_dat$species_i,
                            pc_dist_dat$species_j,
                            sep = '_')


# Calculate Gower distance for all traits
library(FD)
library(tidyr)
trait_mat = trait[,7:ncol(trait)]
trait_mat$Growth = as.factor(trait_mat$Growth)
trait_mat$Span = as.factor(trait_mat$Span)
trait_mat$Pollination = as.factor(trait_mat$Pollination)
trait_mat$Dispersal = as.factor(trait_mat$Dispersal)
trait_mat$Clonality = as.factor(trait_mat$Clonality)
trait_mat$Height..m. = as.numeric(trait_mat$Height..m.)
trait_mat$LDMC..g.g. = as.numeric(trait_mat$LDMC..g.g.)
trait_mat$SLA..cm2.g. = as.numeric(trait_mat$SLA..cm2.g.)
trait_mat$Seedmass..g.1000.seeds. = as.numeric(trait_mat$Seedmass..g.1000.seeds.)
rownames(trait_mat) = trait$Species

fd_gower_mat = as.matrix(gowdis(trait_mat))
save(fd_gower_mat, file = 'code/Phylo_Func/fd_gower_mat.rdata')

fd_gower_dat = as.data.frame(fd_gower_mat) %>% 
  mutate(species_i = rownames(fd_gower_mat),
         .before = 'Abutilon_theophrasti') %>%
  pivot_longer(!species_i,
               names_to = "species_j",
               values_to = "Gower_dis")

fd_gower_dat$sp_pair = paste(fd_gower_dat$species_i,
                             fd_gower_dat$species_j,
                             sep = '_')



# Calculate single traits distance
library(FD)
library(tidyr)
trait_mat = trait[,7:ncol(trait)]
trait_mat$Growth = as.factor(trait_mat$Growth)
trait_mat$Span = as.factor(trait_mat$Span)
trait_mat$Pollination = as.factor(trait_mat$Pollination)
trait_mat$Dispersal = as.factor(trait_mat$Dispersal)
trait_mat$Clonality = as.factor(trait_mat$Clonality)
trait_mat$Height..m. = as.numeric(trait_mat$Height..m.)
trait_mat$LDMC..g.g. = as.numeric(trait_mat$LDMC..g.g.)
trait_mat$SLA..cm2.g. = as.numeric(trait_mat$SLA..cm2.g.)
trait_mat$Seedmass..g.1000.seeds. = as.numeric(trait_mat$Seedmass..g.1000.seeds.)
rownames(trait_mat) = trait$Species

fd_mat_growth = as.matrix(gowdis(as.data.frame(trait_mat$Growth)))
rownames(fd_mat_growth) = trait$Species
colnames(fd_mat_growth) = trait$Species
fd_mat_span = as.matrix(gowdis(as.data.frame(trait_mat$Span)))
rownames(fd_mat_span) = trait$Species
colnames(fd_mat_span) = trait$Species
fd_mat_pollination = as.matrix(gowdis(as.data.frame(trait_mat$Pollination)))
rownames(fd_mat_pollination) = trait$Species
colnames(fd_mat_pollination) = trait$Species
fd_mat_dispersal = as.matrix(gowdis(as.data.frame(trait_mat$Dispersal)))
rownames(fd_mat_dispersal) = trait$Species
colnames(fd_mat_dispersal) = trait$Species
fd_mat_clonality = as.matrix(gowdis(as.data.frame(trait_mat$Clonality)))
rownames(fd_mat_clonality) = trait$Species
colnames(fd_mat_clonality) = trait$Species
fd_mat_height = as.matrix(gowdis(as.data.frame(trait_mat$Height..m.)))
rownames(fd_mat_height) = trait$Species
colnames(fd_mat_height) = trait$Species
fd_mat_ldmc = as.matrix(gowdis(as.data.frame(trait_mat$LDMC..g.g.)))
rownames(fd_mat_ldmc) = trait$Species
colnames(fd_mat_ldmc) = trait$Species
fd_mat_sla = as.matrix(gowdis(as.data.frame(trait_mat$SLA..cm2.g.)))
rownames(fd_mat_sla) = trait$Species
colnames(fd_mat_sla) = trait$Species
fd_mat_seedmass = as.matrix(gowdis(as.data.frame(trait_mat$Seedmass..g.1000.seeds.)))
rownames(fd_mat_seedmass) = trait$Species
colnames(fd_mat_seedmass) = trait$Species

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

fd_dat_all = fd_dat_single %>% left_join(pc_dist_dat[,c('dist', 'sp_pair')],
                                         by = 'sp_pair') %>% 
  left_join(fd_gower_dat[,c('Gower_dis', 'sp_pair')],
            by = 'sp_pair') %>% 
  rename(Multi_conti_traits = dist) %>% 
  rename(Multi_traits = Gower_dis)

#### Save the pd fd
save(pd_dat, file = 'code/Phylo_Func/pd.rdata')
save(fd_dat_all, file = 'code/Phylo_Func/fd_dat_all.rdata')


