###spatio-temporal data
###data sorting
library(spacetime)
library(sp)
library(foreign)

##### Data loading
sp_cover = read.csv('data/original data/BSS_community_332.csv',
                    header = T)
field = read.csv("data/original data/FIELDS.csv")

trait = read.csv('D:/R projects/BSS/data/original data/traits332.csv',
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

C3 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'C3',]$FIELDID)
C4 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'C4',]$FIELDID)
C5 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'C5',]$FIELDID)
C6 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'C6',]$FIELDID)
C7 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'C7',]$FIELDID)
D1 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'D1',]$FIELDID)
D2 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'D2',]$FIELDID)
D3 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'D3',]$FIELDID)
E1 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'E1',]$FIELDID)
E2 = sp_cover %>% filter(Field %in% field[field$FIELDNAME == 'E2',]$FIELDID)

### Exclude trees for fitting
sp_cover = cbind(sp_cover[,c(1:7)],
                 sp_cover %>% select(intersect(trait[trait$Span %in% c("Annual", "Perennial", "Biennial"),]$Species, colnames(sp_cover)[8:ncol(sp_cover)])))
sort(unique(sp_cover$Absolute_year))
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

###data of 50 most frequent species 
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
save(sp_cover_f2_top50, file = 'D:/R projects/BSS/code/Yin_ele_2022/sp_cover_f2_top50.rdata')

require(tidyr)
sp_cover_f2_top50 = sp_cover_f2_top50 %>% pivot_longer(cols = Lonicera_japonica:Aster_ericoides,
                                   names_to = 'Species',
                                   values_to = 'COVER') %>% 
                                   left_join(trait, by = 'Species') %>% 
                                   left_join(field, by = join_by(Field == FIELDID))

sp_cover_f2_top50.annual = subset(sp_cover_f2_top50,sp_cover_f2_top50$Span=="Annual")
length(unique(sp_cover_f2_top50.annual$Species))  ##4 annual sp
sp_cover_f2_top50.biennial=subset(sp_cover_f2_top50,sp_cover_f2_top50$Span=="Biennial")
length(unique(sp_cover_f2_top50.biennial$Species))  ###7 biennial sp 
sp_cover_f2_top50.perennial=subset(sp_cover_f2_top50,sp_cover_f2_top50$Span=="Perennial")
length(unique(sp_cover_f2_top50.perennial$Species))  ###39 perennial sp 

Field=c(rep("C3",48),rep("C4",48),rep("C5",48),rep("C6",48),rep("C7",48),rep("D1",48),rep("D2",48),rep("D3",48),rep("E1",48),rep("E2",48))
Plotid=rep(seq(1:48),10)
Fieldnum=c(rep(1,48),rep(2,48),rep(3,48),rep(4,48),rep(5,48),rep(6,48),rep(7,48),rep(8,48),rep(9,48),rep(10,48))

###creat 0-1 data of each species
Fieldname=unique(sp_cover_f2_top50$FIELDNAME)
Plot=sort(unique(sp_cover_f2_top50$Plot))


BSSzmat.dat=list()
BSScover.dat=list()

BSSothercover.dat=list()
BSSotherrich.dat=list()

###BSSzmat.dat 1/0 data of each of the 100 species across plots and field Ages
###BSScover.dat cover of each of the 100 species across plots and field Ages
###BSSothercover.dat sum of cover of other  species across plots and field Ages
###BSSotherrich.dat  richness of other  species across plots and field Ages
splist = names(freq_sps_top_50)
Ageset = unique(sp_cover_f2_top50$fake_age)
for (kk in 1:length(splist)){
  
  sp=subset(sp_cover_f2_top50,sp_cover_f2_top50$Species==splist[kk]) 
  BSSother=subset(sp_cover_f2_top50,sp_cover_f2_top50$Species!=splist[kk])
  
  exthead=rep(NA,length(Ageset))
  other.cover=rep(NA,length(Ageset))
  other.rich=rep(NA,length(Ageset))
  
  for (i in 1:length(Fieldname)){
    zmat=matrix(NA,nrow=48,ncol=length(Ageset))
    for (j in 1:length(Plot)){
      dat1=sp[which(sp$FIELDNAME==Fieldname[i]&sp$Plot==Plot[j]),]
      for (t in 1:length(Ageset)){
        dat2=subset(dat1,dat1$Age==Ageset[t])
        if (length(dat2$COVER)==0) zmat[j,t]=NA
        if (length(dat2$COVER)==1) zmat[j,t]=dat2$COVER 
      }
      
    }
    exthead=rbind(exthead,zmat)  
    
  }
  
  
  ###other species cover and richness
  
  for (i in 1:length(Fieldname)){
    zcovermat=matrix(NA,nrow=48, ncol = length(Ageset))
    zrichmat=matrix(NA,nrow = 48,ncol = length(Ageset))
    
    for (j in 1:length(Plot)){
      dat1=subset(BSSother, BSSother$FIELDNAME==Fieldname[i]&BSSother$Plot==Plot[j])
      for (t in 1:length(Ageset)){
        dat2=subset(dat1,dat1$Age==Ageset[t])
        
        if (length(dat2$COVER)==0) zcovermat[j,t]=NA
        if (length(dat2$COVER)>0) zcovermat[j,t]=sum(dat2$COVER)
        
        if (length(unique(dat2$Species))==0) zrichmat[j,t]=NA
        if (length(unique(dat2$Species))>0) zrichmat[j,t]=length(unique(dat2$Species))
      }
      
    }
    other.cover=rbind(other.cover,zcovermat)
    other.rich=rbind(other.rich,zrichmat)
  }
  
  
  exthead=exthead[-1,]
  colnames(exthead)=c(paste0("Age",Ageset))
  rownames(exthead)=as.character(seq(1:480))
  
  extheadcover=exthead
  extheadcover[is.na(extheadcover)==1]=0
  
  exthead[exthead>0]=1
  exthead[is.na(exthead)==1]=0
  
  
  other.cover=other.cover[-1,]
  other.rich=other.rich[-1,]
  colnames(other.cover)=c(paste0("Age",Ageset))
  rownames(other.cover)=as.character(seq(1:480))
  colnames(other.rich)=c(paste0("Age",Ageset))
  rownames(other.rich)=as.character(seq(1:480))
  
  other.cover[is.na(other.cover)==1]=0
  other.rich[is.na(other.rich)==1]=0
  
  
  BSSzmat.dat[[kk]]=exthead
  BSScover.dat[[kk]]=extheadcover
  BSSothercover.dat[[kk]]=other.cover
  BSSotherrich.dat[[kk]]=other.rich
}


###growth.list

cover.list=list()
for (i in 1:length(splist)){
  dat=BSScover.dat[[i]]
  dat[dat==0]=1e-6
  cover.list[[i]]=dat
}



growth.list=list()
for (i in 1:length(splist)){
  i = 5
  dat=cover.list[[i]]
  growth.dat=rep(NA,480)
  for (j in 2:dim(dat)[2]){
    #j = 2
    colnames(dat)
    t_diff = as.numeric(gsub('Age','',colnames(dat)[j]))-as.numeric(gsub('Age','',colnames(dat)[j-1]))
    growth=(log(dat[,j])-log(dat[,(j-1)]))/t_diff
    growth.dat=data.frame(growth.dat,growth)
    
  }
  growth.dat=growth.dat[,-1]
  colnames(growth.dat)=paste0("Age",Ageset[-42])
  growth.list[[i]]=growth.dat
}

fit_yin_ele = list(BSSzmat.dat, BSScover.dat, BSSothercover.dat, BSSotherrich.dat,
                   growth.list)
save(fit_yin_ele, file = 'D:/R projects/BSS/code/Yin_ele_2022/fit_yin_ele.rdata')
