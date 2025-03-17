###spatio-temporal data
###data sorting
library(spacetime)
library(sp)
library(foreign)

##### Data loading
sp_cover = read.table('D:/R projects/BSS/data/Li_data/samplewhole.txt',
                      header = T)
field = read.csv('D:/R projects/BSS/data/original data/FIELDS.csv',
                   header = T)
trait = read.csv('D:/R projects/BSS/data/original data/traits332.csv',
                 header = T)

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




Ageset=seq(1,47,by=2)

Ageset%in%C3$Age  
C3new=C3[which(C3$Age%in%Ageset==1),]

Ageset%in%C4$Age
if (sum(Ageset%in%C4$Age)==length(Ageset)){
  C4new=C4[which(C4$Age%in%Ageset==1),]
}
if (sum(Ageset%in%C4$Age)!=length(Ageset)){
  AgeC4<-sort(unique(C4$Age))
  
  Agesetnew<-c(intersect(Ageset,AgeC4),setdiff(Ageset,AgeC4)+1)
  
  C4new=C4[which(C4$Age%in%Agesetnew==1),]

}


Ageset%in%C5$Age
if (sum(Ageset%in%C5$Age)==length(Ageset)){
  C5new=C5[which(C5$Age%in%Ageset==1),]
}
if (sum(Ageset%in%C5$Age)!=length(Ageset)){
  AgeC5<-sort(unique(C5$Age))

  Agesetnew<-c(intersect(Ageset,AgeC5),setdiff(Ageset,AgeC5)+1)
  
  C5new=C5[which(C5$Age%in%Agesetnew==1),]
}

C5new[which(C5new$Age%in%Ageset==0),"Age"]<-C5new[which(C5new$Age%in%Ageset==0),"Age"]-1



Ageset%in%C6$Age
if (sum(Ageset%in%C6$Age)==length(Ageset)){
  C6new=C6[which(C6$Age%in%Ageset==1),]
}
if (sum(Ageset%in%C6$Age)!=length(Ageset)){
  AgeC6<-sort(unique(C6$Age))
  
  Agesetnew<-c(intersect(Ageset,AgeC6),setdiff(Ageset,AgeC6)+1)
  
  C6new=C6[which(C6$Age%in%Agesetnew==1),]
}

C6new[which(C6new$Age%in%Ageset==0),"Age"]<-C6new[which(C6new$Age%in%Ageset==0),"Age"]-1
sum(C6new$Age%in%Ageset)==dim(C6new)[1]

Ageset%in%C7$Age
if (sum(Ageset%in%C7$Age)==length(Ageset)){
  C7new=C7[which(C7$Age%in%Ageset==1),]
}
if (sum(Ageset%in%C7$Age)!=length(Ageset)){
  AgeC7<-sort(unique(C7$Age))
  
  Agesetnew<-c(intersect(Ageset,AgeC7),setdiff(Ageset,AgeC7)+1)
  
  C7new=C7[which(C7$Age%in%Agesetnew==1),]
}
C7new[which(C7new$Age%in%Ageset==0),"Age"]<-C7new[which(C7new$Age%in%Ageset==0),"Age"]-1
sum(C7new$Age%in%Ageset)==dim(C7new)[1]


Ageset%in%D1$Age
if (sum(Ageset%in%D1$Age)==length(Ageset)){
  D1new=D1[which(D1$Age%in%Ageset==1),]
}
if (sum(Ageset%in%D1$Age)!=length(Ageset)){
  AgeD1<-sort(unique(D1$Age))
  
  Agesetnew<-c(intersect(Ageset,AgeD1),setdiff(Ageset,AgeD1)+1)
  
  D1new=D1[which(D1$Age%in%Agesetnew==1),]
}

D1new[which(D1new$Age%in%Ageset==0),"Age"]<-D1new[which(D1new$Age%in%Ageset==0),"Age"]-1
sum(D1new$Age%in%Ageset)==dim(D1new)[1]


Ageset%in%D2$Age
if (sum(Ageset%in%D2$Age)==length(Ageset)){
  D2new=D2[which(D2$Age%in%Ageset==1),]
}
if (sum(Ageset%in%D2$Age)!=length(Ageset)){
  AgeD2<-sort(unique(D2$Age))
  
  Agesetnew<-c(intersect(Ageset,AgeD2),setdiff(Ageset,AgeD2)+1)
  
  D2new=D2[which(D2$Age%in%Agesetnew==1),]
}

D2new[which(D2new$Age%in%Ageset==0),"Age"]<-D2new[which(D2new$Age%in%Ageset==0),"Age"]-1
sum(D2new$Age%in%Ageset)==dim(D2new)[1]


Ageset%in%D3$Age
if (sum(Ageset%in%D3$Age)==length(Ageset)){
  D3new=D3[which(D3$Age%in%Ageset==1),]
}
if (sum(Ageset%in%D3$Age)!=length(Ageset)){
  AgeD3<-sort(unique(D3$Age))
  
  Agesetnew<-c(intersect(Ageset,AgeD3),setdiff(Ageset,AgeD3)+1)
  
  D3new=D3[which(D3$Age%in%Agesetnew==1),]
}

D3new[which(D3new$Age%in%Ageset==0),"Age"]<-D3new[which(D3new$Age%in%Ageset==0),"Age"]-1
sum(D3new$Age%in%Ageset)==dim(D3new)[1]


Ageset%in%E1$Age
if (sum(Ageset%in%E1$Age)==length(Ageset)){
  E1new=E1[which(E1$Age%in%Ageset==1),]
}
if (sum(Ageset%in%E1$Age)!=length(Ageset)){
  AgeE1<-sort(unique(E1$Age))
  
  Agesetnew<-c(intersect(Ageset,AgeE1),setdiff(Ageset,AgeE1)+1)
  
  E1new=E1[which(E1$Age%in%Agesetnew==1),]
}

E1new[which(E1new$Age%in%Ageset==0),"Age"]<-E1new[which(E1new$Age%in%Ageset==0),"Age"]-1
sum(E1new$Age%in%Ageset)==dim(E1new)[1]

Ageset%in%E2$Age
if (sum(Ageset%in%E2$Age)==length(Ageset)){
  E2new=E2[which(E2$Age%in%Ageset==1),]
}
if (sum(Ageset%in%E2$Age)!=length(Ageset)){
  AgeE2<-sort(unique(E2$Age))
  
  Agesetnew<-c(intersect(Ageset,AgeE2),setdiff(Ageset,AgeE2)+1)
  
  E2new=E2[which(E2$Age%in%Agesetnew==1),]
}
E2new[which(E2new$Age%in%Ageset==0),"Age"]<-E2new[which(E2new$Age%in%Ageset==0),"Age"]-1
sum(E2new$Age%in%Ageset)==dim(E2new)[1]

BSS=rbind(C3new,C4new,C5new,C6new,C7new,D1new,D2new,D3new,E1new,E2new)
sp<-paste(BSS$GENUS,BSS$EPITHET,sep="_")
BSS<-cbind(BSS,sp)

###exclude woody species
#BSS=subset(BSS,BSS$GENERAL.FORM=="Annual"|BSS$GENERAL.FORM=="Perennial"|BSS$GENERAL.FORM=="Biennial")
BSS_1 = sp_cover %>% select(intersect(trait[trait$Span %in% c("Annual", "Perennial", "Biennial"),]$Species, colnames(sp_cover)[7:ncol(sp_cover)]))
BSsp=unique(BSS$sp)  ###273 species, 97486 lines

sp = colnames(BSS)[7:ncol(BSS)]
frequen_1<-data.frame(sp = NA, freq = NA)
frequen<-data.frame(sp = NA, freq = NA)
for (i in 1:length(sp)){
 dat<-BSS %>% select(sp[i])
 frequen_1$freq=length(dat[dat>0])
 frequen_1$sp=sp[i]
 frequen = rbind(frequen_1, frequen)
}

splist <- arrange(frequen, -freq)$sp[1:100]

###data of 100 most frequent species 

dt=NULL
for (i in 1:length(splist)){
  dat=subset(BSS,BSS$sp==splist[i])
  dt=rbind(dt,dat)
}

dim(dt)[1] ###94923 lines

BSS100=dt

Origin=NULL
for (i in 1:length(splist)){
  origin.dt=subset(BSS,BSS$sp==splist[i])
  if (length(unique(origin.dt$ORIGIN))==1) Origin[i]=unique(origin.dt$ORIGIN)
  
}

sp.inf=data.frame(splist,frq,Origin)

BSS100.native=subset(BSS100,BSS100$ORIGIN=="N")
length(unique(BSS100.native$sp))  ##58 native sp
BSS100.exotic=subset(BSS100,BSS100$ORIGIN=="E")
length(unique(BSS100.exotic$sp))  ###40 extic sp ##2 unknown


BSS100.annual=subset(BSS100,BSS100$GENERAL.FORM=="Annual")
length(unique(BSS100.annual$sp))  ##20 annual sp
BSS100.biennial=subset(BSS100,BSS100$GENERAL.FORM=="Biennial")
length(unique(BSS100.biennial$sp))  ###15 biennial sp 
BSS100.perennial=subset(BSS100,BSS100$GENERAL.FORM=="Perennial")
length(unique(BSS100.perennial$sp))  ###65 perennial sp 

Field<-c(rep("C3",48),rep("C4",48),rep("C5",48),rep("C6",48),rep("C7",48),rep("D1",48),rep("D2",48),rep("D3",48),rep("E1",48),rep("E2",48))
Plotid<-rep(seq(1:48),10)
Fieldnum<-c(rep(1,48),rep(2,48),rep(3,48),rep(4,48),rep(5,48),rep(6,48),rep(7,48),rep(8,48),rep(9,48),rep(10,48))

###creat 0-1 data of each species
Fieldname<-unique(BSSexotic$FIELDNAME)
Plot<-sort(unique(BSSexotic$PLOTID))


BSSzmat.dat<-list()
BSScover.dat<-list()

BSSothercover.dat<-list()
BSSotherrich.dat<-list()



###BSSzmat.dat 1/0 data of each of the 100 species across plots and field Ages
###BSScover.dat cover of each of the 100 species across plots and field Ages
###BSSothercover.dat sum of cover of other  species across plots and field Ages
###BSSotherrich.dat  richness of other  species across plots and field Ages

for (kk in 1:length(splist)){
  
  sp<-subset(BSS100,BSS100$sp==splist[kk]) 
  BSSother<-subset(BSS100,BSS100$sp!=splist[kk])
  
  exthead=rep(NA,length(Ageset))
  other.cover=rep(NA,length(Ageset))
  other.rich=rep(NA,length(Ageset))
  
  for (i in 1:length(Fieldname)){
    zmat=matrix(NA,nrow=48,ncol=24)
    for (j in 1:length(Plot)){
      dat1=sp[which(sp$FIELDNAME==Fieldname[i]&sp$PLOTID==Plot[j]),]
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
    zcovermat=matrix(NA,nrow=48,ncol=24)
    zrichmat=matrix(NA,nrow = 48,ncol = 24)
    
    for (j in 1:length(Plot)){
      dat1=subset(BSSother,BSSother$FIELDNAME==Fieldname[i]&BSSother$PLOTID==Plot[j])
      for (t in 1:length(Ageset)){
        dat2=subset(dat1,dat1$Age==Ageset[t])
        
        if (length(dat2$COVER)==0) zcovermat[j,t]=NA
        if (length(dat2$COVER)>0) zcovermat[j,t]=sum(dat2$COVER)
        
        if (length(unique(dat2$sp))==0) zrichmat[j,t]=NA
        if (length(unique(dat2$sp))>0) zrichmat[j,t]=length(unique(dat2$sp))
      }
      
    }
    other.cover=rbind(other.cover,zcovermat)
    other.rich=rbind(other.rich,zrichmat)
  }
  
  
  exthead=exthead[-1,]
  colnames(exthead)<-c(paste0("Age",1:24))
  rownames(exthead)<-as.character(seq(1:480))
  
  extheadcover=exthead
  extheadcover[is.na(extheadcover)==1]<-0
  
  exthead[exthead>0]<-1
  exthead[is.na(exthead)==1]<-0
  
  
  other.cover=other.cover[-1,]
  other.rich=other.rich[-1,]
  colnames(other.cover)<-c(paste0("Age",1:24))
  rownames(other.cover)<-as.character(seq(1:480))
  colnames(other.rich)<-c(paste0("Age",1:24))
  rownames(other.rich)<-as.character(seq(1:480))
  
  other.cover[is.na(other.cover)==1]<-0
  other.rich[is.na(other.rich)==1]<-0
  
  
  BSSzmat.dat[[kk]]<-exthead
  BSScover.dat[[kk]]<-extheadcover
  BSSothercover.dat[[kk]]<-other.cover
  BSSotherrich.dat[[kk]]<-other.rich
}


###growth.list

cover.list=list()
for (i in 1:length(splist)){
  dat=BSScover.dat[[i]]
  dat[dat==0]<-1e-6
  cover.list[[i]]=dat
}

growth.list=list()
for (i in 1:length(splist)){
  dat=cover.list[[i]]
  growth.dat=rep(NA,480)
  for (j in 2:dim(dat)[2]){
    growth=(log(dat[,j])-log(dat[,(j-1)]))/2
    growth.dat=data.frame(growth.dat,growth)
    
  }
  growth.dat=growth.dat[,-1]
  colnames(growth.dat)<-paste0("Age",1:23)
  growth.list[[i]]=growth.dat
}


