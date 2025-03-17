###spatio-temporal data
###data sorting
library(spacetime)
library(sp)
library(foreign)

C3<-read.csv("C3.csv")
C4<-read.csv("C4.csv")
C5<-read.csv("C5.csv")
C6<-read.csv("C6.csv")
C7<-read.csv("C7.csv")
D1<-read.csv("D1.csv")
D2<-read.csv("D2.csv")
D3<-read.csv("D3.csv")
E1<-read.csv("E1.csv")
E2<-read.csv("E2.csv")




Ageset=seq(1,47,by=2)

Ageset%in%C3$AGE  
C3new=C3[which(C3$AGE%in%Ageset==1),]

Ageset%in%C4$AGE
if (sum(Ageset%in%C4$AGE)==length(Ageset)){
  C4new=C4[which(C4$AGE%in%Ageset==1),]
}
if (sum(Ageset%in%C4$AGE)!=length(Ageset)){
  AgeC4<-sort(unique(C4$AGE))
  
  Agesetnew<-c(intersect(Ageset,AgeC4),setdiff(Ageset,AgeC4)+1)
  
  C4new=C4[which(C4$AGE%in%Agesetnew==1),]

}


Ageset%in%C5$AGE
if (sum(Ageset%in%C5$AGE)==length(Ageset)){
  C5new=C5[which(C5$AGE%in%Ageset==1),]
}
if (sum(Ageset%in%C5$AGE)!=length(Ageset)){
  AgeC5<-sort(unique(C5$AGE))

  Agesetnew<-c(intersect(Ageset,AgeC5),setdiff(Ageset,AgeC5)+1)
  
  C5new=C5[which(C5$AGE%in%Agesetnew==1),]
}

C5new[which(C5new$AGE%in%Ageset==0),"AGE"]<-C5new[which(C5new$AGE%in%Ageset==0),"AGE"]-1



Ageset%in%C6$AGE
if (sum(Ageset%in%C6$AGE)==length(Ageset)){
  C6new=C6[which(C6$AGE%in%Ageset==1),]
}
if (sum(Ageset%in%C6$AGE)!=length(Ageset)){
  AgeC6<-sort(unique(C6$AGE))
  
  Agesetnew<-c(intersect(Ageset,AgeC6),setdiff(Ageset,AgeC6)+1)
  
  C6new=C6[which(C6$AGE%in%Agesetnew==1),]
}

C6new[which(C6new$AGE%in%Ageset==0),"AGE"]<-C6new[which(C6new$AGE%in%Ageset==0),"AGE"]-1
sum(C6new$AGE%in%Ageset)==dim(C6new)[1]

Ageset%in%C7$AGE
if (sum(Ageset%in%C7$AGE)==length(Ageset)){
  C7new=C7[which(C7$AGE%in%Ageset==1),]
}
if (sum(Ageset%in%C7$AGE)!=length(Ageset)){
  AgeC7<-sort(unique(C7$AGE))
  
  Agesetnew<-c(intersect(Ageset,AgeC7),setdiff(Ageset,AgeC7)+1)
  
  C7new=C7[which(C7$AGE%in%Agesetnew==1),]
}
C7new[which(C7new$AGE%in%Ageset==0),"AGE"]<-C7new[which(C7new$AGE%in%Ageset==0),"AGE"]-1
sum(C7new$AGE%in%Ageset)==dim(C7new)[1]


Ageset%in%D1$AGE
if (sum(Ageset%in%D1$AGE)==length(Ageset)){
  D1new=D1[which(D1$AGE%in%Ageset==1),]
}
if (sum(Ageset%in%D1$AGE)!=length(Ageset)){
  AgeD1<-sort(unique(D1$AGE))
  
  Agesetnew<-c(intersect(Ageset,AgeD1),setdiff(Ageset,AgeD1)+1)
  
  D1new=D1[which(D1$AGE%in%Agesetnew==1),]
}

D1new[which(D1new$AGE%in%Ageset==0),"AGE"]<-D1new[which(D1new$AGE%in%Ageset==0),"AGE"]-1
sum(D1new$AGE%in%Ageset)==dim(D1new)[1]


Ageset%in%D2$AGE
if (sum(Ageset%in%D2$AGE)==length(Ageset)){
  D2new=D2[which(D2$AGE%in%Ageset==1),]
}
if (sum(Ageset%in%D2$AGE)!=length(Ageset)){
  AgeD2<-sort(unique(D2$AGE))
  
  Agesetnew<-c(intersect(Ageset,AgeD2),setdiff(Ageset,AgeD2)+1)
  
  D2new=D2[which(D2$AGE%in%Agesetnew==1),]
}

D2new[which(D2new$AGE%in%Ageset==0),"AGE"]<-D2new[which(D2new$AGE%in%Ageset==0),"AGE"]-1
sum(D2new$AGE%in%Ageset)==dim(D2new)[1]


Ageset%in%D3$AGE
if (sum(Ageset%in%D3$AGE)==length(Ageset)){
  D3new=D3[which(D3$AGE%in%Ageset==1),]
}
if (sum(Ageset%in%D3$AGE)!=length(Ageset)){
  AgeD3<-sort(unique(D3$AGE))
  
  Agesetnew<-c(intersect(Ageset,AgeD3),setdiff(Ageset,AgeD3)+1)
  
  D3new=D3[which(D3$AGE%in%Agesetnew==1),]
}

D3new[which(D3new$AGE%in%Ageset==0),"AGE"]<-D3new[which(D3new$AGE%in%Ageset==0),"AGE"]-1
sum(D3new$AGE%in%Ageset)==dim(D3new)[1]


Ageset%in%E1$AGE
if (sum(Ageset%in%E1$AGE)==length(Ageset)){
  E1new=E1[which(E1$AGE%in%Ageset==1),]
}
if (sum(Ageset%in%E1$AGE)!=length(Ageset)){
  AgeE1<-sort(unique(E1$AGE))
  
  Agesetnew<-c(intersect(Ageset,AgeE1),setdiff(Ageset,AgeE1)+1)
  
  E1new=E1[which(E1$AGE%in%Agesetnew==1),]
}

E1new[which(E1new$AGE%in%Ageset==0),"AGE"]<-E1new[which(E1new$AGE%in%Ageset==0),"AGE"]-1
sum(E1new$AGE%in%Ageset)==dim(E1new)[1]

Ageset%in%E2$AGE
if (sum(Ageset%in%E2$AGE)==length(Ageset)){
  E2new=E2[which(E2$AGE%in%Ageset==1),]
}
if (sum(Ageset%in%E2$AGE)!=length(Ageset)){
  AgeE2<-sort(unique(E2$AGE))
  
  Agesetnew<-c(intersect(Ageset,AgeE2),setdiff(Ageset,AgeE2)+1)
  
  E2new=E2[which(E2$AGE%in%Agesetnew==1),]
}
E2new[which(E2new$AGE%in%Ageset==0),"AGE"]<-E2new[which(E2new$AGE%in%Ageset==0),"AGE"]-1
sum(E2new$AGE%in%Ageset)==dim(E2new)[1]

BSS=rbind(C3new,C4new,C5new,C6new,C7new,D1new,D2new,D3new,E1new,E2new)
sp<-paste(BSS$GENUS,BSS$EPITHET,sep="_")
BSS<-cbind(BSS,sp)

###exclude woody species
BSS=subset(BSS,BSS$GENERAL.FORM=="Annual"|BSS$GENERAL.FORM=="Perennial"|BSS$GENERAL.FORM=="Biennial")
BSsp=unique(BSS$sp)  ###273 species, 97486 lines


frquen<-NULL
for (i in 1:length(BSsp)){
 dat<-BSS[which(BSS$sp==BSsp[i]),]
 frquen[i]=dim(dat)[1]
}

splist<-BSsp[order(frquen,decreasing=T)[1:100]]
frq=frquen[order(frquen,decreasing=T)][1:100]

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



###BSSzmat.dat 1/0 data of each of the 100 species across plots and field ages
###BSScover.dat cover of each of the 100 species across plots and field ages
###BSSothercover.dat sum of cover of other  species across plots and field ages
###BSSotherrich.dat  richness of other  species across plots and field ages

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
        dat2=subset(dat1,dat1$AGE==Ageset[t])
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
        dat2=subset(dat1,dat1$AGE==Ageset[t])
        
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


