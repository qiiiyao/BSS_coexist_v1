###pairwise interaction effects
###BSSzmat.dat 1/0 data of 25 native and 25 exotic species
###BSScover.dat cover of each species

###probablity of occurrence~cover of other invader

`pairoccurrence.fn` <-
  function(){
    # this script fits the modelto the BSS pairwise biotic interactions
    
    library(R2jags)
    
 sim.mod<-function(){
        
        psi~dunif(0,10)
        b ~ dnorm(0,.1) 
        
        
        sigmaplot~dunif(0,.1)
        tauplot<-sqrt(1/(sigmaplot*sigmaplot))
        
        for(i in 1:10){
        c[i] ~ dnorm(0,tauplot)
        
        }
        
        for (i in 1:(nyear-1)){
          a[i] ~ dnorm(0,.1)
        }
        
        for(i in 1:nsite){
        
          for(t in 2:nyear){
        logit(muZ[i,t])<-a[t-1]+b*extcover[i,t-1]+c[Fieldnum[i]]
        
        z[i,t]~dbern(muZ[i,t])
        }
        }
        
        }

    
    
    z=z
    nsite<-dim(z)[1]
    
    nyear<-dim(z)[2]
    
   
    
    extcover=extcover
    
    Fieldnum<-Fieldnum
    
    sim.dat <- list ()
    sim.dat$z<-z
    sim.dat$nsite<-nsite
    sim.dat$nyear<-nyear
    sim.dat$extcover<-extcover
    sim.dat$Fieldnum<-Fieldnum
    
   sim.par<-c("a","b","c") 
    
  sim.fit<-jags(sim.dat,model.file=sim.mod,parameters.to.save=sim.par,n.chains=1,n.iter=1e4)

    return(sim.fit)
    
    
    
  }

pairwiseoccurrence.res<-matrix(NA,50,50)
occurrence.sd<-matrix(NA,50,50)
occurrence.25<-matrix(NA,50,50)
occurrence.975<-matrix(NA,50,50)

for (i in 1:50){
  
  z<-BSSzmat.dat[[i]]
  for (j in 1:50){
    
    extcover<-as.matrix(scale(BSScover.dat[[j]],center=TRUE))
    extcover[is.na(extcover)]=0
    res=pairoccurrence.fn()
    pairwiseoccurrence.res[i,j]= res$BUGSoutput$summary[2,1]
    occurrence.sd[i,j]<- res$BUGSoutput$summary[2,2]
    occurrence.25[i,j]<- res$BUGSoutput$summary[2,3]
    occurrence.975[i,j]<- res$BUGSoutput$summary[2,7]
  }
  
}

pairoccurrence.list=list(pairwiseoccurrence.res,occurrence.sd,occurrence.25,occurrence.975)





###colonization and survival of invading species 
###function of cover of other invading speices
###pairwise colonization and survival effect


`paircolsurv.fn` <-
  function(){
    # this script fits the model  to the survival and colonization of BSS invasive species
    
    library(R2jags)
    
    sim.mod<-function(){
        
        psi~dunif(0,1)
    
        b2 ~ dnorm(0,.1) 
        b3 ~ dnorm(0,.1)
        
        
        for(i in 1:(nyear-1)){
          a [i] ~ dnorm(0,.1)
        b1[i] ~ dnorm(0,.1) 
        
        }
        
        sigmaplot~dunif(0,10)
        tauplot<-sqrt(1/(sigmaplot*sigmaplot))
        
        for(i in 1:10){
        c[i] ~ dnorm(0,tauplot)
        
        
        }
        
        
        for(i in 1:nsite){
        z[i,1]~dbern(psi)
        for(t in 2:nyear){
        logit(muZ[i,t])<- a[t-1]+b1[t-1]*z[i,t-1]+b2*extcover[i,t-1]*z[i,t-1]+b3*extcover[i,t-1]*(1-z[i,t-1])+c[Fieldnum[i]]
        
        z[i,t]~dbern(muZ[i,t])
        }
        }
        
        
        }

    

    z=z
    nsite<-dim(z)[1]
    
    nyear<-dim(z)[2]
    
    
    
    extcover=extcover
    
    Fieldnum<-Fieldnum
    
    sim.dat <- list ()
    sim.dat$z<-z
    sim.dat$nsite<-nsite
    sim.dat$nyear<-nyear
    sim.dat$extcover<-extcover
    sim.dat$Fieldnum<-Fieldnum
    
    sim.par<-c("a","b1","b2","b3","c") 
    
    sim.fit<-jags(sim.dat,model.file=sim.mod,parameters.to.save=sim.par,n.chains=1,n.iter=1e4)
    
    return(sim.fit)
    

    
  }

pairsurv.res<-matrix(NA,50,50)
pairsurv.sd<-matrix(NA,50,50)
pairsurv.25<-matrix(NA,50,50)
pairsurv.975<-matrix(NA,50,50)

paircol.res<-matrix(NA,50,50)
paircol.sd<-matrix(NA,50,50)
paircol.25<-matrix(NA,50,50)
paircol.975<-matrix(NA,50,50)



for (i in 1:50){
  
  z<-BSSzmat.dat[[i]]
  for (j in 1:50){
    
    extcover<-as.matrix(scale(BSScover.dat[[j]],center=TRUE))
    extcover[is.na(extcover)]=0
    res=paircolsurv.fn()
    pairsurv.res[i,j]= res$BUGSoutput$summary[47,1]
    pairsurv.sd[i,j]=res$BUGSoutput$summary[47,2]
    pairsurv.25[i,j]=res$BUGSoutput$summary[47,3]
    pairsurv.975[i,j]=res$BUGSoutput$summary[47,7]
    
    paircol.res[i,j]<-res$BUGSoutput$summary[48,1]
    paircol.sd[i,j]<- res$BUGSoutput$summary[48,2]
    paircol.25[i,j]<- res$BUGSoutput$summary[48,3]
    paircol.975[i,j]<- res$BUGSoutput$summary[48,7]
  }
  
}

paircolsurv.list=list(pairsurv.res,pairsurv.sd,pairsurv.25,pairsurv.975,paircol.res,paircol.sd,paircol.25,paircol.975)





###growth rate 
`pairgrowth.fn` <-
  function(){
    # this script fits the model  to the BSS invading species by evaluate the invader-invader effects 
    
    library(R2jags)
    
    sim.mod<-function(){
      
     for (i in 1: nyear){
       a[i] ~ dnorm(0,.001) 
     }
 
      b ~ dnorm(0,.001) 
      tau ~ dgamma(1.0e-3,1.0e-3)
      
      sigmaplot~dunif(0,.1)
      tauplot<-sqrt(1/(sigmaplot*sigmaplot))
      
      for(i in 1:10){
        c[i] ~ dnorm(0,tauplot)
        
      }
      
      for(i in 1:nsite){
        
        for(t in 1:nyear){
          growth.dat[i,t]~dnorm(mu[i,t],tau)
          
          mu[i,t]<-a[t]+b*extcover[i,t]+c[Fieldnum[i]]
        }
      }
      
    }
    
    
    
    growth.dat=growth.dat
    nsite<-dim(growth.dat)[1]
    
    nyear<-dim(growth.dat)[2]
    
    
    
    extcover=extcover
    
    Fieldnum<-Fieldnum
    
    sim.dat <- list ()
    sim.dat$growth.dat<-growth.dat
    sim.dat$nsite<-nsite
    sim.dat$nyear<-nyear
    sim.dat$extcover<-extcover
    sim.dat$Fieldnum<-Fieldnum
    
    sim.par<-c("a","b","c") 
    
    sim.fit<-jags(sim.dat,model.file=sim.mod,parameters.to.save=sim.par,n.chains=1,n.iter=1e4)
    
    return(sim.fit)
  }



pairgrowth.res<-matrix(NA,50,50)
growth.sd<-matrix(NA,50,50)
growth.25<-matrix(NA,50,50)
growth.975<-matrix(NA,50,50)

for (i in 1:50){
  
  growth.dat<-growth.list[[i]]
  for (j in 1:50){
    
    extcover<-as.matrix(scale(BSScover.dat[[j]],center=TRUE))
    extcover[is.na(extcover)]=0
    res=pairgrowth.fn()
    pairgrowth.res[i,j]= res$BUGSoutput$summary[24,1]
    growth.sd[i,j]<- res$BUGSoutput$summary[24,2]
    growth.25[i,j]<- res$BUGSoutput$summary[24,3]
    growth.975[i,j]<- res$BUGSoutput$summary[24,7]
  }
  
}

pairgrowth.list=list(pairgrowth.res,growth.sd,growth.25,growth.975)





