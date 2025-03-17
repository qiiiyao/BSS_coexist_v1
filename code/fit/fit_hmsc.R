library(Hmsc)
library(R.matlab)

### Sample for reduced-rank regression
load("D:/R projects/BSS/code/data preparation/transformed data/fit_fp_same_ages_top50.RData")

G = fit_fp_same_ages_top50[[1]][[3]][[1]]$G
X = fit_fp_same_ages_top50[[1]][[3]][[1]]$X
T = fit_fp_same_ages_top50[[1]][[3]][[1]]$time
write.csv(cbind(G, T, X), 'D:/Paper_collection/Ovaskainen_et_al_2017_PORS.B/61328-hmsc_matlab_2.1/BSS/data/data.csv')

m = Hmsc(Y = G, XData = X, XFormula = ~1,
         XRRRData = X, XRRRFormula = ~.-1, ncRRR=1,
         distr = "normal")
m = sampleMcmc(m, thin = 10, samples = 500,
               transient = 0, nChains = 4,
               verbose = T)
m$X
m$XInterceptInd



#set working directory to the project folder
#setwd("")

localDir = "."
data.directory = file.path(localDir, "data")
model.directory = file.path(localDir, "models")

abutype = 2 #1: presence-absence; 2: DNA-abundance
climtype = 2 #1: summer; 2: winter

load(file = file.path(data.directory,"hmscData.RData"))

prev = colSums(Y)
sel = which(prev>9)
Y = Y[,sel]
TrData = data.frame(group = TrData$group[sel])
XDatasel = list()
XDatasel.abu = list()
for(i in 1:length(sel)){
  XDatasel[[i]] = XData[[sel[[i]]]]
  XDatasel.abu[[i]] = XData.abu[[sel[[i]]]]
}
if(abutype==1){
  XData = XDatasel
  Xz = Xz[,c(1,2,3,4,4+sel)]
} 
if(abutype==2){
  XData = XDatasel.abu
  Xz = Xzabu2[,c(1,2,3,4,4+sel)]
}

XFormula.1s = ~ temperature + clim.summer + poly(week,degree = 2,raw = TRUE)
XFormula.2s = ~ temperature + clim.summer + poly(week,degree = 2,raw = TRUE) + within.species
XFormula.3s = ~ temperature + clim.summer + poly(week,degree = 2,raw = TRUE) + within.species + Herbivore + Parasitoid + Predator + Saprophage_Detritivore
XFormula.4s = as.formula(paste("~ temperature + clim.summer + poly(week,degree = 2,raw = TRUE) + ",
                               paste(colnames(Y), collapse= " + ")))

XFormula.1w = ~ temperature + clim.winter + poly(week,degree = 2,raw = TRUE)
XFormula.2w = ~ temperature + clim.winter + poly(week,degree = 2,raw = TRUE) + within.species
XFormula.3w = ~ temperature + clim.winter + poly(week,degree = 2,raw = TRUE) + within.species + Herbivore + Parasitoid + Predator + Saprophage_Detritivore
XFormula.4w = as.formula(paste("~ temperature + clim.winter + poly(week,degree = 2,raw = TRUE) + ",
                               paste(colnames(Y), collapse= " + ")))

if(climtype==1){
  XFormulas = c(XFormula.1s, XFormula.2s, XFormula.3s, XFormula.4s)
} else {
  XFormulas = c(XFormula.1w, XFormula.2w, XFormula.3w, XFormula.4w)
}

ns = dim(Y)[2]
XSelect = list()
XSelect2 = list()
XSelect3 = list()
XSelect4 = list()
qq = 0.01
qq2 = 0.1
qq3 = 0.5
qq4 = 1.0
for (k in 1:ns){
  covGroup = k+4
  spGroup = 1:ns
  q = rep(qq, max(spGroup))
  q2 = rep(qq2, max(spGroup))
  q3 = rep(qq3, max(spGroup))
  q4 = rep(qq4, max(spGroup))
  XSelect[[k]] =
    list(covGroup = covGroup,
         spGroup = spGroup,q = q)
  XSelect2[[k]] =
    list(covGroup = covGroup,
         spGroup = spGroup,q = q2)
  XSelect3[[k]] =
    list(covGroup = covGroup,
         spGroup = spGroup,q = q3)
  XSelect4[[k]] =
    list(covGroup = covGroup,
         spGroup = spGroup,q = q4)
}

TrFormula = ~ group

rL.trap = HmscRandomLevel(units = levels(studyDesign$trap))
rL.year = HmscRandomLevel(units = levels(studyDesign$year))

nChains = 30
samples = 1000
for (thin in c(1,10,100)){
  transient = round(0.5*samples*thin)
  
  for (i in 7){
    if(i<4){
      m = Hmsc(Y=1*(Y>0), XData = XData, XFormula = XFormulas[[i]],
               TrData = TrData, TrFormula = TrFormula,
               studyDesign = studyDesign, ranLevels=list(trap=rL.trap,year=rL.year),
               distr="probit")
    } else {
      XSel = XSelect
      if(i==5) XSel = XSelect2
      if(i==6) XSel = XSelect3
      if(i==7) XSel = XSelect4
      m = Hmsc(Y=1*(Y>0), XData = Xz, XFormula = XFormulas[[4]],
               XSelect = XSel,
               TrData = TrData, TrFormula = TrFormula,
               studyDesign = studyDesign, ranLevels=list(trap=rL.trap,year=rL.year),
               distr="probit")
    }
    m = sampleMcmc(m, thin = thin, samples = samples, adaptNf = rep(round(0.4*thin*samples),2), 
                   transient = transient, nChains = nChains, nParallel = nChains)
    filename=file.path(model.directory, paste0("model_",as.character(i),
                                               "_climtype_",as.character(climtype),
                                               "_abutype_",as.character(abutype),
                                               "_chains_",as.character(nChains),
                                               "_samples_",as.character(samples),
                                               "_thin_",as.character(thin)))
    save(m,file=filename)
  }
}
