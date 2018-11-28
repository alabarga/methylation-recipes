source("~/NavarraBiomed/analysis/parkinson/code/experiment.r")

dataDir <- "~/NavarraBiomed/analysis/parkinson/data/2_BS_oxBS"
workDir <- "~/NavarraBiomed/analysis/parkinson/results/2_BS_oxBS"
infoDir <- "~/NavarraBiomed/analysis/parkinson/sampleInfo/"
setwd(workDir)

bsDir <- "~/NavarraBiomed/analysis/parkinson/data/2_BS"
oxDir <- "~/NavarraBiomed/analysis/parkinson/data/2_oxBS"

BS.file <- "BS.Rdata"
if(!file.exists(BS.file)){

  targetsBS <- read.metharray.sheet(bsDir, pattern="BS.csv$")  
  bsLoad <- champ.load(directory = bsDir, arraytype = "EPIC")
  save(bsLoad, targetsBS, file=BS.file)
  

  # rownames(targetsBS) <- targetsBS$Index  
  # BS.RG <- read.metharray.exp(base = dataDir, targets = targetsBS)
  # BS.RG <- read.metharray.exp(targets = targetsBS, extended=T)
  # save(BS.RG, targetsBS, file=BS.file)
  
  
} else {
  load(BS.file)
}

oxBS.file <- "oxBS.Rdata"
if(!file.exists(oxBS.file)){

  targetsOX <- read.metharray.sheet(oxDir, pattern="OX.csv$")
  oxLoad <- champ.load(directory = oxDir, arraytype = "EPIC")
  save(oxLoad, targetsOX, file=oxBS.file)

  
  # rownames(targetsOX) <- targetsOX$Index
  # #OX.RG <- read.metharray.exp(base = dataDir, targets = targetsOX)
  # OX.RG <- read.metharray.exp(targets = targetsOX, extended=T)
  # save(OX.RG, targetsOX, file=oxBS.file)
  
  
} else {
  load(oxBS.file)
}



bs1 <- preprocessFunnorm(BS.RG)
ox <- preprocessFunnorm(OX.RG)
bs2 <- champ.norm(beta=BS.RG$beta,mset=BS.RG$mset,rgSet=BS.RG$rgSet,method="FunctionalNormalization",arraytype="EPIC")


save(bs,ox,file="myNorm.Rdata")
load("norm.Rdata")

### load file

myNorm.file <- "myNorm.Rdata"

if(!file.exists(myNorm.file)){
  
  bs <- champ.norm(beta=bsLoad$beta,mset=bsLoad$mset,rgSet=bsLoad$rgSet,method="BMIQ",arraytype="EPIC")
  ox <- champ.norm(beta=oxLoad$beta,mset=oxLoad$mset,rgSet=v$rgSet,method="BMIQ",arraytype="EPIC")
  
  cpgs <- rownames(bs)[rownames(bs) %in% rownames(ox)]
  
  bs <- bs[cpgs,]
  ox <- ox[cpgs,]
  
} else {
  load(myNorm.file)
}


### BS vs OX
myDMP <- champ.DMP(beta = myNorm, pheno = pd$Sample_Group, arraytype = "EPIC", adjPVal=0.001)
hmc <- myDMP[[1]][abs(myDMP[[1]]$deltaBeta)>0,]
sel_dmp <- rownames(hmc)
p <- innerheatmap(myNorm[sel_dmp, order(pd$Sample_Group)])
p

beta.BS <- myNorm[sel_dmp, which(pd$Sample_Group == 'BS')]
beta.oxBS <- myNorm[sel_dmp, which(pd$Sample_Group == 'OX')]

M.BS <- logit2(beta.BS)
M.oxBS <- logit2(beta.oxBS)

## Naive estimates

## The naive approach to obtain 5-hmC levels is $\beta_{BS} - \beta_{OxBS}$. This approach results in negative values for the 5-hmC levels.

## hmC_naive <- beta.BS - beta.oxBS


## Using Ms
M_naive <- M.BS - M.oxBS
hmC_naive <- 2^M_naive / (2^M_naive + 1)


