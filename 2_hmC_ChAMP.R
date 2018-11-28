source("~/NavarraBiomed/analysis/parkinson/code/experiment.r")

dataDir <- "~/NavarraBiomed/analysis/parkinson/data/1_no_outliers"
workDir <- "~/NavarraBiomed/analysis/parkinson/results/2_hmC"
setwd(workDir)

myLoad.file <- "myLoad.Rdata"
if(!file.exists(myLoad.file)){
  myLoad <- champ.load(directory = dataDir, arraytype = "EPIC", method="minfi")
  myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
  save(myLoad, file=myLoad.file)
} else {
  load(myLoad.file)
}

pd <- myLoad$pd

champ.QC(pheno=myLoad$pd$Status)

myNorm.file <- "myNorm.Rdata"
if(!file.exists(myNorm.file)){
  myNorm <- champ.norm(beta=myLoad$beta,mset=myLoad$mset,rgSet=myLoad$rgSet,method="BMIQ",arraytype="EPIC")
  save(myNorm, file=myNorm.file)
} else {
  load(myNorm.file)
}

myDMP.file <- "myDMP.Rdata"
if(!file.exists(myDMP.file)){
  myDMP <- champ.DMP(beta = myNorm, pheno = pd$Status, arraytype = "EPIC", adjPVal=1)
  save(myDMP, file=myDMP.file)
} else {
  load(myDMP.file)
}


