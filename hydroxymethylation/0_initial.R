source("~/NavarraBiomed/analysis/parkinson/code/experiment.r")

dataDir <- "~/NavarraBiomed/analysis/parkinson/data/0_initial"
workDir <- "~/NavarraBiomed/analysis/parkinson/results/0_initial"
setwd(workDir)

myLoad.file <- "myLoad.Rdata"
if(!file.exists(myLoad.file)){
  myLoad <- champ.load(directory = dataDir, arraytype = "EPIC", method="minfi")
  save(myLoad, file=myLoad.file)
} else {
  load(myLoad.file)
}

CpG.GUI(CpG=rownames(myLoad$beta),arraytype="EPIC")
champ.QC()
for (samplename in colnames(myLoad$beta)){
  plotBetasByType(myLoad$beta[,samplename], annotationEPIC, main=samplename)
}

myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
pd <- myLoad$pd
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Sample_Group))), main="MDS (BSvsOX)")
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Slide))), main="MDS (Slide)")
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Status))), main="MDS (Disease)")

png('original_data_MDSPlot.png')
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Sample_Group))))
dev.off()

remove.outlier<-function(myLoad, samples.to.remove=c()){
  myLoad$pd <- myLoad$pd[!(myLoad$pd$Sample_Name %in% samples.to.remove),]
  myLoad$beta <- myLoad$beta[, !(colnames(myLoad$beta) %in% samples.to.remove)]
  myLoad$intensity <- myLoad$intensity[, !(colnames(myLoad$intensity) %in% samples.to.remove)]
}

samples.to.remove <- c("BK1_K13_BS")
myLoad<-remove.outlier(myLoad,samples.to.remove)
pd <- myLoad$pd

norm.method="SWAN"
myNorm.file <- paste("myNorm",norm.method,"Rdata", sep=".")
if(!file.exists(myNorm.file)){
  myNorm <- champ.norm(arraytype="EPIC", method=norm.method)
  save(myNorm, file=myNorm.file)
} else {
  load(myNorm.file)
}
champ.QC(beta=myNorm,pheno = pd$Sample_Group, PDFplot=FALSE)

png('normalized_data_MDSPlot.png')
par(mfrow=c(2,2))
plotMDS(myNorm, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Sample_Group))), main="MDS (BSvsOX)")
plotMDS(myNorm, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Slide))), main="MDS (Slide)")
plotMDS(myNorm, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Status))), main="MDS (Disease)")
plotMDS(myNorm, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Array))), main="MDS (Array)")
dev.off()

p <- innerheatmap(myNorm[, order(myLoad$pd$Sample_Group)])
p

champ.SVD(beta=myNorm)

myNorm.combat.file <- "myNorm.combat.Rdata"
if(!file.exists(myNorm.combat.file)){
  mod.combat = model.matrix( ~ 1 + myLoad$pd$Sample_Group + myLoad$pd$Status + myLoad$pd$Sex + myLoad$pd$Age)
  myNorm.cb <- sva::ComBat(dat = as.matrix(myNorm) , batch = myLoad$pd$Slide, mod=mod.combat, par.prior = T)
  save(myNorm.combat, file=myNorm.combat.file)
} else {
  load(myNorm.combat.file)
}

champ.SVD(beta=myNorm.cb)

png('normalized_data_MDSPlot.png')
plotMDS(myNorm.cb, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Slide))))
dev.off()

p <- innerheatmap(myNorm.cb[, order(myLoad$pd$Sample_Group)])
p

library(org.Hs.eg.db)
require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

detectCores <- function() {
  return(4)
}
