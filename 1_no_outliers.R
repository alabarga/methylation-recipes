source("~/NavarraBiomed/analysis/parkinson/code/experiment.r")

dataDir <- "~/NavarraBiomed/analysis/parkinson/data/1_no_outliers"
workDir <- "~/NavarraBiomed/analysis/parkinson/results/1_no_outliers"

setwd(workDir)

myLoad.file <- "myLoad.Rdata"
if(!file.exists(myLoad.file)){
  myLoad <- champ.load(directory = dataDir, arraytype = "EPIC", method="minfi")
  save(myLoad, file=myLoad.file)
} else {
  load(myLoad.file)
}

champ.QC()
plotBetasByType(myLoad$mset)
plotBetasByType(myNorm[,1],probeTypes=annotationEPIC[,c('Name','Type')])

myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
pd <- myLoad$pd

mdsPlot(myLoad$beta, numPositions = 3000, sampGroups = pd$Status, sampNames = pd$Sample_Name)

plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Sample_Group))), main="MDS (BSvsOX)")
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Slide))), main="MDS (Slide)")
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Status))), main="MDS (Disease)")

png('original_data_MDSPlot.png')
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Sample_Group))))
dev.off()

pd <- myLoad$pd
champ.QC()

norm.method="BMIQ"
myNorm.file <- paste("myNorm",norm.method,"Rdata", sep=".")
if(!file.exists(myNorm.file)){
  myNorm <- champ.norm(arraytype="EPIC", method=norm.method)
  save(myNorm, file=myNorm.file)
} else {
  load(myNorm.file)
}
champ.QC(beta=myNorm, pheno = pd$Sample_Group, PDFplot=FALSE)

mdsPlot(myNorm, numPositions = 3000, sampNames = pd$Sample_Name, sampGroups = pd$Status)
mdsPlot(myNorm, numPositions = 3000, sampNames = pd$Sample_Name, sampGroups = pd$Sample_Group)
densityPlot(myNorm, sampGroups = pd$Sample_Group, main = paste("Density plot of raw data (", nrow(myNorm), " probes)", sep = ""), xlab = "Beta")


plotMDS(myNorm, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Sample_Group))), main="MDS (BSvsOX)")
plotMDS(myNorm, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Slide))), main="MDS (Slide)")
plotMDS(myNorm, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Status))), main="MDS (Disease)")
plotMDS(myNorm, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Array))), main="MDS (Array)")

p <- innerheatmap(myNorm[, order(myLoad$pd$Sample_Group)])
p

