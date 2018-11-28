##
## BS_PD_vs_CTRL for BS only
##

source("~/NavarraBiomed/analysis/parkinson/code/experiment.r")

dataDir <- "~/NavarraBiomed/analysis/parkinson/data/3_BS_PD_vs_CTRL"
workDir <- "~/NavarraBiomed/analysis/parkinson/results/3_BS_PD_vs_CTRL"
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
  myNorm <- champ.norm(beta=myLoad$beta,mset=myLoad$mset,rgSet=myLoad$rgSet,method="FunctionalNormalization",arraytype="EPIC")
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

DMP.GUI(DMP=myDMP$PD_to_CTRL, beta=myNorm, pheno=pd$Status )

plotVolcano(myDMP$PD_to_CTRL)
plotManhattan(myDMP$PD_to_CTRL)

p <- 0.0001
adjP <- 0.4
deltaBeta <- 0.1

dmps <- subset(myDMP[[1]],abs(myDMP[[1]]$deltaBeta)>deltaBeta & (myDMP[[1]]$P.Value)<p )
sel_dmp <- rownames(dmps)

g <- volcanoPlot(myDMP[[1]], pval=p, logFC=deltaBeta, N=1000)

myDM <- myNorm[sel_dmp,]
plotMDS(myDM, labels=(pd$Sample_Name), col=(as.factor((pd$Status))), main="DMP (Disease)")

plotHeat(myNorm[sel_dmp, order(pd$Status)])
plotHeat(myLoad$beta[sel_dmp, order(pd$Status)])

write.table(EPIC.hg19[sel_dmp,c('seqnames','start','end','cpgname')], file='selected.cpgs.BS.PDvsCTRL.001.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)

plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Sample_Group))), main="RAW (BSvsOX)")
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Slide))), main="RAW (Slide)")
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Status))), main="RAW (Disease)")
plotMDS(myNorm, labels=(pd$Sample_Name), col=(as.factor((pd$Status))), main="After Funnorm (Disease)")

png('original_data_MDSPlot.png')
plotMDS(myLoad$beta, labels=(pd$Sample_Name), col=as.numeric(as.factor((pd$Sample_Group))))
dev.off()

p <- innerheatmap(myNorm[sel_dmp, order(pd$Status)])
p
