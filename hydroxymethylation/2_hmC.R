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
myDMP <- champ.DMP(beta = myNorm, pheno = pd$Sample_Group, arraytype = "EPIC", compare.group = c("BS","OX"))
hmC.selected.cpgs <- rownames(myDMP[[1]][myDMP[[1]]$deltaBeta<0,])

plot(myNorm["cg07307802",order(pd$Sample_Group)], color = c("red","blue")[pd$Sample_Group], pcx=22)
plot(myNorm["cg07307802",order(pd$Sample_Group)], col = c("red","blue")[sort(as.factor(pd$Sample_Group))], pch=20)
legend('topright', legend=levels(as.factor(pd$Sample_Group)),pch=20,col=c("red","blue"))

library(ggplot2)

for (cpgname in head(hmC.selected.cpgs)) {
  print(cpgname)
  df = data.frame(samplename=pd$Sample_ID[order(pd$Sample_Group)], beta=myNorm[cpgname, order(pd$Sample_Group)], ox_bs=sort(pd$Sample_Group))
  g = ggplot(df, aes(samplename, beta, col=ox_bs)) + 
          geom_point() + scale_color_discrete(name="OX/BS") + ggtitle(cpgname)
  print(g)
}

g = ggplot(df, aes(ox_bs, beta)) + 
  geom_boxplot(colour = "black", fill = "#56B4E9")  + geom_jitter() + ggtitle(cpgname) 
print(g )

DMP.GUI(DMP=myDMP[[1]][myDMP[[1]]$deltaBeta<0,],
        beta=myNorm[,order(pd$Sample_Group)],
        pheno=sort(myLoad$pd$Sample_Group),
        cutgroupnumber=4)

myDMP.file <- "myDMP.Rdata"
if(!file.exists(myDMP.file)){
  myDMP <- champ.DMP(beta = myNorm, pheno = pd$Sample_Group, arraytype = "EPIC", compare.group = c("OX","BS"))
  save(myDMP, file=myDMP.file)
} else {
  load(myDMP.file)
}

bs <-getBeta(mset.noob)[, mset$Sample_Group=="BS"]
clases <- mset$Status[mset$Sample_Group=="BS"]
densityPlot(bs, sampGroups = clases, main = paste("Density plot of raw data (", nrow(bs), " probes)", sep = ""), xlab = "Beta")

bs.all <- myNorm[, pd$Sample_Group=="BS"]
ox.all <- myNorm[, pd$Sample_Group=="OX"]
hmc.all <- bs.all-ox.all
quantile(hmc.all[hmc.all<0],probs=c(0.05))

bs.all <- subset(myNorm, subset=(rownames(myNorm) %in% rownames(myDMP$BS_to_OX)), pd$Sample_Group=="BS")
ox.all <- subset(myNorm, subset=(rownames(myNorm) %in% rownames(myDMP$BS_to_OX)), pd$Sample_Group=="OX")
hmc.all <- bs.all-ox.all
quantile(hmc.all[hmc.all<0],probs=c(0.05))


densityPlot(hmc.all, sampGroups = clases, main = paste("Density plot of raw data (", nrow(hmc.all), " probes)", sep = ""), xlab = "Beta")

bs <- subset(myNorm,subset=(rownames(myNorm) %in% hmC.selected.cpgs), pd$Sample_Group=="BS")
ox <- subset(myNorm,subset=(rownames(myNorm) %in% hmC.selected.cpgs), pd$Sample_Group=="OX")
colnames(bs) <- subset(pd,pd$Sample_Group=="BS",c("Sample_ID"))$Sample_ID
colnames(ox) <- subset(pd,pd$Sample_Group=="OX",c("Sample_ID"))$Sample_ID
hmc <- bs-ox
clases <- subset(pd,pd$Sample_Group=="BS",c("Status"))$Status


hmc.file <- "5hmc.Rdata"
save(bs, ox, hmc, clases, file=hmc.file)
load(hmc.file)

densityPlot(hmc, sampGroups = (clases), main = paste("Density plot of raw data (", nrow(hmc), " probes)", sep = ""), xlab = "Beta")

p <- innerheatmap(hmc[,order(clases)])
p

M <- minfi::getMeth(myLoad$mset)
U <- minfi::getUnmeth(myLoad$mset)

mset.noob<- preprocessNoob(myLoad$rgSet)
M <- minfi::getMeth(mset.noob)
U <- minfi::getUnmeth(mset.noob)

mset.ENmix <- preprocessENmix(myLoad$rgSet)
M <- minfi::getMeth(mset.ENmix)
U <- minfi::getUnmeth(mset.ENmix)

save(mset.noob,mset.ENmix, file='msets.Rdata')
load('msets.Rdata')

MethylatedBS <- subset(M,select=(pd$Sample_Group=="BS"))
UnMethylatedBS <- subset(U,select=(pd$Sample_Group=="BS"))
MethylatedOxBS <- subset(M,select=(pd$Sample_Group=="OX"))
UnMethylatedOxBS <- subset(U,select=(pd$Sample_Group=="OX"))
colnames(MethylatedBS) <- mset$Sample_ID[seq(1,18,2)]
colnames(MethylatedOxBS)<- mset$Sample_ID[seq(1,18,2)]
colnames(UnMethylatedOxBS)<- mset$Sample_ID[seq(1,18,2)]
colnames(UnMethylatedBS)<- mset$Sample_ID[seq(1,18,2)]

save(MethylatedBS, UnMethylatedBS, MethylatedOxBS, UnMethylatedOxBS, file="MLML.ENmix.Rdata" )
load("MLML.ENmix.Rdata")
results_exact <- MLML(T.matrix = MethylatedBS , U.matrix = UnMethylatedBS,
                      L.matrix = UnMethylatedOxBS, M.matrix = MethylatedOxBS)

densityPlot(results_exact$hmC, sampGroups = (clases), main = paste("Density plot of raw data (", nrow(results_exact$hmC), " probes)", sep = ""), xlab = "Beta")
p <- innerheatmap(results_exact$hmC[,order(clases)])

results_em <- MLML(T.matrix = MethylatedBS , U.matrix = UnMethylatedBS,
                   L.matrix = UnMethylatedOxBS, M.matrix = MethylatedOxBS, iterative = TRUE)


beta.BS <- MethylatedBS / (MethylatedBS + UnMethylatedBS)
beta.oxBS <- MethylatedBS / (MethylatedOxBS + UnMethylatedOxBS)
N.BS <- (MethylatedBS + UnMethylatedBS)
N.oxBS <- (MethylatedOxBS + UnMethylatedOxBS)

colnames(beta.BS) <- mset$Sample_ID[seq(1,18,2)]
colnames(beta.oxBS)<- mset$Sample_ID[seq(1,18,2)]
colnames(N.BS)<- mset$Sample_ID[seq(1,18,2)]
colnames(N.oxBS )<- mset$Sample_ID[seq(1,18,2)]

mle_estimation <- oxBS.MLE(beta.BS,beta.oxBS,N.BS,N.oxBS)
densityPlot(mle_estimation$"5hmC", sampGroups = (clases), main = paste("Density plot of raw data (", nrow(mle_estimation$"5hmC"), " probes)", sep = ""), xlab = "Beta")

hmC <- results_exact$hmC #[rownames(hmc),]

hmc.mean <- rowMeans(hmc)
hmC.mean <- rowMeans(hmC)

plot(hmc.mean, hmC.mean)

hmc.diff <- rowMeans(hmc[,clases=="PD"]) - rowMeans(hmc[,clases=="CTRL"]) 
hmC.diff <- rowMeans(hmC[,clases=="PD"]) - rowMeans(hmC[,clases=="CTRL"]) 

plot(hmc.mean[hmC.pd.selected.cpgs], hmC.mean[hmC.pd.selected.cpgs])

dmp.hmc <- champ.DMP(beta = hmc, pheno = clases, arraytype = "EPIC", compare.group = c("PD","CTRL"))
dmp.hmC <- champ.DMP(beta = hmC, pheno = clases, arraytype = "EPIC", compare.group = c("PD","CTRL"))

hmC[hmC <= 0] <- 10^-12

M.limma <- log2(hmC/(1-hmC))
M.limma <- lumi::beta2m(hmC)

mod.limma = model.matrix( ~ 0 + clases)
colnames(mod.limma) <- c('CTRL','PD')
cont.matrix <- makeContrasts(
  PDvsCTRL=PD-CTRL,
  levels=mod.limma)
fit.2 <- lmFit(M.limma, mod.limma)
fit.2 <- contrasts.fit(fit.2, cont.matrix)
fit.2 <- eBayes(fit.2)
top <- topTable(fit.2,coef="PDvsCTRL",sort.by="p",number=nrow(hmC))
head(top)

hmc[hmc <= 0] <- 10^-12
M2.limma <- lumi::beta2m(hmc)
mod.limma = model.matrix( ~ 0 + clases)
colnames(mod.limma) <- c('CTRL','PD')
cont.matrix <- makeContrasts(
  PDvsCTRL=PD-CTRL,
  levels=mod.limma)
fit.2 <- lmFit(M2.limma, mod.limma)
fit.2 <- contrasts.fit(fit.2, cont.matrix)
fit.2 <- eBayes(fit.2)
top2 <- topTable(fit.2,coef="PDvsCTRL",sort.by="p",number=nrow(hmc))
head(top2)


















##########################################################################################################33
myLoad.file <- "myLoad.Rdata"
if(!file.exists(myLoad.file)){
  targets <- read.metharray.sheet(dataDir)
  
  RGset <- read.metharray.exp(base = dataDir, targets = targets)
  RGsetEx <- read.metharray.exp(targets = targets, extended=T)
  
  mset <- preprocessRaw(RGset)
  
  save(RGsetEx, targets, mset, file=myLoad.file)
} else {
  load(myLoad.file)
}

RGset <- RGsetEx
qcReport(RGset, sampNames= targets$Sample_Name, sampGroups= targets$Sample_Group, pdf= "qcReport.rgset.pdf")

myNorm.file <- "myNorm.Rdata"
if(!file.exists(myNorm.file)){
  detP<- detectionP(RGset)
  failed<- detP > 0.01
  ## Keep probes which failed in at most this many arrays (0 = the probe passed in all arrays)
  maxFail<- 0
  mset <- preprocessRaw(RGset)
  mset<- mset[rowSums(failed) <= maxFail, ]
  MSet.swan<- preprocessSWAN(RGset, mSet = mset)
  
  save(MSet.swan, file=myNorm.file)
} else {
  load(myNorm.file)
}

myNorm.file <- "myNorm.Rdata"
if(!file.exists(myNorm.file)){
  myNorm <- champ.norm(beta=myLoad$beta,mset=myLoad$mset,rgSet=myLoad$rgSet,method="FunctionalNormalization",arraytype="EPIC")
  save(myNorm, file=myNorm.file)
} else {
  load(myNorm.file)
}


pdf("densityPlot.msetSwan.pdf")
densityPlot(MSet.swan, sampGroups= targets$Sample_Group, main=sprintf("Beta values for filtered probes (n= %s)", nrow(MSet.swan)))
dev.off()

b<- getBeta(MSet.swan)
pcaResult<-prcomp(t(b))

plot(pcaResult$x,
     main= 'Principal components from most variable beta values',
     xlab= sprintf('PC1 (sd: %s%%)',
                     round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
     ylab= sprintf('PC2 (sd: %s%%)',
                    round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
     pch= 19, col= ifelse(targets$Sample_Group=="BS", "blue", "red")
)
grid(col= 'grey30')
text(x= pcaResult$x[, "PC1"] + 4,
     y= pcaResult$x[, "PC2"],
     labels= targets$Sample_Name, cex= 1)


M<- getM(MSet.swan, type= "beta", betaThreshold = 0.001)

dmp_BS <- dmpFinder(M, pheno=targets$Sample_Group, type="categorical", shrinkVar= TRUE)

mapped <- mapToGenome(MSet.swan)
probPos<- as.data.frame(mapped@rowData)
probPos$probe_id<- rownames(probPos)

cpgs.hmC <- read.table('/home/alabarga/NavarraBiomed/analysis/parkinson/results/2_hmC/selected.cpgs.5hmC.txt')
write.table(EPIC.hg19[cpgs.hmC$V1,c('seqnames','start','end','cpgname')], file='selected.cpgs.5hmC.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)

cpgs.hmC.pd <- read.table("~/NavarraBiomed/analysis/parkinson/results/3_hmC_PD_vs_CTRL/selected.cpgs.5hmC.PDvsCTRL.hg19.bed")

bedfile <- read.table('/home/alabarga/NavarraBiomed/analysis/parkinson/results/2_hmC/selected.cpgs.5hmc.bed')
my.granges <- GRanges(seqnames = bedfile$V1, ranges = IRanges(start = bedfile$V2, end = bedfile$V3))
o <- findOverlaps(my.granges,EPIC.hg19.manifest)
hmC.selected.cpgs <- names(EPIC.hg19.manifest)[unique(subjectHits(o))]

write.table(EPIC.hg19[hmC.selected.cpgs,c('seqnames','start','end','cpgname')], file='selected.cpgs.5hmC.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)

cpgs.bs.pd <- read.table("~/NavarraBiomed/analysis/parkinson/results/3_BS_PD_vs_CTRL/selected.cpgs.BS.PDvsCTRL.001.hg19.bed")

hmC.selected.cpgs # 47904
hmC.pd.selected.cpgs # 2450
bs.pd.selected.cpgs # 1429
common.selected.cpgs # 66
common.pd.selected.cpgs # 9

save(hmC.selected.cpgs, hmC.pd.selected.cpgs, bs.pd.selected.cpgs, common.selected.cgps, common.pd.selected.cpgs, file="selected.cpgs.lists.Rdata")
load("selected.cpgs.lists.Rdata")

