
if ((!conf.minfi) & ((conf.norm == "SWAN") | (conf.norm == "FunctionalNormalization"))) {
    stop("You need old minfi style data loading for 'SWAN' and 'FunctionalNormalization'")
}

### LOAD DATA
targets <- read.metharray.sheet(dataDir)
if (conf.minfi) {
    myLoad <- champ.load(dataDir, method="minfi", arraytype="EPIC") # loaded as minfi way OLD
} else {
    myLoad <- champ.load(directory = dataDir, arraytype = "EPIC")     
}

pd <- myLoad$pd
champ.QC()
CpG.GUI(arraytype="EPIC")
QC.GUI()

### NORMALIZATION
# four different normalization mehods
if (conf.norm == "BMIQ") myNorm <- champ.norm(arraytype="EPIC")
if (conf.norm == "PBC") myNorm <- champ.norm(method="PBC", arraytype = "EPIC")
if (conf.norm == "SWAN") myNorm <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="SWAN",arraytype="EPIC")
if (conf.norm == "FunctionalNormalization") myNorm <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="FunctionalNormalization",arraytype="EPIC")

QC.GUI(beta=myNorm, arraytype="EPIC")
p <- innerheatmap(myNorm[, order(myLoad$pd$Sample_Group)])
p

save(myNorm, file="myNorm.BMIQ.Rdata")
load("myNorm.Rdata")

### BATCH CORRECTION
myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
champ.SVD(beta=myNorm)

mod.combat = model.matrix( ~ 1 + myLoad$pd$Sample_Group + myLoad$pd$Sex + myLoad$pd$Age)
myNorm <- sva::ComBat(dat = as.matrix(myNorm) , batch = myLoad$pd$Slide, mod=mod.combat, par.prior = T)
save(myNorm, file="myNorm.CB.Rdata")
load("myNorm.CB.Rdata")


### BS vs OX
myDMP <- champ.DMP(beta = myNorm, arraytype = "EPIC", adjPVal=0.001)
hmc <- myDMP[[1]][abs(myDMP[[1]]$deltaBeta)>0,]
sel_dmp <- rownames(hmc)
p <- innerheatmap(myNorm[sel_dmp, order(myLoad$pd$Sample_Group)])
p

mvalues.limma <- myNorm
mod.limma = model.matrix( ~ -1 + myLoad$pd$Sample_Group + myLoad$pd$Sex + myLoad$pd$Status + myLoad$pd$Age)
colnames(mod.limma) <- c('BS','OX','Sex','Disease','Age')
cont.matrix <- makeContrasts(
OXvsBS=OX-BS,
levels=mod.limma)
fit.2 <- lmFit(mvalues.limma, mod.limma)
fit.2 <- contrasts.fit(fit.2, cont.matrix)
fit.2 <- eBayes(fit.2)
top <- topTable(fit.2,coef="OXvsBS",sort.by="logFC",number=nrow(mvalues.limma))
sel_probes <- rownames(top[(top$adj.P.Val < 0.001) & (top$logFC > 0),])

hist(top$adj.P.Val)
hist(top$logFC)
heatmap(myNorm[sel_probes, order(myLoad$pd$Sample_Group)], Rowv=NA, Colv=NA)
plotMDS(myNorm[sel_probes, order(myLoad$pd$Sample_Group)], labels=sort(myLoad$pd$Sample_Group), col=as.numeric(as.factor(sort(myLoad$pd$Sample_Group))))

save.image("~/NavarraBiomed/analysis/parkinson/data/session2.RData")
load("~/NavarraBiomed/analysis/parkinson/data/session2.RData")


### 5hmC

targets.BS <- targets[targets$Sample_Group == 'BS',]
targets.OX <- targets[targets$Sample_Group == 'OX',]
beta.BS <- myNorm[,targets.BS$Sample_Name]
beta.oxBS<- myNorm[,targets.OX$Sample_Name]

hmC_naive <- beta.BS[sel_probes, ] - beta.oxBS[sel_probes, ]
C_naive <- 1 - beta.BS
mC_naive <- beta.oxBS

### PD vs CTRL

mvalues.limma <- hmC_naive
mod.limma = model.matrix( ~ -1 + targets.BS$Status + targets.BS$Sex +  + targets.BS$Age)
colnames(mod.limma) <- c('CTRL','PD','Sex','Age')
cont.matrix <- makeContrasts(
  PDvsCTRL=PD-CTRL,
  levels=mod.limma)

fit.2 <- lmFit(mvalues.limma, mod.limma)
fit.2 <- contrasts.fit(fit.2, cont.matrix)
fit.2 <- eBayes(fit.2)

top.pd <- topTable(fit.2,coef="PDvsCTRL",sort.by="logFC",number=nrow(mvalues.limma))
sel_probe.pd <- rownames(top.pd[(top.pd$adj.P.Val < 0.2) & (top.pd$logFC > 0),])
cpgs <- as.matrix(hmC_naive[sel_probe.pd, order(targets.BS$Status)])

hist(top.pd$adj.P.Val)
hist(top.pd$logFC)
heatmap(cpgs, Rowv=NA, Colv=NA)
plotMDS(cpgs, labels=sort(targets.BS$Status), col=as.numeric(as.factor(sort(targets.BS$Status))))

save.image("~/NavarraBiomed/analysis/parkinson/data/session3.RData")
load("~/NavarraBiomed/analysis/parkinson/data/session3.RData")

### GREAT

EPIC.manifest.hg37 <- read.csv('/home/alabarga/NavarraBiomed/R/meth/EPIC/EPIC.hg19.manifest.tsv', sep = "\t")
rownames(EPIC.manifest.hg37) <- EPIC.manifest.hg37$probeID


write.table(EPIC.manifest.hg37[sel_probe.pd, c('CpG_chrm','CpG_beg','CpG_end')], file='selected.cpgsPDvsCTRL.bed', quote=F, sep="\t", row.names=F, col.names=F)

hmc.probes <- rownames(hmC_naive)
write.table(EPIC.manifest.hg37[hmc.probes, c('CpG_chrm','CpG_beg','CpG_end')], file='selected.cpgs.5hmc.bed', quote=F, sep="\t", row.names=F, col.names=F)

### 
write.table(EPIC.manifest.hg38[cpgSel,c('seqnames','start','end')], file='selected.cpgsPDvsCTRL.bed', quote=F, sep="\t", row.names=F, col.names=F)


### DMR

design <- model.matrix(~ 0 + targets.BS$Status + targets.BS$Sex +  + targets.BS$Age)
colnames(design) <- c("CTRL", "PD", 'SEX', 'AGE')
design[,'AGE'] <- as.numeric(design[,'AGE'])

cont.matrix <- makeContrasts(PDvsCTRL=PD-CTRL, levels=design)
myannotation <- cpg.annotate("array", cpgs, analysis.type="differential", cont.matrix=cont.matrix, design=design, contrasts = TRUE, fdr = 1, coef="PDvsCTRL", what="Beta", arraytype = "EPIC")

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, pcutoff = 1)
dmrs <- dmrcoutput$results[abs(dmrcoutput$results$maxbetafc)>0.1,]
dim(dmrs)


# three DMR methods

conf.group <- targets.BS$Status
conf.compare <- c("BS","OX")

if (conf.dmr == "Bumphunter") myDMR <- champ.DMR(beta=cpgs, pheno=conf.group, arraytype = "EPIC")
if (conf.dmr == "DMRcate") myDMR <- champ.DMR(beta=cpgs, pheno=conf.group, arraytype = "EPIC", method="DMRcate", cores=1)
if (conf.dmr == "ProbeLasso") myDMR <- champ.DMR(beta=cpgs, pheno=conf.group, arraytype = "EPIC", method="ProbeLasso", compare.group=conf.compare)

DMP.GUI()
DMP.GUI(beta=myNorm[, order(myLoad$pd$Sample_Group)], pheno=sort(myLoad$pd$Sample_Group))

sel_probes <- rownames(myDMP$PD_to_CTRL)
p <- innerheatmap(myNorm[sel_probes, order(myLoad$pd$Sample_Group)])
p

save(myDMP, file="myDMP.Rdata")
load("myDMP.Rdata")

# three DMR methods

if (conf.dmr == "Bumphunter") myDMR <- champ.DMR(beta= arraytype = "EPIC")
if (conf.dmr == "DMRcate") myDMR <- champ.DMR(arraytype = "EPIC", method="DMRcate", cores=1)
if (conf.dmr == "ProbeLasso") myDMR <- champ.DMR(arraytype = "EPIC", method="ProbeLasso", compare.group=conf.compare)

DMR.GUI(DMR=myDMR, arraytype="EPIC", compare.group=conf.compare)

save(myDMR, file="myDMR")
load("myDMR")

myGSEA <- champ.GSEA(DMP=myDMP[[1]],arraytype = "EPIC")

myBlock <- champ.Block(arraytype = "EPIC")
Block.GUI(arraytype="EPIC",compare.group=c("PrEC_cells","LNCaP_cells"))

myEpiMod <- champ.EpiMod(arraytype="EPIC")
myrefbase <- champ.refbase(arraytype = "EPIC")
myCNA <- champ.CNA(control = F,arraytype = "EPIC")




