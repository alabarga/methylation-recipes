
### load file

myNorm.file <- "myNorm.Rdata"
load(myNorm.file)
load('sampleInfo.dat')

### BS vs OX
myDMP <- champ.DMP(beta = myNorm, arraytype = "EPIC", adjPVal=0.001)
hmc <- myDMP[[1]][abs(myDMP[[1]]$deltaBeta)>0,]
sel_dmp <- rownames(hmc)
p <- innerheatmap(myNorm[sel_dmp, order(myLoad$pd$Sample_Group)])
p

beta.BS <- myNorm[sel_dmp, which(pd$Sample_Group == 'BS')]
beta.oxBS <- myNorm[sel_dmp, which(pd$Sample_Group == 'OX')]

M.BS <- logit2(beta.BS)
M.oxBS <- logit2(beta.oxBS)

## Naive estimates

## The naive approach to obtain 5-hmC levels is $\beta_{BS} - \beta_{OxBS}$. This approach results in negative values for the 5-hmC levels.

hmC_naive <- beta.BS - beta.oxBS


## Using Ms
M_naive <- M.BS - M.oxBS
hmC_naive <- 2^M_naive / (2^M_naive + 1)


## OxBS.MLE: jointly estimation of 5mC and 5hmC levels from the paired experimental data
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5181539/
  
mle_estimation <- oxBS.MLE(beta.BS, beta.oxBS, N.BS, N.oxBS)

## MLML: consistent simultaneous estimates of DNA methylation and hydroxymethylation
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789553/

mlml_estimation <- MLML(T = M.BS , U = U.BS, L = U.oxBS, M = M.oxBS, exact=TRUE)



beta.OX <- beta.OX[, colnames(beta.BS)]
load("raw.Rdata")
load("BS.Rdata")
load("oxBS.Rdata")

targetsBS <- read.metharray.sheet(infoDir, pattern="BS.csv$")
rownames(targetsBS) <- targetsBS$Index

targetsOX <- read.metharray.sheet(dataDir, pattern="OX.csv$")
rownames(targetsOX) <- targetsOX$Sample_Name

cols <- targetsBS[colnames(beta.BS),'Sample_Name']
common <-  rownames(beta.BS)[ rownames(beta.BS) %in% rownames(beta.oxBS) ]

hmC_naive <- beta.BS[common, ] - beta.oxBS[common, targetsOX[cols, 'Index'] ]
C_naive <- 1-beta.BS
mC_naive <- beta.oxBS

## remove sample

cols <- cols[1:9]
hmC_naive <- hmC_naive[,1:9]
targets <- targetsBS[colnames(hmC_naive),]
colnames(hmC_naive) <- cols

dmp <- champ.DMP(beta = hmC_naive, pheno=targets$Sample_Group, adjPVal = 1)
plot(dmp$PD_to_CTRL$deltaBeta,-log10(dmp$PD_to_CTRL$adj.P.Val), pch=20)
multidensity(hmC_naive)
hist(dmp$PD_to_CTRL$adj.P.Val)

DMP.GUI(DMP=dmp$PD_to_CTRL, beta=hmC_naive, pheno=targets$Sample_Group )

p <- 0.00001
adjP <- 0.4
deltaBeta <- 0.1

dmpSel <- dmp$PD_to_CTRL[(abs(dmp$PD_to_CTRL$deltaBeta) > deltaBeta) & (dmp$PD_to_CTRL$adj.P.Val < adjP) ,]
cpgSel <- rownames(dmpSel)

## GREAT

write.table(EPIC.manifest.hg38[cpgSel,c('seqnames','start','end')], file='selected.cpgsPDvsCTRL.bed', quote=F, sep="\t", row.names=F, col.names=F)

cols_no <- cols[1:9]
clases <- as.factor(c(rep("BS",9), rep("oxBS",9)))
colnames(M.BS) <- targetsBS$Sample_Name
colnames(M.oxBS) <- targetsOX$Sample_Name

colnames(beta.BS) <- targetsBS$Sample_Name
colnames(beta.oxBS) <- targetsOX$Sample_Name

Mtotal <- cbind(M.BS[common, cols_no], M.oxBS[common, cols_no])
beta <- cbind(beta.BS[common, cols_no], beta.oxBS[common, cols_no])
myDMP <- champ.DMP(beta = beta, pheno=clases, adjPVal = 1)

plot(myDMP[[1]]$logFC,-log10(myDMP[[1]]$adj.P.Val), pch=20)

DMP.GUI(DMP=myDMP[[1]], beta=beta, pheno=clases)

save(hmC_naive, myDMP, file = "5hmC.Rdata")
load("5hmC.Rdata")

##### mCSEA

#infoDir <- "/home/alabarga/NavarraBiomed/analysis/parkinson/sampleInfo"
targetsBS <- read.metharray.sheet(infoDir, pattern="BS.csv$")
pheno <- targetsBS[,c('Sample_Group', 'Age','Sex')]
#colnames(pheno) <- c('expla', 'cov1', 'cov2')
rownames(pheno) <- targetsBS[,'Sample_Name']
pheno$Age <- as.numeric(pheno$Age)

myRank <- rankProbes(hmC_naive, pheno, refGroup = "CTRL", typeInput = "M", covariates = c('Sex'), continuous=c('Age'))

hmC_beta = 2^hmC_naive/(2^hmC_naive + 1)

myResults <- mCSEATest(myRank, hmC_beta, pheno, regionsTypes = "promoters", platform = "EPIC")
proms <- myResults[["promoters"]][-7]
attach(proms)

head(proms[order(padj),],20)
subset(proms, pval <0.001 & size > 10)

for (dmrName in rownames(subset(proms, pval <0.001 & size > 10))) {

  mCSEAPlot(myResults, regionType = "promoters", 
            dmrName = dmrName, genes=TRUE,
            transcriptAnnotation = "symbol", makePDF = FALSE)
  
}

detach(proms)

head(myResults[["promoters"]][,-7],20)

mCSEAPlotGSEA(myRank, myResults, regionType = "promoters", dmrName = "ACER3")
#####

cols <-colnames(hmC_naive)
cols_no <- cols[1:9]
selDMP <- myDMP[[1]][(myDMP[[1]]$logFC > 0) & (myDMP[[1]]$adj.P.Val < 0.2), ]
beta <- hmC_naive[rownames(selDMP), cols_no]

myDMP <- champ.DMP(beta = beta, pheno=targetsOX[cols_no, 'Sample_Group'], adjPVal = 1, arraytype = "EPIC")
DMP.GUI(DMP=myDMP[[1]], beta=beta, pheno=targetsOX[cols_no, 'Sample_Group'])

myDMPsel <- myDMP[[1]][myDMP[[1]]$P.Value < 0.05,]
myDMR <- champ.DMR(beta = beta, pheno=targetsOX[cols_no, 'Sample_Group'], method = "Bumphunter", B=100, minProbes=3, arraytype="EPIC")
myGSEA <- champ.GSEA(beta = hmC_naive[rownames(selDMP),], DMP=myDMPsel, DMR=NULL, arraytype = "EPIC")

myDMR <- champ.DMR(beta = hmC_naive, pheno=targetsOX[cols, 'Sample_Group'], method = "Bumphunter", B=100, adjPvalDmr=0.1, adjPvalProbe=0.8, arraytype="EPIC")


DMP.GUI(DMP=myDMP.2[[1]], beta=-hmC_naive[rownames(selDMP),], pheno=targetsOX[cols, 'Sample_Group'])

rgBS <-  preprocessFunnormRedGreen(RGsetBS, sex= 1*(targetsBS$Sex == "M"))
rgOxBS <-  preprocessFunnormRedGreen(RGsetOX, sex= 1*(targetsOX$Sex == "M"))


qc <- QCinfo(RGsetEx)

## BS vs OX

pD <- pData(mdat)
remove_samples <- c('201096090140_R07C01') 
mdat <- mdat[,!(pD$Index %in% remove_samples)]

pD <- pData(mdat)

BS <- pD$Index[pD$Methylation == "5mC"]
oxBS <- pD$Index[pD$Methylation == "5hmC"]

M.BS <- minfi::getMeth(mdat)[,BS]
U.BS <- getUnmeth(mdat)[,BS]

M.oxBS<- getMeth(mdat)[,oxBS]
U.oxBS <- getUnmeth(mdat)[,oxBS]

N.BS <- M.BS + U.BS
N.oxBS <- M.oxBS + U.oxBS

beta.BS <- beta[, BS]
beta.oxBS <- beta[, oxBS]

colnames(beta.BS) <- pD$Sample_Name[pD$Type == 'BS']
colnames(beta.oxBS) <- pD$Sample_Name[pD$Type == 'OX']
colnames(N.BS) <- pD$Sample_Name[pD$Type == 'BS']
colnames(N.oxBS) <- pD$Sample_Name[pD$Type == 'OX']


## Naive estimates

## The naive approach to obtain 5-hmC levels is $\beta_{BS} - \beta_{OxBS}$. This approach results in negative values for the 5-hmC levels.

hmC_naive <- beta.BS - beta.oxBS
C_naive <- 1-beta.BS
mC_naive <- beta_oxBS

naive_estimation <- 0

mle_estimation <- oxBS.MLE(beta.BS, beta.oxBS, N.BS, N.oxBS)

## MLML: consistent simultaneous estimates of DNA methylation and hydroxymethylation
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789553/

mlml_estimation <- MLML(T = M.BS , U = U.BS, L = U.oxBS, M = M.oxBS, exact=TRUE)


par(mfrow =c(1,3)) 

multidensity(mlml_estimation$mC, main="5-mC using MLML")
multidensity(mlml_estimation$hmC, main="5-hmC using MLML")
multidensity(mlml_estimation$C, main="5-C using MLML")


hmc <- mlml_estimation$hmC
hmC_naive <- M.BS - M.oxBS
pd.hmc <- pd[colnames(hmc),]
hmc <- hmC_naive



mvalues.limma <- logit2(hmc)
mvalues.limma <-  hmC_naive

betas <- 2^hmC_naive / (2^hmC_naive + 1)
betas[betas <= 0] <- 10^-10
betas[betas >= 1] <- 1- 10^-10

mvalues.limma <- betas
mod.limma = model.matrix( ~ -1 + pd.hmc$Status + pd.hmc$Sex + pd.hmc$Age)
colnames(mod.limma) <- c('CTRL','PD','Sex','Age')
cont.matrix <- makeContrasts(
  PDvsCTRL=PD-CTRL,
  levels=mod.limma)
fit.2 <- lmFit(mvalues.limma, mod.limma)
fit.2 <- contrasts.fit(fit.2, cont.matrix)
fit.2 <- eBayes(fit.2)
top <- topTable(fit.2,coef="PDvsCTRL",sort.by="p",number=nrow(mvalues.limma))
head(top)

sel_limma <- rownames(top)[1:1000]

p <- innerheatmap(hmC_naive[, order(clases)])
p



