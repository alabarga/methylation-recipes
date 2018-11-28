# pipeline minfi

targetsBS <- read.metharray.sheet(dataDir, pattern="BS.csv$") # BS only
targetsOX <- read.metharray.sheet(dataDir, pattern="OX.csv$") # OX only
source('~/NavarraBiomed/R/meth/EPIC/pipeline_minfi.R')

# RGset <- read.metharray.exp(targets = targets)
RGsetEx <- read.metharray.exp(targets = targetsBS, extended=T)

M.BS <- minfi::getMeth(mdat)
U.BS <- minfi::getUnmeth(mdat)
N.BS <- M.BS + U.BS
beta.BS <- getBeta(mdat)

M.oxBS <- minfi::getMeth(mdat)
U.oxBS <- minfi::getUnmeth(mdat)
N.oxBS <- M.oxBS + U.oxBS
beta.oxBS <- getBeta(mdat)

save(M.BS, U.BS, N.BS, beta.BS, file="BS.Rdata" )
save(M.oxBS, U.oxBS, N.oxBS, beta.oxBS, file="oxBS.Rdata" )
save(RGsetEx, GMsetEx, mraw, mdat, qc, anno, file="raw.Rdata")

rm(RGsetEx)
rm(GMsetEx)
rm(anno)
rm(qc)
gc()



















cc <- estimateCellCounts(RGsetEx, returnAll=TRUE)
cc$counts
library(reshape)
counts <- melt(cc$counts)
colnames(counts) <-c('Sample','CellType','Percentage')

ggplot(data = counts, aes(x = Sample, y = CellType)) +
  geom_tile(aes(fill = Percentage)) 

source("~/NavarraBiomed/R/meth/EPIC.R")
load("~/NavarraBiomed/parkinson/session.RData")

## # Principal component regression analysis plot
cov<-data.frame(group=pData(mdat)$Sample_Group, slide=factor(pData(mdat)$Slide))
pcrplot(beta, cov, npc=6)

## #filter out low quality and outlier data points for each probe;
## #rows and columns with too many missing value can be removed 
## #if specify; Do imputation to fill missing data if specify.
beta <- rm.outlier(beta, qcscore=qc, rmcr=TRUE, impute=TRUE)

## #Non-negative control surrogate variables
sva<-ctrlsva(RGsetEx)

## ANALYSIS

rownames(targets) <- targets$Index
sampleNames <- pData(mdat)$Sample_Name
clases <- as.factor(pData(mdat)$Sample_Group)
colnames(beta) <- sampleNames

plotMDS(beta, labels=sampleNames, col=as.integer(clases))

QC.GUI(beta=beta,
       pheno=pData(mdat)$Type, arraytype = arraytype,
       arraytype="EPIC")