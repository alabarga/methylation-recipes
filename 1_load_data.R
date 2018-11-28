source("~/NavarraBiomed/analysis/parkinson/code/experiment.r")

dataDir <- "~/NavarraBiomed/analysis/parkinson/data/1_no_outliers"
workDir <- "~/NavarraBiomed/analysis/parkinson/results/2_hmC"

targets <- read.metharray.sheet(dataDir, pattern="*.csv")
RGsetEx <- read.metharray.exp(targets = targets, extended=T)

## step-by-step
mraw <- preprocessRaw(RGsetEx)
b <-getBeta(mraw, "Illumina") # type="Illumina" sets offset=100 as per Genome Studio

M <- minfi::getMeth(mraw)
U <- minfi::getUnmeth(mraw)
N <- M + U

myNorm <- champ.norm(beta=b, arraytype="EPIC")

qc <- QCinfo(RGsetEx, outlier=FALSE)
mdat<-preprocessENmix(RGsetEx, bgParaEst="oob", dyeCorr="RELIC", QCinfo=qc, nCores=4)
b <-getBeta(mdat, "Illumina")
myNorm <- champ.norm(beta=b, arraytype="EPIC")


