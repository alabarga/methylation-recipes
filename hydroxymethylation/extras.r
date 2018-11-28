# minfi

baseDir <- "~/NavarraBiomed/"
projectDir <- "~/NavarraBiomed/analysis/parkinson/"




dataDir <- "~/NavarraBiomed/analysis/parkinson/data"
infoDir <- "/home/alabarga/NavarraBiomed/analysis/parkinson/sampleInfo"
workDir <- "~/NavarraBiomed/analysis/parkinson/results"
setwd(workDir)

conf.minfi <- FALSE
conf.norm <- "BMIQ"
conf.dmr <- "DMRcate"
conf.adjPVal = 0.5

# pipeline ChAMP
conf.group <- ""
conf.compare <- c("BS","OX")
workDir <- "~/NavarraBiomed/analysis/parkinson/results/results_5hmC"
setwd(workDir)

conf.group <- ""
conf.compare <- c("CTRL","PD")
workDir <- "~/NavarraBiomed/analysis/parkinson/results"
setwd(workDir)

source('~/NavarraBiomed/R/meth/EPIC/pipeline_ChAMP.R')