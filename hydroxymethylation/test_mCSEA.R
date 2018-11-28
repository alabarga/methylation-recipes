data(mcseadata)
library(mCSEA)


myRank <- rankProbes(betaTest, phenoTest, refGroup = "Control")
myResults <- mCSEATest(myRank, betaTest, phenoTest, 
                       regionsTypes = "promoters", platform = "EPIC")
mCSEAPlot(myResults, regionType = "promoters", 
          dmrName = "CLIC6",
          transcriptAnnotation = "symbol", makePDF = FALSE)

library(mCSEA)
data(mcseadata)
data(precomputedmCSEA)
mCSEAPlot(myResults, "promoters", "CLIC6", transcriptAnnotation = "symbol", makePDF = FALSE)
