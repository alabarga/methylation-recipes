
# Create new preprocessFunNorm() function
preprocessFunnormRedGreen <- function (rgSet, nPCs = 2, sex = NULL, verbose = TRUE)
{
  #minfi:::.isRG(rgSet)
  rgSet <- updateObject(rgSet)
  if (verbose)
    cat("[preprocessFunnorm] Mapping to genome\n")
  gmSet <- mapToGenome(rgSet)
  subverbose <- max(as.integer(verbose) - 1L, 0)
  if (verbose)
    cat("[preprocessFunnorm] Quantile extraction\n")
  extractedData <- minfi:::.extractFromRGSet450k(rgSet)
  if (is.null(sex)) {
    gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
    sex <- rep(1L, length(gmSet$predictedSex))
    sex[gmSet$predictedSex == "F"] <- 2L
  }
  rm(rgSet)
  if (verbose)
    cat("[preprocessFunnorm] Normalization\n")
  CN <- getCN(gmSet)
  minfi:::.normalizeFunnorm450k(object = gmSet, extractedData = extractedData,
                                sex = sex, nPCs = nPCs, verbose = subverbose)
  
}

detectCores <- function(){
  4
}

gr2bed<- function(cpgs, filename='foo.bed') {
  
  write.table(EPIC.manifest.hg38[cpgs,c('seqnames','start','end')], file=filename, quote=F, sep="\t", row.names=F, col.names=F)
  # to write that to foo.bed. The only trick is remembering the BED uses 0-based coordinates. 
  # If you have a GRangesList rather than a GRanges object, just use unlist(gr) in place of gr (things should still be in the same order).
  
}
