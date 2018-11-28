
if(!file.exists("EPIC.manifest.hg38.bed")){
  
  data(EPIC.manifest.hg38)
  write.table(EPIC.manifest.hg38[,c('seqnames','start','end')], file='EPIC.manifest.hg38.bed', quote=F, sep="\t", row.names=F, col.names=F)
}

EPIC.hg19 <- as.data.frame(EPIC.hg19.manifest)
EPIC.hg19$cpgname <- rownames(EPIC.hg19)
annotationEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotationEPIC$cpgname <- rownames(annotationEPIC)

if(!file.exists("EPIC.hg19.manifest.bed")){
  # df <- data.frame(chr=seqnames((EPIC.hg19.manifest)),start=start(ranges(EPIC.hg19.manifest)), end=end(ranges(EPIC.hg19.manifest)))
  # write.table(df, file='EPIC.manifest.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)
  
  write.table(EPIC.hg19[,c('seqnames','start','end','cpgname')], file='EPIC.manifest.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)
  
}

#ez<-getMappedEntrezIDs(cpgs, array.type = "EPIC")
#unlist(mget(ez$sig.eg, org.Hs.egSYMBOL))

