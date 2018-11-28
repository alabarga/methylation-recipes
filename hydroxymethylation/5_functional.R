
workDir <- "~/NavarraBiomed/analysis/parkinson/results/5_functional_analysis"
setwd(workDir)

load("selected.cpgs.lists.Rdata")

hmC.selected.cpgs # 47904
hmC.pd.selected.cpgs # 2450
bs.pd.selected.cpgs # 1429
common.selected.cpgs # 66
common.pd.selected.cpgs # 9
all.probes <- EPIC.hg19$cpgname 

table(annotationEPIC[hmC.pd.selected.cpgs,]$chr)
table(annotationEPIC[hmC.selected.cpgs,]$chr)
table(annotationEPIC[bs.pd.selected.cpgs,]$chr)

hmC.selected.cpgs <- sort(as.character(hmC.selected.cpgs ))
hmC.pd.selected.cpgs <- sort(as.character(hmC.pd.selected.cpgs ))
bs.pd.selected.cpgs <-  sort(as.character(bs.pd.selected.cpgs ))

save(all.probes, hmC.selected.cpgs, hmC.pd.selected.cpgs, bs.pd.selected.cpgs, common.selected.cpgs, common.pd.selected.cpgs, file="selected.cpgs.lists.Rdata" )

EPIC.hg19 <- as.data.frame(EPIC.hg19.manifest)
EPIC.hg19$cpgname <- rownames(EPIC.hg19)

EPIC.hg19 <- EPIC.hg19[sort(rownames(EPIC.hg19)),]

write.table(EPIC.hg19[,c('seqnames','start','end','cpgname')], file='EPIC.manifest.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)

write.table(EPIC.hg19[hmC.selected.cpgs,c('seqnames','start','end','cpgname')], file='selected.cpgs.5hmC.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)
write.table(EPIC.hg19[hmC.pd.selected.cpgs,c('seqnames','start','end','cpgname')], file='selected.cpgs.5hmC.pd.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)
write.table(EPIC.hg19[bs.pd.selected.cpgs,c('seqnames','start','end','cpgname')], file='selected.cpgs.bs.pd.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)
write.table(EPIC.hg19[common.selected.cpgs,c('seqnames','start','end','cpgname')], file='selected.common.cpgs.bed', quote=F, sep="\t", row.names=F, col.names=F)
write.table(EPIC.hg19[common.pd.selected.cpgs,c('seqnames','start','end','cpgname')], file='selected.common.pd.hg19.bed', quote=F, sep="\t", row.names=F, col.names=F)

write.table(annotationEPIC[hmC.selected.cpgs,c("cpgname","UCSC_RefGene_Name")], file='hmC.selected.cpgs.annot.csv', quote=F, sep="\t", row.names=F)
write.table(annotationEPIC[hmC.pd.selected.cpgs,c("cpgname","UCSC_RefGene_Name")], file='hmC.pd.selected.cpgs.annot.csv', quote=F, sep="\t", row.names=F)
write.table(annotationEPIC[bs.pd.selected.cpgs,c("cpgname","UCSC_RefGene_Name")], file='bs.pd.selected.cpgs.annot.csv', quote=F, sep="\t", row.names=F)
write.table(annotationEPIC[common.selected.cpgs,c("cpgname","UCSC_RefGene_Name")], file='common.selected.cpgs.annot.csv', quote=F, sep="\t", row.names=F)
write.table(annotationEPIC[common.pd.selected.cpgs,c("cpgname","UCSC_RefGene_Name")], file='common.pd.selected.cpgs.annot.csv', quote=F, sep="\t", row.names=F)

df <- annotationEPIC[hmC.selected.cpgs,c("cpgname","UCSC_RefGene_Name")]
s <- strsplit(df$UCSC_RefGene_Name, split = ";")
mapping <- data.frame(cpgname = rep(df$cpgname, sapply(s, length)), genename = unlist(s))
mapping.nodup <- mapping[!duplicated(mapping), ]
t1 <- table(mapping.nodup$genename)

df <- annotationEPIC[hmC.pd.selected.cpgs,c("cpgname","UCSC_RefGene_Name")]
s <- strsplit(df$UCSC_RefGene_Name, split = ";")
mapping <- data.frame(cpgname = rep(df$cpgname, sapply(s, length)), genename = unlist(s))
mapping.nodup <- mapping[!duplicated(mapping), ]
t2 <- table(mapping.nodup$genename)

df <- annotationEPIC[bs.pd.selected.cpgs,c("cpgname","UCSC_RefGene_Name")]
s <- strsplit(df$UCSC_RefGene_Name, split = ";")
mapping <- data.frame(cpgname = rep(df$cpgname, sapply(s, length)), genename = unlist(s))
mapping.nodup <- mapping[!duplicated(mapping), ]
t3 <- table(mapping.nodup$genename)

df <- annotationEPIC[,c("cpgname","UCSC_RefGene_Name")]
s <- strsplit(df$UCSC_RefGene_Name, split = ";")
mapping <- data.frame(cpgname = rep(df$cpgname, sapply(s, length)), genename = unlist(s))
mapping.nodup <- mapping[!duplicated(mapping), ]
t0 <- table(mapping.nodup$genename)

m <- merge(t1,t2, by='Var1',all=T)
colnames(m) <- c('gene', 'hmC','hmC.pd')

m <- merge(m,t3, by.x='gene', by.y='Var1', all=T)
colnames(m) <- c('gene', 'hmC','hmC.pd','bs.pd')

m <- merge(t0, m, by.x='Var1', by.y='gene')
colnames(m) <- c('gene', 'total', 'hmC','hmC.pd','bs.pd')

m[is.na(m)] <- 0

write.csv(m,file='gene.table.csv',row.names=F)

listMarts()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

# mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

getBM(attributes = c("affy_hg_u95av2", "hgnc_symbol", "chromosome_name", "band"),
      filters    = "affy_hg_u95av2",
      values     = c("1939_at","1503_at","1454_at"), 
      mart       = mart)

panther.assoc.file <- '/home/alabarga/NavarraBiomed/analysis/SequenceAssociationPathway3.6.1.txt'
panther <- read.panther(panther.assoc.file)
P00049<- subset(panther,(pathway.acc=="P00049") & (organism =="HUMAN"))

P00049.genes <- getBM(attributes = c("uniprotswissprot", "hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                      filters    = "uniprotswissprot",
                      values     = s$uniprot.id, 
                      mart       = mart)

chrs <- c(as.character(1:23), "X","Y")

P00049.genes <- subset(P00049.genes, chromosome_name %in% chrs)

my.granges <- GRanges(seqnames = paste("chr",P00049.genes$chromosome_name,sep=""), ranges = IRanges(start = P00049.genes$start_position, end = P00049.genes$end_position))
o <- findOverlaps(my.granges, EPIC.hg19.manifest)
P00049.cpgs <- names(EPIC.hg19.manifest)[unique(subjectHits(o))]

