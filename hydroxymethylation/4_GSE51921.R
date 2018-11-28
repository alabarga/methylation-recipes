# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Feb 20 08:55:38 EST 2018
# DNA methylation analysis in idiopathic and LRRK2-associated Parkinson's disease (PD)
# https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE51921
#
################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

source("~/NavarraBiomed/analysis/parkinson/code/experiment.r")

dataDir <- "~/NavarraBiomed/analysis/parkinson/data/4_GSE51921"
workDir <- "~/NavarraBiomed/analysis/parkinson/results/4_GSE51921"
setwd(workDir)

# load series and platform data from GEO

# gset <- getGEO("GSE51921", GSEMatrix =TRUE, AnnotGPL=FALSE)

gset <- getGEO(filename=paste(dataDir,"GSE51921_series_matrix.txt.gz",sep="/"))
gset <- getGEO(filename=paste(dataDir,"GSE51921_series_matrix.txt.gz",sep="/"))

if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
sampleInfo <- pData(gset)
summary(sampleInfo)

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "0000111111111100011110001111"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")    # set group names

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# differential expression using LIMMA
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

# load NCBI platform annotation
gpl <- annotation(gset)
gpl <- "GPL13534"

platf <- getGEO(gpl, AnnotGPL = TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))
# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
tT <- merge(tT, ncbifd, by = "ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order
colnames(tT)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","CHR","MAPINFO","RANGE_START","RANGE_END","RANGE_GB","SPOT_ID","UCSC_RefGene_Name"))
write.table(tT, file="GSE51921_Top250.csv", row.names=F, sep="\t")


# differential expression using LIMMA
fl <- as.factor(sml)
gset$description <- fl
exprs(gset) <- ex_cb_corrected
  
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(ex_cb_corrected, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
sel <- rownames(tT)

colors <- c("red",'blue')[fl]
heatmap(ex_cb_corrected[sel,order(fl)], ColSideColors=colors[order(fl)], Colv=NA )

# copy the ex object created above and transpose it to get rows as samples and columns as probes
pca.data.sel <- ex_cb_corrected[sel,order(fl)]

# rename samples
grps <- labels[fl]
grpcol <- c("red",'blue')[fl]
colnames(pca.data) <- paste(grps, colnames(pca.data), sep="-")

# remove NAs
pca.data.sel <- na.omit(as.matrix(pca.data.sel))

# compute PCA
pca.sel <- prcomp(t(pca.data.sel), scale=TRUE)

color_by <- 'disease'
fl <- as.factor(sml)
grpcol <- rainbow(nlevels(fl))[fl]
labels <- levels(fl)
plot(pca.sel$x[,1], pca.sel$x[,2], xlab="PCA1", ylab="PCA2", main=paste("PCA for components 1&2 (", color_by, ")"), type="p", pch=16, col=grpcol)
legend("bottomright", labels, fill=rainbow(nlevels(fl)), bty="n")

################################################################
#   Boxplot for selected GEO samples
################################################################

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("CTRL","PD")

# set parameters and draw the plot
palette(c("#f4dfdf","#dfeaf4", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE51921", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

cell_type <- as.factor(subset(sampleInfo[colnames(ex),], select=c('characteristics_ch1.10'))[,1])
ex_cb_corrected <- sva::ComBat(ex, cell_type)

# copy the ex object created above and transpose it to get rows as samples and columns as probes
pca.data <- ex_cb_corrected

# rename samples
grps <- labels[fl]
grpcol <- c("red",'blue')[fl]
colnames(pca.data) <- paste(grps, colnames(pca.data), sep="-")

# remove NAs
pca.data <- na.omit(as.matrix(pca.data))

# compute PCA
pca <- prcomp(t(pca.data), scale=TRUE)

# identify variance in components
summary(pca)

# the first 2 component group 64% of the total variance
# the first 3 component group 79% of the total variance
# the first 4 component group 90% of the total variance

color_by <- 'disease'
fl <- as.factor(sml)


color_by <- 'cell type'
fl <- as.factor(subset(sampleInfo[colnames(ex),], select=c('characteristics_ch1.10'))[,1])

color_by <- 'sex'
fl <- as.factor(subset(sampleInfo[colnames(ex),], select=c('characteristics_ch1.2'))[,1])

grpcol <- rainbow(nlevels(fl))[fl]
labels <- levels(fl)

# components #1 and #2
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main=paste("PCA for components 1&2 (", color_by, ")"), type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,2], rownames(pca.data), cex=0.75)
legend("bottomright", labels, fill=palette(), bty="n")

subset(sampleInfo, select=c('source_name_ch1','characteristics_ch1.10','characteristics_ch1.3'))



# show other component pairs for the example
# components #1 and #3
plot(pca$x[,1], pca$x[,3], xlab="PCA1", ylab="PCA3", main=paste("PCA for components 1&3 (", color_by, ")"), type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,3], rownames(pca.data), cex=0.75)

# components #2 and #3
plot(pca$x[,2], pca$x[,3], xlab="PCA2", ylab="PCA3", main="PCA for components 2&3", type="p", pch=10, col=grpcol)
text(pca$x[,2], pca$x[,3], rownames(pca.data), cex=0.75)
