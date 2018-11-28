# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Feb 20 06:21:57 EST 2018
# Expression data from human patients with slow or rapid Parkinson's Disease progression
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80599
#
################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE80599", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13667", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "0001111011000001101101111100111100000100110001100110011011010111010"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names

################################################################
# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("Rapid","Slow")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE80599", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

# log2 transform
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
write.table(tT, file="GSE80599_Top250.csv", row.names=F, sep="\t")

sel <- rownames(tT)
colors <- c("red",'blue')[fl]
heatmap(ex[sel,order(fl)], ColSideColors=colors[order(fl)], Colv=NA )

# copy the ex object created above and transpose it to get rows as samples and columns as probes
pca.data.sel <- ex[sel,order(fl)]

# rename samples
grps <- labels[fl]
grpcol <- c("red",'blue')[fl]
colnames(pca.data) <- paste(grps, colnames(pca.data), sep="-")

# remove NAs
pca.data.sel <- na.omit(as.matrix(pca.data.sel))

# compute PCA
pca.sel <- prcomp(t(pca.data.sel), scale=TRUE)

color_by <- 'progression'
grpcol <- rainbow(nlevels(fl))[fl]
labels <- levels(fl)
plot(pca.sel$x[,1], pca.sel$x[,2], xlab="PCA1", ylab="PCA2", main=paste("PCA for components 1&2 (", color_by, ")"), type="p", pch=16, col=grpcol)
legend("bottomright", labels, fill=rainbow(nlevels(fl)), bty="n")



# load NCBI platform annotation
gpl <- annotation(gset)
gpl <- "GPL13667"

platf <- getGEO(gpl, AnnotGPL = TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))
# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
tT <- merge(tT, ncbifd, by = "ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.Symbol","GB_LIST","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


