library(qqman)

# Manhattan-plot p-value vs chr position higlighting the DMPs adjp<0.05 (FDR 5%)

require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotationEPIC <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

plotMDS <- function(x, labels=NA, col, main="MDS 2D plot of dissimilarity matrix", medoidIDs=NA){
  require(stats)
  require(ggplot2)
  
  x.scale <- scale(t(x), center = TRUE, scale = TRUE)
  
  if (any(is.na(labels))) labels <- colnames(X)
  
  #generate dissimilarity matrix
  x.diss <- daisy(x.scale, stand="FALSE")
  
  # Classical (Metric) Multidimensional Scaling
  # Also known as principal coordinates analysis (Gower, 1966)
  # 2D transposition of dissimilarity object
  mdsObj <- cmdscale(d = x.diss, k = 2, eig = TRUE) #slow
  mdsMat <- mdsObj$points
  mdsEig <- mdsObj$eig
  
  mdsDF  <- data.frame(uniqueID=labels, cluster = factor(col), 
                       var1=mdsMat[,1], var2=mdsMat[,2])
  
  if(is.na(medoidIDs)){
    plot   <- ggplot(data = mdsDF, aes(x=var1, y=var2, colour=cluster, label=labels)) + 
      geom_text() +
      ggtitle(main)
    print(plot)
  }else{
    mdsDFmed <- mdsDF[mdsDF$uniqueID %in% medoidIDs, ]
    
    plot   <- ggplot(data = mdsDF, aes(x=var1, y=var2, colour=cluster)) + 
      geom_point(position="jitter") +
      geom_point(data=mdsDFmed, size=6) +
      ggtitle(main)
    print(plot)
  }
}

volcanoPlot <- function(x, pval=0.0001, logFC=0.25, N=1000) {

  x <- x[1:N,]
  x$threshold <- (x$P.Value < pval) | (abs(x$logFC) > logFC)
  x$threshold2 <- (x$P.Value < pval) & (abs(x$logFC) > logFC)
  x$labels <- rownames(x)
  
  g = ggplot(data=x, aes(x=logFC, y=-log10(P.Value), colour = threshold)) +
    geom_point(alpha=0.4, size=1) +
    theme(legend.position = "none") +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    geom_text_repel(aes(x=logFC, y=-log10(P.Value), label = ifelse(threshold2, labels,"")))
  
 return(g)
}

plotHeat <- function(x, main="Heatmap (ggplot)") {
  
  # scale data by "columns" 
  x = data.frame(x) 
  x <- x[order(rowSums(x),decreasing = T),]
  x$names = rownames(x) 
  x.m = melt(x, id.vars = "names") 
  
  # ======================================================== 
  # Plot 
  # ======================================================== 
  ggplot(data = x.m, aes(x = variable, y = names)) + 
    geom_tile(aes(fill = value), color = "white", size = 0) + 
    scale_fill_gradient(low = "white", high = "steelblue") + 
    theme_grey(base_size = 10) + 
    ggtitle(main) + 
    xlab("Samples") + 
    theme(axis.ticks = element_blank(), 
          panel.background = element_blank(), 
          plot.title = element_text(size = 12, colour = "gray50")) 

}

# ggplot(res_tableOE_ordered) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
#   geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels == T, rownames(res_tableOE_ordered),""))) +
#   ggtitle("Mov10 overexpression") +
#   xlab("log2 fold change") + 
#   ylab("-log10 adjusted p-value") +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25))) 

#####################################
# MANHATTAN PLOT WITH QQMAN PACKAGE #
#####################################

plotManhattan <- function(limmaResult, annotationEPIC=annotationEPIC) {

# filter out bad QC and other probes from annotationEPIC -> annotationEPIC.flt

# annotationEPIC$chr is class character --> "chr1", for the manhattan() we need numeric chr1 -> 1

DMresult.plot = data.frame(SNP=annotationEPIC.flt$Name, CHR=chr, BP=annotationEPIC$pos, P=limmaResult$P.Value, stringsAsFactors=FALSE)

# I also include the deltabeta value to use the same dataset for the volcano plot, and also corrected/adjusted pvalues

#pdf("manhattan.pdf", width=8, height=6)
manhattan(DMresult.plot, suggestiveline = F, genomewideline = F, cex=0.3, col=c("azure4", "aquamarine3"))
#par(fig=c(0,1,0,1), new=TRUE)
#manhattan(subset(DMresult.plot, adjP<0.05),suggestiveline = F, genomewideline = F, cex=0.3, col=c("black", "aquamarine4")) #destacar los puntos que son significativos según la FDR que pongas
abline(h=2, col="red", lty=2, lwd=1.5) #pval=0.01
#dev.off()

}

################################
# VOLCANO PLOT WITH R GRAPHICS #
################################

plotVolcano <- function(limmaResult){
  #pdf("volcano34.pdf", width=5, height=5)
  deltabeta <- limmaResult$logFC
  P <- limmaResult$P.Value
  
  with(DMresult.plot, plot(deltabeta, -log10(P), pch=".", col="darkgray", cex=0.3, main="Volcano plot", xlim=c(-1,1))) #representar en gris todos los puntos
  with(subset(DMresult.plot, P<.01 & deltabeta>.05), points(deltabeta, -log10(P), pch=".", col="darksalmon", cex=0.3)) #a partir de aquí ir coloreando los que quieras destacar por significancia o foldchange
  with(subset(DMresult.plot, P<.01 & deltabeta>.1), points(deltabeta, -log10(P), pch=".", col="red", cex=0.3))
  with(subset(DMresult.plot, P<.01 & (-deltabeta)>.05), points(deltabeta, -log10(P), pch=".", col="darkturquoise", cex=0.3))
  with(subset(DMresult.plot, P<.01 & (-deltabeta)>.1), points(deltabeta, -log10(P), pch=".", col="blue", cex=0.3))
  abline(v=0.05, col="darkslategray", lty=2)
  abline(v=-0.05, col="darkslategray", lty=2)
  abline(v=0.1, col="darkslategray", lty=2)
  abline(v=-0.1, col="darkslategray", lty=2)
  abline(h=2, col="darkslategray", lty=2, lwd=1.5)
  #dev.off()
}
