
####################################################################
###  Analyzing RNA-seq data with DESeq2 
### Michael I. Love, Simon Anders, and Wolfgang Huber
### 01/04/2019
### cite - Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8
### #https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
####################################################################

####################################################################
### 1. Include Packages 
####################################################################


library(pheatmap)
library(ggplot2)
library(VennDiagram)
library(limma)
library(RColorBrewer)
library(UniProt.ws)
library(igraph)
library(piano)  
library(corrplot)
library(DESeq2)


####################################################################
### 2. Set Working Directory , Set Parameters 
####################################################################

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all")

thr.fc=2 #thresold_FoldChange
thr.pv =0.05 #thresold_pValue1
swati.color = c(brewer.pal(8,"Dark2")) #TA
swati.color2 = c(brewer.pal(9, "Greens")) #SM
swati.color.Blues = c(brewer.pal(8,"Blues"))


####################################################################
#3. Reading SampleTable and Combined Read Counts 
####################################################################


#Reading SampleTable
sampleTable = read.csv("./filt_txt_all_2/sampleTable_ED.csv", header = T)
sampleTable_ED_M9 = sampleTable[10:18,]



#Reading Combined Read Counts 
counts.all = read.csv("./filt_txt_all_2/counts.all_ED_ S1583.csv", header = T) #change it 
rownames(counts.all) = counts.all$X
counts.all = counts.all[,-1]
colnames(counts.all) = gsub(".txt","",colnames(counts.all))
colnames(counts.all) = sampleTable$sampleName
counts.all_ED_M9 = counts.all[,10:18]

#pheatmap(counts.all_ED_M9,
#  scale="row",
#  cluster_rows=F, cluster_cols=T,
#  show_rownames = F, show_colnames = T, 
#  fontsize=10, legend=TRUE,
#  main = "HeatMap \nTotal Read Counts - M9\nTranscriptomic"
#)




####################################################################
## PCA Analysis  - M9
####################################################################

#pca = prcomp(t(counts.all_ED_M9), center = TRUE, scale = FALSE)
#percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
#pca2 = cbind(as.data.frame(pca$x), sampleTable_ED_M9)   #check the order of ko and WT

##pdf("../Figures/M9/PCA_Plot_M9.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
#ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype)) + geom_point(size=3) + 
#  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
#  ylab(paste0("PC2: ",percent.var[2], "% variance")) + 
#  labs(colour = "Genotype",shape="Media") +
#  ggtitle("PCA Plot - M9 \nTranscriptomic Analysis") + #CHANGE THE NAME HERE
#  theme(
#    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
#    axis.title.x = element_text(color="black", size=12, face="bold"),
#    axis.title.y = element_text(color="black", size=12, face="bold")
#  ) + 
#  #geom_text(aes(PC1, PC2, colour=genotype),label=pca2$genotype) +    #adds labels to each data point
#  coord_fixed()
##dev.off()                  #turn this OFF if just want to see the picture in the Plots




################################################################
#5. DESeq2 - Count matrix input
################################################################

cts.M9 = counts.all_ED_M9
##  Note: 
#In order to benefit from the default settings of the package, 
#you should put the variable of interest at the end of the formula and make sure the control level is the first level.
cts.M9 = cts.M9[,cbind("E_M9_WT1_S2  ", "E_M9_WT2_S2  ",  "E_M9_WT3_S2  ",
                       "E_M9_del1_S2  ", "E_M9_del2_S2  ", "E_M9_del3_S2  ",
                       "E_M9_S1_S2  ",   "E_M9_S2_S2  ",   "E_M9_S3_S2  ")]
cts.M9_ko = cts.M9[,1:6]
cts.M9_dc = cts.M9[,c(1:3,7:9)]
cts.M9_ko_vs_dc = cts.M9[, c(7:9,4:6)]

colData.M9 = sampleTable_ED_M9    #Match order with column order of cts.M9
colData.M9_ko = colData.M9[c(1:3,7:9),]
colData.M9_dc = colData.M9[c(4:9),]
colData.M9_ko_vs_dc = colData.M9[c(1:6),]


#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. 
rownames(colData.M9) = colData.M9$sampleName
rownames(colData.M9_ko) = colData.M9_ko$sampleName
rownames(colData.M9_dc) = colData.M9_dc$sampleName
rownames(colData.M9_ko_vs_dc) = colData.M9_ko_vs_dc$sampleName


all(rownames(colData.M9) %in% colnames(cts.M9)) # should be TRUE
all(rownames(colData.M9) == colnames(cts.M9)) # should be FALSE  
cts.M9 <- cts.M9[, rownames(colData.M9)]
all(rownames(colData.M9) == colnames(cts.M9)) # should be TRUE

#for ko vs wt 
all(rownames(colData.M9_ko) %in% colnames(cts.M9_ko)) # should be TRUE
all(rownames(colData.M9_ko) == colnames(cts.M9_ko)) # should be TRUE 
cts.M9_ko <- cts.M9_ko[, rownames(colData.M9_ko)]
all(rownames(colData.M9_ko) == colnames(cts.M9_ko)) # should be TRUE


#for dc vs wt 
all(rownames(colData.M9_dc) %in% colnames(cts.M9_dc)) # should be TRUE
all(rownames(colData.M9_dc) == colnames(cts.M9_dc)) # should be TRUE  
cts.M9_dc <- cts.M9_dc[, rownames(colData.M9_dc)]
all(rownames(colData.M9_dc) == colnames(cts.M9_dc)) # should be TRUE


#for ko vs dc
all(rownames(colData.M9_ko_vs_dc) %in% colnames(cts.M9_ko_vs_dc)) # should be TRUE
all(rownames(colData.M9_ko_vs_dc) == colnames(cts.M9_ko_vs_dc)) # should be TRUE, if FALSE, next line applies but i have done it to all to avoid errors
cts.M9_ko_vs_dc <- cts.M9_ko_vs_dc[, rownames(colData.M9_ko_vs_dc)]
all(rownames(colData.M9_ko_vs_dc) == colnames(cts.M9_ko_vs_dc)) # should be TRUE, makes sure the names of the rows and the columns are the same


library("DESeq2")
dds.M9 <- DESeqDataSetFromMatrix(countData = cts.M9,
                                 colData = colData.M9,
                                 design = ~ genotype)
dds.M9

#for ko vs wt 
dds.M9_ko <- DESeqDataSetFromMatrix(countData = cts.M9_ko,
                                 colData = colData.M9_ko,
                                 design = ~ genotype)
dds.M9_ko



#for dc vs wt 
dds.M9_dc <- DESeqDataSetFromMatrix(countData = cts.M9_dc,
                                    colData = colData.M9_dc,
                                    design = ~ genotype)
dds.M9_dc


#for ko vs dc 
dds.M9_ko_vs_dc <- DESeqDataSetFromMatrix(countData = cts.M9_ko_vs_dc,
                                    colData = colData.M9_ko_vs_dc,
                                    design = ~ genotype)
dds.M9_ko_vs_dc


################################################################
#5. DESeq2 - Pre-filtering
################################################################

keep.M9 = rowSums(counts(dds.M9)) >= 20
dds.M9 = dds.M9[keep.M9,]

#for ko vs wt 
keep.M9_ko = rowSums(counts(dds.M9_ko)) >= 20
dds.M9_ko = dds.M9_ko[keep.M9_ko,]

#for dc vs wt 
keep.M9_dc = rowSums(counts(dds.M9_dc)) >= 20
dds.M9_dc = dds.M9_dc[keep.M9_dc,]

#for ko vs dc 
keep.M9_ko_vs_dc = rowSums(counts(dds.M9_ko_vs_dc)) >= 20
dds.M9_ko_vs_dc = dds.M9_ko_vs_dc[keep.M9_ko_vs_dc,]


#### Define factor levels

# Method 1 - Using Factor
#dds.M9$genotype <- factor(dds.M9$genotype, levels = c("wt","ko"))
dds.M9_ko$genotype = factor(dds.M9_ko$genotype, levels = c("wt","ko"))
dds.M9_dc$genotype = factor(dds.M9_dc$genotype, levels = c("wt","dc"))
dds.M9_ko_vs_dc$genotype = factor(dds.M9_ko_vs_dc$genotype, levels = c("dc","ko"))



### OR ###

# Method 2 - Using relevel
#dds.M9$genotype = relevel(dds.M9$genotype, ref = "wt")


################################################################
#5. DESeq2 - Differential expression analysis
################################################################

dds.M9 <- DESeq(dds.M9)      # DESeq is DONE at this step!! 
res.M9 <- results(dds.M9)
res.M9                    ## WORKS! ko vs wt = i.e., condition treated vs untreated


#for ko vs wt 
dds.M9_ko <- DESeq(dds.M9_ko)      # DESeq is DONE at this step!! 
res.M9_ko <- results(dds.M9_ko)
res.M9_ko                    ## WORKS! ko vs wt = i.e., condition treated vs untreated
#log2 fold change (MLE): genotype ko vs wt 
#Wald test p-value: genotype ko vs wt 
#DataFrame with 6100 rows and 6 columns


#for dc vs wt 
dds.M9_dc <- DESeq(dds.M9_dc)      # DESeq is DONE at this step!! 
res.M9_dc <- results(dds.M9_dc)
res.M9_dc                    ## WORKS! dc vs wt = i.e., condition treated vs untreated

#log2 fold change (MLE): genotype ko vs wt 
#Wald test p-value: genotype ko vs wt 
#log2 fold change (MLE): genotype dc vs wt 
#Wald test p-value: genotype dc vs wt 
#DataFrame with 6143 rows and 6 columns

#for ko vs dc 
dds.M9_ko_vs_dc <- DESeq(dds.M9_ko_vs_dc)      # DESeq is DONE at this step!! 
res.M9_ko_vs_dc <- results(dds.M9_ko_vs_dc)
res.M9_ko_vs_dc                    ## WORKS! dc vs wt = i.e., condition treated vs untreated
#log2 fold change (MLE): genotype ko vs dc 
#Wald test p-value: genotype ko vs dc 
#DataFrame with 6111 rows and 6 columns




################################################################
#5. Log fold change shrinkage for visualization and ranking
################################################################

#resultsNames(dds.M9)
#resLFC.M9 <- lfcShrink(dds.M9, coef="genotype_ko_vs_wt", type="apeglm") #apeglm method for effect size shrinkage
#resLFC.M9

#for ko vs wt 
resultsNames(dds.M9_ko)
resLFC.M9_ko <- lfcShrink(dds.M9_ko, coef="genotype_ko_vs_wt", type="apeglm") #apeglm method for effect size shrinkage
resLFC.M9_ko

#for dc vs wt 
resultsNames(dds.M9_dc)
resLFC.M9_dc <- lfcShrink(dds.M9_dc, coef="genotype_dc_vs_wt", type="apeglm") #apeglm method for effect size shrinkage
resLFC.M9_dc

#for ko vs dc 
resultsNames(dds.M9_ko_vs_dc)
resLFC.M9_ko_vs_dc <- lfcShrink(dds.M9_ko_vs_dc, coef="genotype_ko_vs_dc", type="apeglm") #apeglm method for effect size shrinkage
resLFC.M9_ko_vs_dc
#log2 fold change (MAP): genotype ko vs dc 
#Wald test p-value: genotype ko vs dc 
#DataFrame with 6111 rows and 5 columns




################################################################
# p-values and adjusted p-values
################################################################


#resOrdered.M9 <- res.M9[order(res.M9$pvalue),]
#summary(res.M9)
#sum(res.M9$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
#res05 <- results(dds.M9, alpha=0.05) ## adjusted p-value < 0.05
#summary(res05)
#sum(res05$padj < 0.05, na.rm=TRUE)

#for ko vs wt 
resOrdered.M9_ko <- res.M9_ko[order(res.M9_ko$padj),]
summary(res.M9_ko)
sum(res.M9_ko$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
sum(res.M9_ko$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
sum(res.M9_ko$padj < 0.01, na.rm=TRUE) #How many adjusted p-values were less than 0.01?
sum(res.M9_ko$padj < 0.001, na.rm=TRUE) #How many adjusted p-values were less than 0.001?
sum(res.M9_ko$padj < 0.0001, na.rm=TRUE) #How many adjusted p-values were less than 0.0001?

res05_ko <- results(dds.M9_ko, alpha=0.05) ## adjusted p-value < 0.05
summary(res05_ko)
sum(res05_ko$padj < 0.05, na.rm=TRUE)

#for dc vs wt 
resOrdered.M9_dc <- res.M9_dc[order(res.M9_dc$padj),]
summary(res.M9_dc)
sum(res.M9_dc$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
sum(res.M9_ko$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
sum(res.M9_ko$padj < 0.01, na.rm=TRUE) #How many adjusted p-values were less than 0.01?
sum(res.M9_ko$padj < 0.001, na.rm=TRUE) #How many adjusted p-values were less than 0.001?
sum(res.M9_ko$padj < 0.0001, na.rm=TRUE) #How many adjusted p-values were less than 0.0001?

res05_dc <- results(dds.M9_dc, alpha=0.05) ## adjusted p-value < 0.05
summary(res05_dc)
sum(res05_dc$padj < 0.05, na.rm=TRUE)
sum(res05_dc$padj < 0.01, na.rm=TRUE)
sum(res05_dc$padj < 0.001, na.rm=TRUE)
sum(res05_dc$padj < 0.0001, na.rm=TRUE)

#for ko vs dc 
resOrdered.M9_ko_vs_dc <- res.M9_ko_vs_dc[order(res.M9_ko_vs_dc$padj),]
summary(res.M9_ko_vs_dc)
sum(res.M9_ko_vs_dc$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
sum(res.M9_ko_vs_dc$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
sum(res.M9_ko_vs_dc$padj < 0.01, na.rm=TRUE) #How many adjusted p-values were less than 0.01?
sum(res.M9_ko_vs_dc$padj < 0.001, na.rm=TRUE) #How many adjusted p-values were less than 0.001?
sum(res.M9_ko_vs_dc$padj < 0.0001, na.rm=TRUE) #How many adjusted p-values were less than 0.0001?

res05_ko_vs_dc <- results(dds.M9_ko_vs_dc, alpha=0.05) ## adjusted p-value < 0.05
summary(res05_ko_vs_dc)
sum(res05_ko_vs_dc$padj < 0.05, na.rm=TRUE)




################################################################
# Exporting results to CSV files
################################################################


#write.csv(as.data.frame(resOrdered), 
#          file="diff.exp.genes_M9.csv")         #prints all the differentially expressed genes
#resSig.M9.05 <- subset(resOrdered.M9, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#resSig.M9
#write.csv(as.data.frame(resSig.M9.05), 
#          file="diff.exp.genes_SIG.05_M9.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#resSig.M9.01 <- subset(resOrdered.M9, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#resSig.M9
#write.csv(as.data.frame(resSig.M9.01), 
#          file="diff.exp.genes_SIG.01_M9.csv")


#for ko vs wt 
resSig.M9.05_ko <- subset(resOrdered.M9_ko, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.05_ko
#write.csv(resSig.M9.05_ko, "./diff.exp.gene/DEG_ED/resSig.M9.05_ko.csv")
##write.csv(as.data.frame(resSig.M9.05_ko), 
##          file="DEG_ko_vs_wt_p0.05_M9_eld.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.01_ko <- subset(resOrdered.M9_ko, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.01_ko
#write.csv(resSig.M9.01_ko, "./diff.exp.gene/DEG_ED/resSig.M9.01_ko.csv")
##write.csv(as.data.frame(resSig.M9.01_ko), 
##          file="DEG_ko_vs_wt_p0.01_M9_eld.csv")



#for dc vs wt 
resSig.M9.05_dc <- subset(resOrdered.M9_dc, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.05_dc
#write.csv(resSig.M9.05_dc, "./diff.exp.gene/DEG_ED/resSig.M9.05_dc.csv")
##write.csv(as.data.frame(resSig.M9.05_dc), 
##          file="DEG_dc_vs_wt_p0.05_M9_eld.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.01_dc <- subset(resOrdered.M9_dc, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.01_dc
#write.csv(resSig.M9.01_dc, "./diff.exp.gene/DEG_ED/resSig.M9.01_dc.csv")
##write.csv(as.data.frame(resSig.M9.01_dc), 
##          file="DEG_dc_vs_wt_p0.01_M9_eld.csv")


#for ko vs dc 
resSig.M9.05_ko_vs_dc <- subset(resOrdered.M9_ko_vs_wt, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.05_ko_vs_dc
#write.csv(resSig.M9.05_ko_vs_dc, "./diff.exp.gene/DEG_ED/M9/DEG_raw_M9_eld/resSig.M9.05_ko_vs_dc.csv")
##write.csv(as.data.frame(resSig.M9.05_ko_vs_dc), 
##          file="DEG_ko_vs_dc_vs_wt_p0.05_M9_eld.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.01_ko_vs_dc <- subset(resOrdered.M9_ko_vs_dc, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.01_ko_vs_dc
#write.csv(resSig.M9.01_ko_vs_dc, "./diff.exp.gene/DEG_ED/M9/DEG_raw_M9_eld/resSig.M9.01_ko_vs_dc.csv")
##write.csv(as.data.frame(resSig.M9.01_ko_vs_dc), 
##          file="DEG_ko_vs_dc_vs_wt_p0.01_M9_eld.csv")




#################

# this gives log2(n + 1)
ntd <- normTransform(dds.M9_ko)
library("vsn")
meanSdPlot(assay(ntd))

#Heatmap of the count matrix
library("pheatmap")
select <- order(rowMeans(counts(dds.M9_ko,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds.M9_ko)[,c("media","genotype")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



#Heatmap of the sample-to-sample distances

vsd <- vst(dds.M9_ko, blind=FALSE)
rld <- rlog(dds.M9_ko, blind=FALSE)
head(assay(vsd), 3)


sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$media, vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#Principal component plot of the samples

plotPCA(vsd, intgroup=c("media", "genotype"))


pcaData <- plotPCA(vsd, intgroup=c("media", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype, shape=media)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


#Approach to count outliers


par(mar=c(8,5,2,2))
boxplot(log10(assays(dds.M9_ko)[["cooks"]]), range=0, las=2)

#Dispersion plot and fitting alternatives
plotDispEsts(dds.M9_ko)


#Independent filtering of results
plot(metadata(res.M9_ko)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res.M9_ko)$lo.fit, col="red")
abline(v=metadata(res.M9_ko)$filterTheta)



#Count outlier detection
W <- res.M9_ko$stat
maxCooks <- apply(assays(dds.M9_ko)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds.M9_ko)
p <- 3
abline(h=qf(.99, p, m - p))



#Filtering criteria
plot(res.M9_ko$baseMean+1, -log10(res05_ko$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))



#Histogram of p values for all tests

use <- res05_ko$baseMean > metadata(res05_ko)$filterThreshold
h1 <- hist(res05_ko$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res05_ko$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


################################################################
# Exploring and exporting results - MA-plot
################################################################

##Points which fall out of the window are plotted as open triangles pointing either up or down.
#plotMA(res.M9, ylim=c(-2,2)) #Points will be colored red if the adjusted p value is less than 0.1
##It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
#plotMA(resLFC.M9, ylim=c(-2,2))
#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
#idx.M9 <- identify(res.M9$baseMean, res.M9$log2FoldChange)
#rownames(res.M9)[idx.M9]  # Doesn't end


#for ko vs wt 
#plotMA(res.M9_ko, ylim=c(-2,2)) #Points will be colored red if the adjusted p value is less than 0.1
#plotMA(resLFC.M9_ko, ylim=c(-2,2))
#idx.M9_ko <- identify(res.M9_ko$baseMean, res.M9_ko$log2FoldChange)
#rownames(res.M9_ko)[idx.M9_ko]  # Doesn't end


#for dc vs wt 
#plotMA(res.M9_dc, ylim=c(-2,2)) #Points will be colored red if the adjusted p value is less than 0.1
#plotMA(resLFC.M9_dc, ylim=c(-2,2))
##idx.M9_dc <- identify(res.M9_dc$baseMean, res.M9_dc$log2FoldChange)
##rownames(res.M9_dc)[idx.M9_dc]  # Doesn't end


################################################################
# Plot counts
################################################################

#plotCounts(dds.M9, gene=which.min(res.M9$padj), intgroup="genotype")
#d.M9 <- plotCounts(dds.M9, gene=which.min(res.M9$padj), intgroup="genotype", 
#                   returnData=TRUE)
#library("ggplot2")
#ggplot(d.M9, aes(x=genotype, y=count)) + 
#  geom_point(position=position_jitter(w=0.1,h=0)) + 
#  scale_y_log10(breaks=c(25,100,400))

#for ko vs wt 
#plotCounts(dds.M9_ko, gene=which.min(res.M9_ko$padj), intgroup="genotype")
#d.M9_ko <- plotCounts(dds.M9_ko, gene=which.min(res.M9_ko$padj), intgroup="genotype", 
#                   returnData=TRUE)
#library("ggplot2")
#ggplot(d.M9_ko, aes(x=genotype, y=count)) + 
#  geom_point(position=position_jitter(w=0.1,h=0)) + 
#  scale_y_log10(breaks=c(25,100,400))

##for dc vs wt 
#plotCounts(dds.M9_dc, gene=which.min(res.M9_dc$padj), intgroup="genotype")
#d.M9_dc <- plotCounts(dds.M9_dc, gene=which.min(res.M9_dc$padj), intgroup="genotype", 
#                      returnData=TRUE)
#library("ggplot2")
#ggplot(d.M9_dc, aes(x=genotype, y=count)) + 
#  geom_point(position=position_jitter(w=0.1,h=0)) + 
#  scale_y_log10(breaks=c(25,100,400))


################################################################
# More information on results columns
################################################################

#mcols(res.M9)$description
mcols(res.M9_ko)$description
mcols(res.M9_dc)$description


#Rich visualization and reporting of results
#ReportingTools
#http://bioconductor.org/packages/release/bioc/html/ReportingTools.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ReportingTools", version = "3.8")


################################################################
# #EnhancedVolcano
################################################################


#https://www.biostars.org/p/335751/
#last update: January 26, 2019
#devtools::install_github('kevinblighe/EnhancedVolcano')
#BiocManager::install('EnhancedVolcano')
#R 3.5.0 (2018-04-23)
#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')
#BiocManager::install('EnhancedVolcano')
#library(EnhancedVolcano)

#Plot the most basic volcano plot

#https://github.com/kevinblighe/EnhancedVolcano#plot-the-most-basic-volcano-plot
#EnhancedVolcano(res.SH5,
#                lab = rownames(res.SH5),
#                x = 'log2FoldChange',
#                y = 'padj',
#                xlim = c(-7, 8), ylim = c(0,26),
#                title = 'del_spoVG vs WT - SH5 \n Transcriptomics - padj < 0.05',
#                pCutoff = (res.SH5$padj)<0.01,
#                #FCcutoff = 1.5,
#                transcriptPointSize = 2.5,
#                transcriptLabSize = 3.0)




#################################################################
###                   Volcano Plots
##################################################################


#par(mar=c(5,5,5,5), #bottom, left, top and right margins
#    cex=1.0, cex.main=1.4, 
#    cex.axis=1.4, cex.lab=1.4)  
#topT <- as.data.frame(res.LB_ko)
##Adjusted P values (FDR Q values)
#with(topT, plot(log2FoldChange, -log10(padj), 
#                pch=20, main="Volcano plot - LB_ko_eld ", cex=1.0, 
#                xlab=bquote(~Log[2]~fold~change), 
#                ylab=bquote(~-log[10]~p.adj.Val)))
#with(subset(topT, padj<0.05 & abs(padj)<0.05), 
#     points(log2FoldChange, -log10(padj), 
#            pch=20, col="red", cex=0.5))
##with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))
##Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
#abline(v=-2, col="black", lty=4, lwd=2.0)
##abline(v=0, col="black", lty=3, lwd=1.0)
##abline(v=2, col="black", lty=4, lwd=2.0)
#abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)







###################################  ALL GOOD SO FAR  ################################### 

