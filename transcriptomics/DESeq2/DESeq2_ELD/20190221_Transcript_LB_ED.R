
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
library(EnhancedVolcano)
library("vsn")


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
sampleTable_ED_LB = sampleTable[1:9,]


#Reading Combined Read Counts 
counts.all = read.csv("./filt_txt_all_2/counts.all_ED_ S1583.csv", header = T) #change it 
rownames(counts.all) = counts.all$X
counts.all = counts.all[,-1]
colnames(counts.all) = gsub(".txt","",colnames(counts.all))
colnames(counts.all) = sampleTable$sampleName
counts.all_ED_LB = counts.all[,1:9]

#pheatmap(counts.all_ED_LB,
#  scale="row",
#  cluster_rows=F, cluster_cols=T,
#  show_rownames = F, show_colnames = T, 
#  fontsize=10, legend=TRUE,
#  main = "HeatMap \nTotal Read Counts - LB\nTranscriptomic"
#)



################################################################
#   PCA Analysis - LB
################################################################

#pca = prcomp(t(counts.all_ED_LB), center = TRUE, scale = FALSE)
#percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
#pca2 = cbind(as.data.frame(pca$x), sampleTable_ED_LB)   #check the order of ko and WT

##pdf("../Figures/LB/PCA_Plot_LB.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
#ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype)) + geom_point(size=3) + 
#  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
#  ylab(paste0("PC2: ",percent.var[2], "% variance")) + 
#  labs(colour = "Genotype",shape="Media") +
#  ggtitle("PCA Plot - LB \nTranscriptomic Analysis") + #CHANGE THE NAME HERE
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

cts.LB = counts.all_ED_LB
##  Note: 
#In order to benefit from the default settings of the package, 
#you should put the variable of interest at the end of the formula and make sure the control level is the first level.
cts.LB = cts.LB[,cbind("E_LB_WT1_S1  ", "E_LB_WT2_S1  ",  "E_LB_WT3_S1  ",
                       "E_LB_del1_S1  ", "E_LB_del2_S1  ", "E_LB_del3_S1  ",
                       "E_LB_S1_S1  ",   "E_LB_S2_S1  ",   "E_LB_S3_S1  ")]
cts.LB_ko = cts.LB[,1:6]
cts.LB_dc = cts.LB[,c(1:3,7:9)]
cts.LB_ko_vs_dc =  cts.LB[, c(7:9, 4:6)]


colData.LB = sampleTable_ED_LB    #Match order with column order of cts.LB
colData.LB_ko = colData.LB[c(7:9, 1:3),]
colData.LB_dc = colData.LB[c(7:9),]
colData.LB_dc = rbind(colData.LB_dc, colData.LB[c(4:6),])
colData.LB_ko_vs_dc = colData.LB[c(4:6, 1:3),] 

#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. 
rownames(colData.LB) = colData.LB$sampleName
rownames(colData.LB_ko) = colData.LB_ko$sampleName
rownames(colData.LB_dc) = colData.LB_dc$sampleName
rownames(colData.LB_ko_vs_dc) = colData.LB_ko_vs_dc$sampleName

#all(rownames(colData.LB) %in% colnames(cts.LB)) # should be TRUE
#all(rownames(colData.LB) == colnames(cts.LB)) # should be FALSE  
#cts.LB <- cts.LB[, rownames(colData.LB)]
#all(rownames(colData.LB) == colnames(cts.LB)) # should be TRUE

#for ko vs wt 
all(rownames(colData.LB_ko) %in% colnames(cts.LB_ko)) # should be TRUE
all(rownames(colData.LB_ko) == colnames(cts.LB_ko)) # should be TRUE  
cts.LB_ko <- cts.LB_ko[, rownames(colData.LB_ko)]
all(rownames(colData.LB_ko) == colnames(cts.LB_ko)) # should be TRUE, makes sure the names of the rows and the columns are the same

#for dc vs wt 
all(rownames(colData.LB_dc) %in% colnames(cts.LB_dc)) # should be TRUE
all(rownames(colData.LB_dc) == colnames(cts.LB_dc)) # should be TRUE, if FALSE, next line applies but i have doe it to all to avoid errors
cts.LB_dc <- cts.LB_dc[, rownames(colData.LB_dc)]
all(rownames(colData.LB_dc) == colnames(cts.LB_dc)) # should be TRUE, makes sure the names of the rows and the columns are the same

#for ko vs dc
all(rownames(colData.LB_ko_vs_dc) %in% colnames(cts.LB_ko_vs_dc)) # should be TRUE
all(rownames(colData.LB_ko_vs_dc) == colnames(cts.LB_ko_vs_dc)) # should be TRUE, if FALSE, next line applies but i have doe it to all to avoid errors
cts.LB_ko_vs_dc <- cts.LB_ko_vs_dc[, rownames(colData.LB_ko_vs_dc)]
all(rownames(colData.LB_ko_vs_dc) == colnames(cts.LB_ko_vs_dc)) # should be TRUE, makes sure the names of the rows and the columns are the same


library("DESeq2")
#dds.LB <- DESeqDataSetFromMatrix(countData = cts.LB,
#                                 colData = colData.LB,
#                                 design = ~ genotype)
#dds.LB

#for ko vs wt 
dds.LB_ko <- DESeqDataSetFromMatrix(countData = cts.LB_ko,
                                 colData = colData.LB_ko,
                                 design = ~ genotype)
dds.LB_ko

#for dc vs wt 
dds.LB_dc <- DESeqDataSetFromMatrix(countData = cts.LB_dc,
                                    colData = colData.LB_dc,
                                    design = ~ genotype)
dds.LB_dc


#for ko vs dc 
dds.LB_ko_vs_dc <- DESeqDataSetFromMatrix(countData = cts.LB_ko_vs_dc,
                                    colData = colData.LB_ko_vs_dc,
                                    design = ~ genotype)
dds.LB_ko_vs_dc




################################################################
#5. DESeq2 - Pre-filtering
################################################################

#keep.LB = rowSums(counts(dds.LB)) >= 20
#dds.LB = dds.LB[keep.LB,]

#for ko vs wt 
keep.LB_ko = rowSums(counts(dds.LB_ko)) >= 20
dds.LB_ko = dds.LB_ko[keep.LB_ko,]

#for dc vs wt 
keep.LB_dc = rowSums(counts(dds.LB_dc)) >= 20
dds.LB_dc = dds.LB_dc[keep.LB_dc,]

#for ko vs dc 
keep.LB_ko_vs_dc = rowSums(counts(dds.LB_ko_vs_dc)) >= 20
dds.LB_ko_vs_dc = dds.LB_ko_vs_dc[keep.LB_ko_vs_dc,]



################################################################
#  #### Define factor levels
################################################################

# Method 1 - Using Factor
#dds.LB$genotype <- factor(dds.LB$genotype, levels = c("wt","ko"))
dds.LB_ko$genotype = factor(dds.LB_ko$genotype, levels = c("wt","ko"))
dds.LB_dc$genotype = factor(dds.LB_dc$genotype, levels = c("wt","dc"))
dds.LB_ko_vs_dc$genotype = factor(dds.LB_ko_vs_dc$genotype, levels = c("dc","ko"))


### OR ###

# Method 2 - Using relevel
#dds.LB$genotype = relevel(dds.LB$genotype, ref = "wt")


################################################################
#5. DESeq2 - Differential expression analysis
################################################################

#dds.LB <- DESeq(dds.LB)      # DESeq is DONE at this step!! 
#res.LB <- results(dds.LB)
#res.LB                    ## WORKS! ko vs wt = i.e., condition treated vs untreated


#for ko vs wt 
dds.LB_ko <- DESeq(dds.LB_ko)      # DESeq is DONE at this step!! 
res.LB_ko <- results(dds.LB_ko)
res.LB_ko                    ## WORKS! ko vs wt = i.e., condition treated vs untreated

#for dc vs wt 
dds.LB_dc <- DESeq(dds.LB_dc)      # DESeq is DONE at this step!! 
res.LB_dc <- results(dds.LB_dc)
res.LB_dc                    ## WORKS! dc vs wt = i.e., condition treated vs untreated

#log2 fold change (MLE): genotype ko vs wt      ### log2 fold change (MAP): condition treated vs untreated 
#Wald test p-value: genotype ko vs wt 

#for ko vs dc 
dds.LB_ko_vs_dc <- DESeq(dds.LB_ko_vs_dc)      # DESeq is DONE at this step!! 
res.LB_ko_vs_dc <- results(dds.LB_ko_vs_dc)
res.LB_ko_vs_dc                    ## WORKS! genotype ko vs dc = i.e., condition treated vs untreated


#log2 fold change (MLE): genotype ko vs dc 
#Wald test p-value: genotype ko vs dc 
#DataFrame with 5977 rows and 6 columns



################################################################
#5. Log fold change shrinkage for visualization and ranking
################################################################

#resultsNames(dds.LB)
#resLFC.LB <- lfcShrink(dds.LB, coef="genotype_ko_vs_wt", type="apeglm") #apeglm method for effect size shrinkage
#resLFC.LB

#for ko vs wt 
resultsNames(dds.LB_ko)
resLFC.LB_ko <- lfcShrink(dds.LB_ko, coef="genotype_ko_vs_wt", type="apeglm") #apeglm method for effect size shrinkage
resLFC.LB_ko

#for dc vs wt 
resultsNames(dds.LB_dc)
resLFC.LB_dc <- lfcShrink(dds.LB_dc, coef="genotype_dc_vs_wt", type="apeglm") #apeglm method for effect size shrinkage
resLFC.LB_dc

#for ko vs dc 
resultsNames(dds.LB_ko_vs_dc)
resLFC.LB_ko_vs_dc <- lfcShrink(dds.LB_ko_vs_dc, coef="genotype_ko_vs_dc", type="apeglm") #apeglm method for effect size shrinkage
resLFC.LB_ko_vs_dc

#log2 fold change (MAP): genotype ko vs dc 
#Wald test p-value: genotype ko vs dc 
#DataFrame with 5977 rows and 5 columns


################################################################
# p-values and adjusted p-values
################################################################


#resOrdered.LB <- res.LB[order(res.LB$pvalue),]
#summary(res.LB)
#sum(res.LB$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
#res05 <- results(dds.LB, alpha=0.05) ## adjusted p-value < 0.05
#summary(res05)
#sum(res05$padj < 0.05, na.rm=TRUE)

#for ko vs wt 
resOrdered.LB_ko <- res.LB_ko[order(res.LB_ko$padj),]
summary(res.LB_ko)
sum(res.LB_ko$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
sum(res.LB_ko$padj < 0.01, na.rm=TRUE) #How many adjusted p-values were less than 0.01?
sum(res.LB_ko$padj < 0.001, na.rm=TRUE) #How many adjusted p-values were less than 0.001?
sum(res.LB_ko$padj < 0.0001, na.rm=TRUE) #How many adjusted p-values were less than 0.0001?

res05_ko <- results(dds.LB_ko, alpha=0.05) ## adjusted p-value < 0.05
summary(res05_ko)
sum(res05_ko$padj < 0.05, na.rm=TRUE)   #1246

#write.csv(resOrdered.LB_ko, "./diff.exp.gene/DEG_ED/resOrdered.LB_ko.csv")
#write.csv(res.LB_ko, "./diff.exp.gene/DEG_ED/res.LB_ko.csv")


#for dc vs wt 
resOrdered.LB_dc <- res.LB_dc[order(res.LB_dc$padj),]
summary(res.LB_dc)
sum(res.LB_dc$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
sum(res.LB_dc$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
sum(res.LB_dc$padj < 0.01, na.rm=TRUE) #How many adjusted p-values were less than 0.01?
sum(res.LB_dc$padj < 0.001, na.rm=TRUE) #How many adjusted p-values were less than 0.001?
sum(res.LB_dc$padj < 0.0001, na.rm=TRUE) #How many adjusted p-values were less than 0.0001?

res05_dc <- results(dds.LB_dc, alpha=0.05) ## adjusted p-value < 0.05
summary(res05_dc)
sum(res05_dc$padj < 0.05, na.rm=TRUE)


#for ko vs dc 
resOrdered.LB_ko_vs_dc <- res.LB_ko_vs_dc[order(res.LB_ko_vs_dc$padj),]
summary(res.LB_ko_vs_dc)
sum(res.LB_ko_vs_dc$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
sum(res.LB_ko_vs_dc$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
sum(res.LB_ko_vs_dc$padj < 0.01, na.rm=TRUE) #How many adjusted p-values were less than 0.01?
sum(res.LB_ko_vs_dc$padj < 0.001, na.rm=TRUE) #How many adjusted p-values were less than 0.001?
sum(res.LB_ko_vs_dc$padj < 0.0001, na.rm=TRUE) #How many adjusted p-values were less than 0.0001?

res05_ko_vs_dc <- results(dds.LB_ko_vs_dc, alpha=0.05) ## adjusted p-value < 0.05
summary(res05_ko_vs_dc)
sum(res05_ko_vs_dc$padj < 0.05, na.rm=TRUE)



################################################################
# Exporting results to CSV files
################################################################


#write.csv(as.data.frame(resOrdered), 
#          file="diff.exp.genes_LB.csv")         #prints all the differentially expressed genes
#resSig.LB.05 <- subset(resOrdered.LB, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#resSig.LB
#write.csv(as.data.frame(resSig.LB.05), 
#          file="diff.exp.genes_SIG.05_LB.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#resSig.LB.01 <- subset(resOrdered.LB, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#resSig.LB
#write.csv(as.data.frame(resSig.LB.01), 
#          file="diff.exp.genes_SIG.01_LB.csv")



resSig.LB_ko <- subset(resOrdered.LB_ko, padj < 1)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.LB_ko
#write.csv(as.data.frame(resSig.LB_ko), 
#          file="DEG_ko_vs_wt_p0.05_LB_eld.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.

#for ko vs wt 
resSig.LB.05_ko <- subset(resOrdered.LB_ko, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.LB.05_ko
dim(resSig.LB.05_ko)
#write(resSig.LB.05_ko, "./diff.exp.gene/DEG_ED/resSig.LB.05_ko.csv")
#write.csv(as.data.frame(resSig.LB.05_ko), 
#          file="DEG_ko_vs_wt_p0.05_LB_eld.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.

resSig.LB.01_ko <- subset(resOrdered.LB_ko, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.LB.01_ko
#write.csv(resSig.LB.01_ko, "./diff.exp.gene/DEG_ED/resSig.LB.01_ko.csv")
#write.csv(as.data.frame(resSig.LB.01_ko), 
#          file="DEG_ko_vs_wt_p0.01_LB_eld.csv")


#for dc vs wt 
resSig.LB.05_dc <- subset(resOrdered.LB_dc, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.LB.05_dc
dim(resSig.LB.05_dc)
#write.csv(resSig.LB.05_dc, "./diff.exp.gene/DEG_ED/resSig.LB.05_dc.csv")
#write.csv(as.data.frame(resSig.LB.05_dc), 
#          file="DEG_dc_vs_wt_p0.05_LB_eld.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.

resSig.LB.01_dc <- subset(resOrdered.LB_dc, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.LB.01_dc
#write.csv(resSig.LB.01_dc, "./diff.exp.gene/DEG_ED/LB/DEG_raw_LB_eld/resSig.LB.01_dc.csv")
#write.csv(as.data.frame(resSig.LB.01_dc), 
#          file="DEG_dc_vs_wt_p0.01_LB_eld.csv")


#for ko vs dc 
resSig.LB.05_ko_vs_dc <- subset(resOrdered.LB_ko_vs_dc, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.LB.05_ko_vs_dc
#log2 fold change (MLE): genotype ko vs dc 
#Wald test p-value: genotype ko vs dc 
#DataFrame with 106 rows and 6 columns
dim(resSig.LB.05_ko_vs_dc)
#write.csv(resSig.LB.05_ko_vs_dc, "./diff.exp.gene/DEG_ED/LB/DEG_raw_LB_eld/resSig.LB.05_ko_vs_dc.csv")
#write.csv(as.data.frame(resSig.LB.05_dc),        #original script
#          file="DEG_dc_vs_wt_p0.05_LB_eld.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.

resSig.LB.01_ko_vs_dc <- subset(resOrdered.LB_ko_vs_dc, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.LB.01_ko_vs_dc
#log2 fold change (MLE): genotype ko vs dc 
#Wald test p-value: genotype ko vs dc 
#DataFrame with 20 rows and 6 columns
#write.csv(resSig.LB.01_ko_vs_dc, "./diff.exp.gene/DEG_ED/LB/DEG_raw_LB_eld/resSig.LB.01_ko_vs_dc.csv")
#write.csv(as.data.frame(resSig.LB.01_dc), 
#          file="DEG_dc_vs_wt_p0.01_LB_eld.csv")





################################################################
# Exploring and exporting results - MA-plot
################################################################

##Points which fall out of the window are plotted as open triangles pointing either up or down.
#plotMA(res.LB, ylim=c(-2,2)) #Points will be colored red if the adjusted p value is less than 0.1
##It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
#plotMA(resLFC.LB, ylim=c(-2,2))
#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
#idx.LB <- identify(res.LB$baseMean, res.LB$log2FoldChange)
#rownames(res.LB)[idx.LB]  # Doesn't end


#for ko vs wt 
#plotMA(res.LB_ko, ylim=c(-2,2)) #Points will be colored red if the adjusted p value is less than 0.1
#plotMA(resLFC.LB_ko, ylim=c(-2,2))
#idx.LB_ko <- identify(res.LB_ko$baseMean, res.LB_ko$log2FoldChange)
#rownames(res.LB_ko)[idx.LB_ko]  # Doesn't end


#for dc vs wt 
#plotMA(res.LB_dc, ylim=c(-2,2)) #Points will be colored red if the adjusted p value is less than 0.1
#plotMA(resLFC.LB_dc, ylim=c(-2,2))
#idx.LB_dc <- identify(res.LB_dc$baseMean, res.LB_dc$log2FoldChange)
#rownames(res.LB_dc)[idx.LB_dc]  # Doesn't end


################################################################
# Plot counts
################################################################

#plotCounts(dds.LB, gene=which.min(res.LB$padj), intgroup="genotype")
#d.LB <- plotCounts(dds.LB, gene=which.min(res.LB$padj), intgroup="genotype", 
#                   returnData=TRUE)
#library("ggplot2")
#ggplot(d.LB, aes(x=genotype, y=count)) + 
#  geom_point(position=position_jitter(w=0.1,h=0)) + 
#  scale_y_log10(breaks=c(25,100,400))

#for ko vs wt 
plotCounts(dds.LB_ko, gene=which.min(res.LB_ko$padj), intgroup="genotype")
d.LB_ko <- plotCounts(dds.LB_ko, gene=which.min(res.LB_ko$padj), intgroup="genotype", 
                   returnData=TRUE)
library("ggplot2")
#ggplot(d.LB_ko, aes(x=genotype, y=count)) + 
#  geom_point(position=position_jitter(w=0.1,h=0)) + 
#  scale_y_log10(breaks=c(25,100,400))

#for dc vs wt 
#plotCounts(dds.LB_dc, gene=which.min(res.LB_dc$padj), intgroup="genotype")
#d.LB_dc <- plotCounts(dds.LB_dc, gene=which.min(res.LB_dc$padj), intgroup="genotype", 
                      returnData=TRUE)
library("ggplot2")
#ggplot(d.LB_dc, aes(x=genotype, y=count)) + 
#  geom_point(position=position_jitter(w=0.1,h=0)) + 
#  scale_y_log10(breaks=c(25,100,400))




################################################################
#   Heatmap of the count matrix - LB
################################################################

#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix

## this gives log2(n + 1)
ntd <- normTransform(dds.LB_ko)
library("vsn")
meanSdPlot(assay(ntd))
library("pheatmap")
#select <- order(rowMeans(counts(dds.LB_ko,normalized=TRUE)),
#                decreasing=TRUE)[1:20]
#df <- as.data.frame(colData(dds.LB_ko)[,c("genotype","media")])
#pheatmap(assay(ntd)[select,], 
#         cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=TRUE, annotation_col=df,
#         scale = "row")


#vsd <- vst(dds.LB_ko, blind=FALSE)
#head(assay(vsd), 3)
#meanSdPlot(assay(vsd))
#pheatmap(assay(vsd)[select,], 
#         cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=TRUE, annotation_col=df,
#         scale = "row")


select <- order(rowMeans(counts(dds.LB_ko,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds.LB_ko)[,c("genotype","media")])
rld <- rlog(dds.LB_ko, blind=FALSE)
meanSdPlot(assay(rld))
pheatmap(assay(rld)[select,], 
         cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df,
         scale = "row",
         color = colorRampPalette(c("white", "goldenrod")) (10))


select <- order(rowMeans(counts(dds.LB_dc,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds.LB_ko)[,c("genotype","media")])
rld <- rlog(dds.LB_ko, blind=FALSE)
meanSdPlot(assay(rld))
pheatmap(assay(rld)[select,], 
         cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df,
         scale = "row",
         color = colorRampPalette(c("white", "blue")) (10))




################################################################
# More information on results columns
################################################################

#mcols(res.LB)$description
mcols(res.LB_ko)$description
mcols(res.LB_dc)$description


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
