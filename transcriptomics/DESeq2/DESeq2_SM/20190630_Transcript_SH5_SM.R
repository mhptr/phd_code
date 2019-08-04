
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

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/filt_txt_all_2/")

thr.fc=2 #thresold_FoldChange
thr.pv =0.05 #thresold_pValue1
swati.color = c(brewer.pal(8,"Dark2")) #TA
swati.color2 = c(brewer.pal(9, "Greens")) #SM
swati.color.Blues = c(brewer.pal(8,"Blues"))


####################################################################
#3. Reading SampleTable and Combined Read Counts 
####################################################################


#Reading SampleTable
sampleTable = read.csv("./sampleTable_SM.csv", header = T)
sampleTable_SM_SH5 = sampleTable[13:18,]

#Reading Combined Read Counts 
counts.all = read.csv("./counts.all_SM_S1583.csv", header = T) #change it 
rownames(counts.all) = counts.all$X
counts.all = counts.all[,-1]
colnames(counts.all) = gsub(".txt","",colnames(counts.all))
colnames(counts.all) = sampleTable$sampleName
counts.all_SM_SH5 = counts.all[,13:18]

pheatmap(counts.all_SM_SH5,
  scale="row",
  cluster_rows=F, cluster_cols=T,
  show_rownames = F, show_colnames = T, 
  fontsize=10, legend=TRUE,
  main = "HeatMap \nTotal Read Counts - SH5\nTranscriptomic"
)


####################################################################
###4.PCA Analysis  
####################################################################

pca = prcomp(t(counts.all_SM_SH5), center = TRUE, scale = FALSE)
percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
pca2 = cbind(as.data.frame(pca$x), sampleTable_SM_SH5)   #check the order of ko and WT

#pdf("../Figures/LB/PCA_Plot_LB.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
  ylab(paste0("PC2: ",percent.var[2], "% variance")) + 
  labs(colour = "Genotype",shape="Media") +
  ggtitle("PCA Plot - SH5 \nTranscriptomic Analysis") + #CHANGE THE NAME HERE
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  ) + 
  #geom_text(aes(PC1, PC2, colour=genotype),label=pca2$genotype) +    #adds labels to each data point
  coord_fixed()
#dev.off()                  #turn this OFF if just want to see the picture in the Plots




################################################################
#5. DESeq2 - Count matrix input
################################################################

cts.SH5 = counts.all_SM_SH5
##  Note: 
#In order to benefit from the default settings of the package, 
#you should put the variable of interest at the end of the formula and make sure the control level is the first level.
cts.SH5 = cts.SH5[,cbind( "S_SH5_W1", "S_SH5_W2", "S_SH5_W3", "S_SH5_DEL1", "S_SH5_DEL2", "S_SH5_DEL3")]
colData.SH5 = sampleTable_SM_SH5    #Match order with column order of cts.SH5
colData.SH5 = colData.SH5[c(4:6,1:3),]

#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. 
rownames(colData.SH5) = colData.SH5$sampleName

all(rownames(colData.SH5) %in% colnames(cts.SH5)) # should be TRUE
all(rownames(colData.SH5) == colnames(cts.SH5)) # should be TRUE  
cts.SH5 <- cts.SH5[, rownames(colData.SH5)]
all(rownames(colData.SH5) == colnames(cts.SH5)) # should be TRUE

library("DESeq2")
dds.SH5 <- DESeqDataSetFromMatrix(countData = cts.SH5,
                                 colData = colData.SH5,
                                 design = ~ genotype)
dds.SH5


#class: DESeqDataSet 
#dim: 6495 6 
#metadata(1): version
#assays(1): counts
#rownames(6495): new_4215473_4215670 new_1_148 ... new_4215437_4215505_c S1583
#rowData names(0):
#  colnames(6): S_SH5_W1 S_SH5_W2 ... S_SH5_DEL2 S_SH5_DEL3
#colData names(4): sampleName fileName genotype media




################################################################
#5. DESeq2 - Pre-filtering
################################################################

keep.SH5 = rowSums(counts(dds.SH5)) >= 20
dds.SH5 = dds.SH5[keep.SH5,]


#### Define factor levels

# Method 1 - Using Factor
#dds.SH5$genotype <- factor(dds.SH5$genotype, levels = c("wt","ko"))

### OR ###

# Method 2 - Using relevel
dds.SH5$genotype = relevel(dds.SH5$genotype, ref = "wt")


################################################################
#5. DESeq2 - Differential expression analysis
################################################################

dds.SH5 <- DESeq(dds.SH5)      # DESeq is DONE at this step!! 
res.SH5 <- results(dds.SH5)
res.SH5                    ## WORKS! ko vs wt = i.e., condition treated vs untreated

#log2 fold change (MLE): genotype ko vs wt 
#Wald test p-value: genotype ko vs wt 
#DataFrame with 6186 rows and 6 columns
#baseMean      log2FoldChange             lfcSE
#<numeric>           <numeric>         <numeric>
#  new_4215473_4215670   12.0954800570096   0.543630055515194  0.71277799677296
#new_1_148              18.504469936835   0.707404164194905 0.617120423794394
#new_24_297_c          25.4710628412446  -0.341173907460714  0.60034176336507
#new_150_409           29.9529454273955    1.20058873950729 0.523850341855708
#new_299_1067_c        245.301527417555  -0.399813651416031  0.27506710533755



################################################################
# Log fold change shrinkage for visualization and ranking
################################################################

resultsNames(dds.SH5)
resLFC.SH5 <- lfcShrink(dds.SH5, coef="genotype_ko_vs_wt", type="apeglm") #apeglm method for effect size shrinkage
resLFC.SH5

#log2 fold change (MAP): genotype ko vs wt 
#Wald test p-value: genotype ko vs wt 
#DataFrame with 6186 rows and 5 columns
#baseMean      log2FoldChange             lfcSE
#<numeric>           <numeric>         <numeric>
#  new_4215473_4215670   12.0954800570096  0.0569205949766604  0.23627827173612
#new_1_148              18.504469936835  0.0998625524490488 0.249880650990958
#new_24_297_c          25.4710628412446 -0.0480900012843495 0.228542104128988
#new_150_409           29.9529454273955   0.400004819723353 0.651463510830655
#new_299_1067_c        245.301527417555   -0.20080307885119 0.226935055724352


################################################################
# p-values and adjusted p-values
################################################################


resOrdered.SH5 <- res.SH5[order(res.SH5$padj),]
summary(res.SH5)

#out of 6186 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 196, 3.2%
#LFC < 0 (down)     : 136, 2.2%
#outliers [1]       : 0, 0%
#low counts [2]     : 1320, 21%
#(mean count < 39)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results


sum(res.SH5$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
#[1] 332


res05 <- results(dds.SH5, alpha=0.05) ## adjusted p-value < 0.05
summary(res05)

#out of 6186 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 132, 2.1%
#LFC < 0 (down)     : 96, 1.6%
#outliers [1]       : 0, 0%
#low counts [2]     : 840, 14%
#(mean count < 21)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

sum(res05$padj < 0.05, na.rm=TRUE)
#[1] 228
sum(res05$padj < 0.01, na.rm=TRUE)
#[1] 121


################################################################
# Exploring and exporting results - MA-plot
################################################################

##Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res.SH5, ylim=c(-2,2)) #Points will be colored red if the adjusted p value is less than 0.1
plotMA(res.SH5, ylim=c(-4,4))
##It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
plotMA(resLFC.SH5, ylim=c(-2,2))
plotMA(resLFC.SH5, ylim=c(-4,4))
#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
#idx.SH5 <- identify(res.SH5$baseMean, res.SH5$log2FoldChange)
#rownames(res.SH5)[idx.SH5]  # Doesn't end



################################################################
# Plot counts
################################################################

plotCounts(dds.SH5, gene=which.min(res.SH5$padj), intgroup="genotype")
d.SH5 <- plotCounts(dds.SH5, gene=which.min(res.SH5$padj), intgroup="genotype", 
                   returnData=TRUE)
library("ggplot2")
ggplot(d.SH5, aes(x=genotype, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))



################################################################
# More information on results columns
################################################################

mcols(res.SH5)$description



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
library(EnhancedVolcano)

#Plot the most basic volcano plot

#https://github.com/kevinblighe/EnhancedVolcano#plot-the-most-basic-volcano-plot
EnhancedVolcano(res.SH5,
                lab = rownames(res.SH5),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-7, 8), ylim = c(0,26),
                title = 'del_spoVG vs WT - SH5 \n Transcriptomics - padj < 0.05',
                pCutoff = (res.SH5$padj)<0.01,
                #FCcutoff = 1.5,
                transcriptPointSize = 2.5,
                transcriptLabSize = 3.0)


################################################################
# Exporting results to CSV files
################################################################


##write.csv(as.data.frame(resOrdered), 
##          file="diff.exp.genes_SH5.csv")         #prints all the differentially expressed genes
#resSig.SH5


resSig.SH5 = subset(resOrdered.SH5)
#write.csv(resSig.SH5, "../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_vs_wt_all_SH5_sm_30062019.csv")


resSig.SH5.05 <- subset(resOrdered.SH5, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#write.csv(resSig.SH5.05, "../diff.exp.gene/DEG_SM/Diff_SIG_sm/DEG_ko_vs_wt_p0.05_SH5_sm.csv")
#write.csv(as.data.frame(resSig.SH5.05), 
#          file="DEG_ko_vs_wt_p0.05_SH5_sm.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.

resSig.SH5.01 <- subset(resOrdered.SH5, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.SH5.01
#write.csv(resSig.SH5.01, "../diff.exp.gene/DEG_SM/Diff_SIG_sm/DEG_ko_vs_wt_p0.01_SH5_sm.csv")
#write.csv(as.data.frame(resSig.SH5.01), 
#          file="DEG_ko_vs_wt_p0.01_SH5_sm.csv")



###################################  ALL GOOD SO FAR  ################################### 


####################################################################
###4.PCA Analysis  
####################################################################

pca = prcomp(t(counts.all_SM_SH5), center = TRUE, scale = FALSE)
percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
pca2 = cbind(as.data.frame(pca$x), sampleTable_SM_SH5)   #check the order of ko and WT

#pdf("../Figures/LB/PCA_Plot_LB.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
  ylab(paste0("PC2: ",percent.var[2], "% variance")) + 
  labs(colour = "Genotype",shape="Media") +
  ggtitle("PCA Plot - SH5 \nTranscriptomic Analysis") + #CHANGE THE NAME HERE
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  ) + 
  #geom_text(aes(PC1, PC2, colour=genotype),label=pca2$genotype) +    #adds labels to each data point
  coord_fixed()
#dev.off()                  #turn this OFF if just want to see the picture in the Plots


