
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

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/filt_txt_all_2")

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
sampleTable_SM_SH2 = sampleTable[7:12,]

#Reading Combined Read Counts 
counts.all = read.csv("./counts.all_SM_S1583.csv", header = T) #change it 
rownames(counts.all) = counts.all$X
counts.all = counts.all[,-1]
colnames(counts.all) = gsub(".txt","",colnames(counts.all))
colnames(counts.all) = sampleTable$sampleName
counts.all_SM_SH2 = counts.all[,7:12]


################################################################
#5. DESeq2 - Count matrix input
################################################################

cts.SH2 = counts.all_SM_SH2
##  Note: 
#In order to benefit from the default settings of the package, 
#you should put the variable of interest at the end of the formula and make sure the control level is the first level.
cts.SH2 = cts.SH2[,cbind( "S_SH2_W1", "S_SH2_W2", "S_SH2_W3",
                        "S_SH2_DEL1", "S_SH2_DEL2", "S_SH2_DEL3")]
colData.SH2 = sampleTable_SM_SH2    #Match order with column order of cts.SH2


#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. 
rownames(colData.SH2) = colData.SH2$sampleName

all(rownames(colData.SH2) %in% colnames(cts.SH2)) # should be TRUE
all(rownames(colData.SH2) == colnames(cts.SH2)) # should be FALSE  
cts.SH2 <- cts.SH2[, rownames(colData.SH2)]
all(rownames(colData.SH2) == colnames(cts.SH2)) # should be TRUE

library("DESeq2")
dds.SH2 <- DESeqDataSetFromMatrix(countData = cts.SH2,
                                 colData = colData.SH2,
                                 design = ~ genotype)
dds.SH2

#output

#class: DESeqDataSet 
#dim: 6495 6 
#metadata(1): version
#assays(1): counts
#rownames(6495): new_4215473_4215670 new_1_148 ... new_4215437_4215505_c S1583
#rowData names(0):
#  colnames(6): S_SH2_DEL1 S_SH2_DEL2 ... S_SH2_W2 S_SH2_W3
#colData names(4): sampleName fileName genotype media


################################################################
#5. DESeq2 - Pre-filtering
################################################################

keep.SH2 = rowSums(counts(dds.SH2)) >= 20
dds.SH2 = dds.SH2[keep.SH2,]


#### Define factor levels

# Method 1 - Using Factor
#dds.SH2$genotype <- factor(dds.SH2$genotype, levels = c("wt","ko"))

### OR ###

# Method 2 - Using relevel
dds.SH2$genotype = relevel(dds.SH2$genotype, ref = "wt")


################################################################
#5. DESeq2 - Differential expression analysis
################################################################

dds.SH2 <- DESeq(dds.SH2)      # DESeq is DONE at this step!! 
res.SH2 <- results(dds.SH2)
res.SH2                    ## WORKS! ko vs wt = i.e., condition treated vs untreated

#output
#log2 fold change (MLE): genotype ko vs wt 
#Wald test p-value: genotype ko vs wt 
#DataFrame with 5778 rows and 6 columns
#baseMean     log2FoldChange             lfcSE               stat             pvalue
#<numeric>          <numeric>         <numeric>          <numeric>          <numeric>
#  new_24_297_c          15.1799393964703  0.123227133804314 0.686916352433933  0.179391760536325  0.857630100586679
#new_150_409           15.6887670582187  0.636274297990861 0.635249857579886   1.00161265744297  0.316530705367445
#new_299_1067_c        65.3776259377723 -0.089949202262934 0.403329604291476 -0.223016612978228  0.823522576809578
#BSU00010               254.62410051172  0.311434879960896 0.336901827991898  0.924408400563461  0.355273722613138
#new_1751_1938         43.7790438759012  0.158259698539004 0.426719133560259  0.370875562149255  0.710730217728052



################################################################
# Log fold change shrinkage for visualization and ranking
################################################################

resultsNames(dds.SH2)
#[1] "Intercept"         "genotype_ko_vs_wt"
resLFC.SH2 <- lfcShrink(dds.SH2, coef="genotype_ko_vs_wt", type="apeglm") #apeglm method for effect size shrinkage
#using 'apeglm' for LFC shrinkage. If used in published research, please cite:
#  Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
#sequence count data: removing the noise and preserving large differences.
#Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
resLFC.SH2
#View(as.data.frame.matrix(resLFC.SH2))
#write.csv(as.data.frame.matrix(resLFC.SH2), "../diff.exp.gene/DEG_SM/Diff_SIG_sm/resLFC.SH2.csv")


#output
#log2 fold change (MAP): genotype ko vs wt 
#Wald test p-value: genotype ko vs wt 
#DataFrame with 5778 rows and 5 columns
#baseMean      log2FoldChange             lfcSE             pvalue              padj
#<numeric>           <numeric>         <numeric>          <numeric>         <numeric>
#  new_24_297_c          15.1799393964703   0.022986655206565 0.293537902620498  0.857630100586679 0.974031880688414
#new_150_409           15.6887670582187   0.142384921439735 0.321708579497533  0.316530705367445 0.770109064365736
#new_299_1067_c        65.3776259377723 -0.0351202342642024 0.253637191581099  0.823522576809578 0.969415454629426
#BSU00010               254.62410051172   0.158780784186938  0.25418877980532  0.355273722613138  0.79689438067828
#new_1751_1938         43.7790438759012  0.0588201209902989  0.26200204306265  0.710730217728052 0.938075438611867
#...                                


################################################################
# p-values and adjusted p-values
################################################################


resOrdered.SH2 <- res.SH2[order(res.SH2$pvalue),]
summary(res.SH2)
#output
#out of 5778 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 128, 2.2%
#LFC < 0 (down)     : 76, 1.3%
#outliers [1]       : 7, 0.12%
#low counts [2]     : 0, 0%
#(mean count < 2)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
sum(res.SH2$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
#[1] 204


## adjusted p-value < 0.05
res05 <- results(dds.SH2, alpha=0.05) ## adjusted p-value < 0.05
summary(res05)
#output
#out of 5778 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 104, 1.8%
#LFC < 0 (down)     : 64, 1.1%
#outliers [1]       : 7, 0.12%
#low counts [2]     : 337, 5.8%
#(mean count < 5)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
sum(res05$padj < 0.05, na.rm=TRUE)
#[1] 168


## adjusted p-value < 0.01
res01 <- results(dds.SH2, alpha=0.01) ## adjusted p-value < 0.01
summary(res01)
#output
#out of 5778 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 50, 0.87%
#LFC < 0 (down)     : 32, 0.55%
#outliers [1]       : 7, 0.12%
#low counts [2]     : 897, 16%
#(mean count < 11)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
sum(res01$padj < 0.01, na.rm=TRUE)
#[1] 82




################################################################
# Exploring and exporting results - MA-plot
################################################################

##Points which fall out of the window are plotted as open triangles pointing either up or down.
#Points will be colored red if the adjusted p value is less than 0.1
##It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
#plotMA(resLFC.SH2, ylim=c(-2,2), main = "MA plot \npadj-LFC less than 0.1 \nTranscriptomics")

#padj<0.05
#plotMA(res05, ylim=c(-2,2), main = "MA plot \npadj less than 0.05 \nTranscriptomics")
plotMA(resLFC.SH2, ylim=c(-2,2), alpha = 0.05, 
       main = "MA plot \nShrunken_LFC_padj<0.05 \nTranscriptomics")

#padj<0.01
#plotMA(res01, ylim=c(-2,2), main = "MA plot \npadj less than 0.01 \nTranscriptomics")
plotMA(resLFC.SH2, ylim=c(-2,2), alpha = 0.01, 
       main = "MA plot \nShrunken_LFC_padj<0.01 \nTranscriptomics")

#padj<0.001
#plotMA(res01, ylim=c(-2,2), main = "MA plot \npadj less than 0.01 \nTranscriptomics")
plotMA(resLFC.SH2, ylim=c(-2,2), alpha = 0.001, 
       main = "MA plot \nShrunken_LFC_padj<0.001 \nTranscriptomics")

#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
idx.SH2 <- identify(resLFC.SH2$baseMean, resLFC.SH2$log2FoldChange)

#Click on the points to get their row numbers
#once done, click finish in Plot viewer to get the row numbers

rownames(resLFC.SH2)[idx.SH2]  
#[1] "new_18866_18974"     "BSU01040"            "BSU01050"            "new_1041930_1041993" "BSU09670"            "BSU10410"           
#[7] "new_1117688_1117833" "new_1117924_1118790" "BSU15300"            "BSU17710"            "BSU21229"            "BSU21610"           
#[13] "BSU24690"            "BSU30560"            "BSU33770"            "BSU34930"            "BSU35609"            "BSU36310_copy2"     
#[19] "BSU37490"            "BSU39320"            "BSU40910"           
 

#> View(as.data.frame.matrix(resLFC.SH2))
#> View(as.data.frame.matrix(resLFC.SH2[5082,]))
#> View(as.data.frame.matrix(resLFC.SH2[4778,]))
#> View(as.data.frame.matrix(resLFC.SH2[2177,]))




###############################################################
# Plot counts
################################################################

#gene with minimum p-val
#plotCounts(dds.SH2, gene=which.min(res.SH2$padj), intgroup="genotype") 
#You can select the gene to plot by rowname or by numeric index.
plotCounts(dds.SH2, gene=which.min(res05$padj), intgroup="genotype") 
plotCounts(dds.SH2, gene="BSU00490", intgroup="genotype")
plotCounts(dds.SH2, gene="BSU00200", intgroup="genotype")
plotCounts(dds.SH2, gene="BSU13190", intgroup="genotype")




#customized plotting

d.SH2 <- plotCounts(dds.SH2, gene=which.min(res05$padj), intgroup="genotype", 
                   returnData=TRUE)
library("ggplot2")
ggplot(d.SH2, aes(x=genotype, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))



################################################################
# More information on results columns
################################################################

mcols(res.SH2)$description



#Rich visualization and reporting of results
#ReportingTools
#http://bioconductor.org/packages/release/bioc/html/ReportingTools.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ReportingTools", version = "3.8")



################################################################
# Exporting results to CSV files
################################################################


##write.csv(as.data.frame(resOrdered), 
##          file="diff.exp.genes_SH2.csv")         #prints all the differentially expressed genes
#resSig.SH2


resSig.SH2 = subset(resOrdered.SH2)
#write.csv(resSig.SH2, "../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_vs_wt_all_SH2_sm.csv")


resSig.SH2.05 <- subset(resOrdered.SH2, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#write.csv(as.data.frame(resSig.SH2.05), 
#          file="DEG_ko_vs_wt_p0.05_SH2_sm.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.SH2.01 <- subset(resOrdered.SH2, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.SH2.01
#write.csv(as.data.frame(resSig.SH2.01), 
#          file="DEG_ko_vs_wt_p0.01_SH2_sm.csv")



###################################  ALL GOOD SO FAR  ################################### 


#pheatmap(counts.all_SM_SH2,
#  scale="row",
#  cluster_rows=F, cluster_cols=T,
#  show_rownames = F, show_colnames = T, 
#  fontsize=10, legend=TRUE,
#  main = "HeatMap \nTotal Read Counts - SH2\nTranscriptomic"
#)


####################################################################
###4.PCA Analysis  
####################################################################

pca = prcomp(t(counts.all_SM_SH2), center = TRUE, scale = FALSE)
percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
pca2 = cbind(as.data.frame(pca$x), sampleTable_SM_SH2)   #check the order of ko and WT

#pdf("../Figures/LB/PCA_Plot_LB.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
  ylab(paste0("PC2: ",percent.var[2], "% variance")) + 
  labs(colour = "Genotype",shape="Media") +
  ggtitle("PCA Plot - SH2 \nTranscriptomic Analysis") + #CHANGE THE NAME HERE
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  ) + 
  #geom_text(aes(PC1, PC2, colour=genotype),label=pca2$genotype) +    #adds labels to each data point
  coord_fixed()
#dev.off()                  #turn this OFF if just want to see the picture in the Plots


