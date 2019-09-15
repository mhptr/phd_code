
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
library(tidyr)
library(EnhancedVolcano)


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
sampleTable_SM_M9 = sampleTable[1:6,]
sampleTable_SM_SH2 = sampleTable[7:12,]
sampleTable_SM_SH5 = sampleTable[13:18,]

#Reading Combined Read Counts 
counts.all = read.csv("./counts.all_SM_S1583.csv", header = T) #change it 
rownames(counts.all) = counts.all$X
counts.all = counts.all[,-1]
colnames(counts.all) = gsub(".txt","",colnames(counts.all))
colnames(counts.all) = sampleTable$sampleName
counts.all_SM_M9 = counts.all[,1:6]
counts.all_SM_SH2 = counts.all[,7:12]
counts.all_SM_SH5 = counts.all[,13:18]


#crude heatMap 

pheatmap(#counts.all, 
  counts.all_SM_M9,
  #counts.all_SM_SH2,
  #counts.all_SM_SH5,
  scale="row",
  cluster_rows=F, cluster_cols=T,
  show_rownames = F, show_colnames = T, 
  fontsize=10, legend=TRUE,
  #main = "HeatMap \nTotal Read Counts - All\nTranscriptomic"
  main = "HeatMap \nTotal Read Counts - M9\nTranscriptomic"
  #main = "HeatMap \nTotal Read Counts - SH2\nTranscriptomic"
  #main = "HeatMap \nTotal Read Counts - SH5\nTranscriptomic"
)


####################################################################
###4.PCA Analysis  
####################################################################


######################################
### 4.1 PCA - All
######################################

pca = prcomp(t(counts.all), center = TRUE, scale = FALSE)
percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
pca2 = cbind(as.data.frame(pca$x), sampleTable)   #check the order of ko and WT

#write.csv(pca2, "./Combined /Filtered /CSV/20181217_PCA_Combined.csv")

#pdf("../Figures/LB/PCA_Plot_LB.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype,shape=pca2$media)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
  ylab(paste0("PC2: ",percent.var[2], "% variance")) +
  labs(colour = "Genotype",shape="Media") +
  ggtitle("PCA Plot - All Conditions \nTranscriptomic Analysis - Raw Data") + #CHANGE THE NAME HERE 
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  ) + 
  #geom_text(aes(PC1, PC2, colour=conditions),label=pca2$genotype) +    #adds labels to each data point
  coord_fixed()
#dev.off()                  #turn this OFF if just want to see the picture in the Plots




################################################################
#4.2 PCA - Conditionwise - M9 
################################################################

pca = prcomp(t(counts.all_SM_M9), center = TRUE, scale = FALSE)
percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
pca2 = cbind(as.data.frame(pca$x), sampleTable_SM_M9)   #check the order of ko and WT

#pdf("../Figures/LB/PCA_Plot_LB.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype)) + geom_point(size=3) + 
  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
  ylab(paste0("PC2: ",percent.var[2], "% variance")) + 
  labs(colour = "Genotype",shape="Media") +
  ggtitle("PCA Plot - M9 \nTranscriptomic Analysis - Raw Data") + #CHANGE THE NAME HERE
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

cts.M9 = counts.all_SM_M9
##  Note: 
#In order to benefit from the default settings of the package, 
#you should put the variable of interest at the end of the formula and make sure the control level is the first level.
cts.M9 = cts.M9[,cbind( "S_M9_W1", "S_M9_W2", "S_M9_W3",
                        "S_M9_DEL1", "S_M9_DEL2", "S_M9_DEL3")]
colData.M9 = sampleTable_SM_M9    #Match order with column order of cts.M9


#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. 
rownames(colData.M9) = colData.M9$sampleName

all(rownames(colData.M9) %in% colnames(cts.M9)) # should be TRUE
all(rownames(colData.M9) == colnames(cts.M9)) # should be FALSE  
cts.M9 <- cts.M9[, rownames(colData.M9)]
all(rownames(colData.M9) == colnames(cts.M9)) # should be TRUE

library("DESeq2")
dds.M9 <- DESeqDataSetFromMatrix(countData = cts.M9,
                                 colData = colData.M9,
                                 design = ~ genotype)
dds.M9


################################################################
#5. DESeq2 - Pre-filtering
################################################################

keep.M9 = rowSums(counts(dds.M9)) >= 20
dds.M9 = dds.M9[keep.M9,]


#### Define factor levels

# Method 1 - Using Factor
#dds.M9$genotype <- factor(dds.M9$genotype, levels = c("wt","ko"))

### OR ###

# Method 2 - Using relevel
dds.M9$genotype = relevel(dds.M9$genotype, ref = "wt")


################################################################
#5. DESeq2 - Differential expression analysis
################################################################

dds.M9 <- DESeq(dds.M9)      # DESeq is DONE at this step!! 
res.M9 <- results(dds.M9)
res.M9                    ## WORKS! ko vs wt = i.e., condition treated vs untreated

#log2 fold change (MLE): genotype ko vs wt 
#Wald test p-value: genotype ko vs wt 
#DataFrame with 6094 rows and 6 columns


################################################################
# Log fold change shrinkage for visualization and ranking
################################################################

resultsNames(dds.M9)
resLFC.M9 <- lfcShrink(dds.M9, coef="genotype_ko_vs_wt", type="apeglm") #apeglm method for effect size shrinkage
resLFC.M9



################################################################
# p-values and adjusted p-values
################################################################


resOrdered.M9 <- res.M9[order(res.M9$pvalue),]
summary(res.M9)
sum(res.M9$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
res05 <- results(dds.M9, alpha=0.05) ## adjusted p-value < 0.05
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)


################################################################
# Exploring and exporting results - MA-plot
################################################################

##Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res.M9, ylim=c(-2,2)) #Points will be colored red if the adjusted p value is less than 0.1
##It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
plotMA(resLFC.M9, ylim=c(-2,2))
#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
#idx.M9 <- identify(res.M9$baseMean, res.M9$log2FoldChange)
#rownames(res.M9)[idx.M9]  # Doesn't end



################################################################
# Plot counts
################################################################

plotCounts(dds.M9, gene=which.min(res.M9$padj), intgroup="genotype")
d.M9 <- plotCounts(dds.M9, gene=which.min(res.M9$padj), intgroup="genotype", 
                   returnData=TRUE)
library("ggplot2")
ggplot(d.M9, aes(x=genotype, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))



################################################################
# More information on results columns
################################################################

mcols(res.M9)$description



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
##          file="DEG_ko_vs_wt_M9_sm.csv")         #prints all the differentially expressed genes
#resSig.M9

resSig.M9 = subset(resOrdered.M9)
#write.csv(resSig.M9, "../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_vs_wt_all_M9_sm.csv")

resSig.M9.05 <- subset(resOrdered.M9, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#write.csv(as.data.frame(resSig.M9.05), 
#          file="DEG_ko_vs_wt_p0.05_M9_sm.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.01 <- subset(resOrdered.M9, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.M9.01
#write.csv(as.data.frame(resSig.M9.01), 
#          file="DEG_ko_vs_wt_p0.01_M9_sm.csv")




##########################################################
## VENN
##########################################################


resSig.M9.05 = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_M9_minus_spoVG_upp.csv", header = T)
rownames(resSig.M9.05) = resSig.M9.05$X
resSig.SH2.05 = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_SH2_minus_spoVG_upp.csv", header = T)
rownames(resSig.SH2.05) = resSig.SH2.05$X
resSig.SH5.05 = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_SH5_minus_spoVG_upp.csv", header = T)
rownames(resSig.SH5.05) = resSig.SH5.05$X

venn.plot = venn.diagram(
  list("M9_Transcriptomic" = rownames(resSig.M9.05),
       "SH2_Transcriptomic" = rownames(resSig.SH2.05),
       "SH5_Transcriptomic" = rownames(resSig.SH5.05)),
  filename = NULL,
  col="transparent",
  fill=swati.color[1:3], height = 3000, width = 3000,resolution =500, imagetype = "tiff",
  main = "Venn Diagram - All Conditions \nTranscriptomic Analysis",main.pos  = c(0.5, 1.05), main.fontface = "bold",
  main.fontfamily = "serif", main.col = "Dark Blue", main.cex = 1.5
)
grid.newpage()
grid.draw(venn.plot)



#Getting intersects of overlaps
#https://stackoverflow.com/questions/23559371/how-to-get-the-list-of-items-in-venn-diagram-in-r

library(reshape)
library(R.utils)

## data
M9_Transcriptomic = data.frame("M9_Transcriptomic" = rownames(resSig.M9.05), "M9_Transcriptomic" = 1)
SH2_Transcriptomic = data.frame("SH2_Transcriptomic" = rownames(resSig.SH2.05), "SH2_Transcriptomic" = 1)
SH5_Transcriptomic = data.frame("SH5_Transcriptomic" = rownames(resSig.SH5.05), "SH5_Transcriptomic" = 1)


## a merged data frame.
Transcripts <- list(M9_Transcriptomic = M9_Transcriptomic, 
                    SH2_Transcriptomic = SH2_Transcriptomic,
                    SH5_Transcriptomic = SH5_Transcriptomic)
merged_Transcripts <- merge_recurse(Transcripts)

###################################  ALL GOOD SO FAR  ################################### 


#14/07/19

#################################################################
###                   Volcano Plots
##################################################################

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/filt_txt_all_2/")
library(ggplot2)
# volcano_trans = function(DEG, tag) {
#   DEG_File_Location = paste("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/", DEG, sep = "")
#   DEG_File = read.csv(DEG_File_Location, header = T)
#   rownames(DEG_File) = DEG_File$X
#   DEG_File = DEG_File[,-1]
#   
#   adj.P.Val =0.05
#   logFC=0.05
#   g = ggplot(data=DEG_File,
#     aes(title, x=logFC, y=-log10(adj.P.Val), colour=adj.P.Val<0.05)) +   
#     #colour=-log10(P.Value) +
#     ggtitle("Volcano Plot - Transcriptomics Analysis") +      #CHANGE THE NAME HERE 
#     theme(
#       plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
#       axis.title.x = element_text(color="black", size=12, face="bold"),
#       axis.title.y = element_text(color="black", size=12, face="bold") 
#     ) + 
#     geom_hline(yintercept=-log10(adj.P.Val), linetype="dashed", color = "red") + #thr.fc is 2
#         geom_point(alpha=0.4, size=1.75) +
#     ylim(c(0,23))+  xlim(c(-12,12))+   #LB, M9, SH2
#     xlab("log2 Fold Change (Del_spoVG vs WT)") + ylab("-log10 adj.P.Val")
#   print(g)
# }
# volcano_trans("DEG_ko_vs_wt_p0.05_M9_minus_spoVG_upp.csv", "M9")


DEG = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)

#12/09/19

thr.adj.pv =0.05

g = ggplot( data = DEG, 
  #data=diffexp.pr.all.SH2, 
  #data=diffexp.pr.all.SH5,
  aes(title, x=log2FoldChange, y=-log10(padj), colour= padj<0.05)) +   #colour=-log10(P.Value)
  
  ggtitle("Volcano Plot - M9 \nProteomic Analysis") + #CHANGE THE NAME HERE
  #ggtitle("Volcano Plot - SH2 \nProteomic Analysis") +
  #ggtitle("Volcano Plot - SH5 \nProteomic Analysis") +
  
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold") 
  ) + 
  geom_hline(yintercept=-log10(thr.adj.pv), linetype="dashed", color = "red") +
  #geom_vline(xintercept=-log2(thr.fc), linetype="dashed", color = "red") +  #thr.fc is 2
  #geom_vline(xintercept=log2(thr.fc), linetype="dashed", color = "red") +
  geom_point(alpha=0.4, size=1.75) +
  #ylim(c(0,23))+  xlim(c(-12,12))+   #LB, M9, SH2
  #ylim(c(0,28))+  xlim(c(-14,14))+   #SH5
  xlab("log2 Fold Change (Del_spoVG vs WT)") + ylab("-log10 adj.P.Val")
print(g)









# https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html#introduction
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('EnhancedVolcano')
# devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
library(ggplot2)

EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-8, 8),
                title = 'Volcano Plot - Transcriptomics Analysis',
                pCutoff = 10e-16,
                FCcutoff = 1.5,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0)


EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-6, 6),
                ylim = c(0,77),
                title = 'Volcano Plot - Transcriptomics Analysis',
                pCutoff = 10e-16,
                FCcutoff = 1.5,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0,
                col=c('black', 'orange', 'green', 'red3'),
                colAlpha = 1)



EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-6, 6),
                pCutoff = 10e-12,
                FCcutoff = 1.5,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                transcriptPointSize = 3.0,
                transcriptLabSize = 4.0,
                colAlpha = 1,
                legend=c('NS','Log (base 2) fold-change','P value',
                         'P value & Log (base 2) fold-change'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)


EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('new_3788276_3788397_c', 'BSU17180', 'BSU00500'),
                xlim = c(-5.5,8),
                ylim = c(0, 25),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-14,
                FCcutoff = 2.0,
                transcriptPointSize = 3.0,
                transcriptLabSize = 5.0,
                transcriptLabCol = 'black',
                transcriptLabFace = 'bold',
                boxedlabels = TRUE,
                colAlpha = 4/5,
                legend=c('NS','Log (base 2) fold-change','padj',
                         'padj & Log (base 2) fold-change'),
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')


##

#https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html#adjust-colour-and-alpha-for-point-shading

DEG = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)
rownames(DEG) = DEG$
DEG = DEG[,-1]


EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-4.5, 4.5),
                ylim = c(0, 50),
                title = 'N061011 versus N61311',
                pCutoff = 10e-16,
                FCcutoff = 1.5,
                #pointSize = c(ifelse(DEG$log2FoldChange>2, 8,1)),
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)




###################

#15/07/19

library(ggplot2)
library(scales)
library(limma)
# # set up an example dataset (see ?lmFit)
# sd <- 0.3*sqrt(4/rchisq(100,df=4))
# y <- matrix(rnorm(100*6,sd=sd),100,6)
# rownames(y) <- paste("Gene",1:100)
# y[1:2,4:6] <- y[1:2,4:6] + 2
# design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
# options(digits=3)

# Ordinary fit

dim(DEG$log2FoldChange)

fit <- lmFit(DEG$log2FoldChange, DEG$padj)
fit <- eBayes(fit)
tt <- topTable(fit,coef=2, n = Inf)

# transformation function for reverse log1p axis
revlog_trans <- function(base = exp(1)) {
  trans <- function(x) -log1p(x)
  inv <- function(x) expm1(-x)
  scales::trans_new("revlog1p", trans, inv, domain = c(0, Inf))
}

ggplot(tt, aes(x = logFC, y = P.Value)) +
  scale_fill_gradient(low = "lightgray", high = "navy") +
  scale_color_gradient(low = "lightgray", high = "navy") +
  scale_y_continuous(trans = revlog_trans(), expand = c(0.005, 0.005)) +
  expand_limits(y = c(0, 1)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon",
                  show.legend = FALSE) +
  geom_point(data = subset(tt, P.Value < 0.05), 
             color = "red", alpha = 0.5) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_linedraw() + 
  theme(panel.grid = element_blank()) +
  xlab("Fold change (log2)") +
  ylab("P-Value") +
  annotate("text", x = min(tt$logFC), y = 1,
           label = "Nominally significant",
           color = "black", hjust = 0) +
  annotate("point", x = min(tt$logFC) - 0.05, y = 1, color = "red")


######################################################################################
### 5. Installing KEGGREST  #No need to run - Already downloaded csv exists  
######################################################################################


## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#pkgs <- rownames(installed.packages())
#biocLite(pkgs, type="source")                     #something worked 


## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("KEGGREST")                              #worked 

#library(KEGGREST)
#bsu.geneName = keggList("bsu")
#bsu.geneName.df = data.frame(geneName = bsu.geneName)
#bsu.geneName.df$geneID = names(bsu.geneName)

#bsu.uniprot = keggConv("bsu", "uniprot")
#bsu.uniprot.df = data.frame(geneID = bsu.uniprot)
#bsu.uniprot.df$uniProt = names(bsu.uniprot)

#bsu.pathways.gene = keggLink("pathway", "bsu")
#bsu.pathways.gene.df = data.frame(keggPathwayID = bsu.pathways.gene)
#bsu.pathways.gene.df$geneID = names(bsu.pathways.gene)

#bsu.keggPathName = keggList("pathway", "bsu")
#bsu.keggPathName.df = data.frame(keggPathName = bsu.keggPathName)
#bsu.keggPathName.df$keggPathwayID = names(bsu.keggPathName)
#bsu.keggPathName.df$keggPathName = gsub(" - Bacillus subtilis subsp. subtilis 168","",bsu.keggPathName.df$keggPathName)

#bsu.kegg.pathway = merge(merge(bsu.geneName.df,bsu.uniprot.df,by="geneID"),bsu.pathways.gene.df,by="geneID")
#bsu.kegg.pathway$keggPathName = unlist(lapply(bsu.kegg.pathway$keggPathwayID, function(x)
#  bsu.keggPathName.df$keggPathName[bsu.keggPathName.df$keggPathwayID==x]))
#bsu.kegg.pathway$geneID = gsub("bsu:","",bsu.kegg.pathway$geneID)
#bsu.kegg.pathway$uniProt = gsub("up:","",bsu.kegg.pathway$uniProt)
#bsu.kegg.pathway$keggPathwayID = gsub("path:","",bsu.kegg.pathway$keggPathwayID)

#write.csv(bsu.kegg.pathway,"./data/bsu.kegg.pathway.csv",row.names = F)

######################################################################################
### 5. KEGG Pathway Enrichment
######################################################################################

library(pheatmap)

DEG_all_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)
rownames(DEG_all_M9) = DEG_all_M9$X
DEG_all_M9 = DEG_all_M9[,-1]

DEG_all_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH2_minus_spoVG_upp.csv", header = T)
rownames(DEG_all_SH2) = DEG_all_SH2$X
DEG_all_SH2 = DEG_all_SH2[,-1]

DEG_all_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH5_minus_spoVG_upp.csv", header = T)
rownames(DEG_all_SH5) = DEG_all_SH5$X
DEG_all_SH5 = DEG_all_SH5[,-1]

bsu.kegg.pathway = read.csv("../bsu.kegg.pathway.csv", header = T)
DEG_all = list(DEG_all_M9,DEG_all_SH2,DEG_all_SH5)
names(DEG_all) = c("M9","SH2","SH5")
kegg.enrich.all = c()                              #make dataframe 
for(i in 1:length(DEG_all)){
  my.de = DEG_all[[i]]
  my.de.genes = rownames(my.de)[my.de$padj<0.05]
  m = length(my.de.genes)                             #m = length of Sig Diff Ex pr
  n = dim(my.de)[1] - m                               #n = length of all diff. protein - length of Sig Diff Ex pr
  my.kegg.enrich = c()   
  for(j in 1:length(unique(bsu.kegg.pathway$keggPathName))) {
    my.path = unique(bsu.kegg.pathway$keggPathName)[j]
    my.path.gene = unique(bsu.kegg.pathway$geneID[bsu.kegg.pathway$keggPathName == my.path])
    k = length(my.path.gene)             # k = Total no of genes in pathway
    x = sum(my.de.genes%in%my.path.gene) # x = no of differentially expressed genes in a particular pathway
    my.de.path.gene = paste(my.de.genes[my.de.genes%in%my.path.gene],collapse = "|")
    pv = phyper(x-1,m,n,k,lower.tail = F) 
    my.kegg.enrich = rbind(my.kegg.enrich,cbind(pv,k,x,my.de.path.gene))
  }
  colnames(my.kegg.enrich) = paste0(names(DEG_all)[i],"_",colnames(my.kegg.enrich))
  kegg.enrich.all = as.data.frame(cbind(kegg.enrich.all,my.kegg.enrich)) 
}
rownames(kegg.enrich.all) = unique(bsu.kegg.pathway$keggPathName)
#kegg.enrich.all = kegg.enrich.all[!apply(kegg.enrich.all, 1, function(x) any(is.na(x))),]

kegg.enrich.all = kegg.enrich.all[as.numeric(as.character(kegg.enrich.all$M9_k))>2,]
kegg.enrich.sig = kegg.enrich.all[(as.numeric(as.character(kegg.enrich.all$M9_pv))<0.05|
                                     as.numeric(as.character(kegg.enrich.all$SH2_pv))<0.05|
                                     as.numeric(as.character(kegg.enrich.all$SH5_pv))<0.05),]
#write.csv(kegg.enrich.sig, "../diff.exp.gene/DEG_SM/DEG_KEGG/kegg.enrich.sig_transcriptomics.csv")

kegg.enrich.sig.pv = cbind(as.numeric(as.character(kegg.enrich.sig$M9_pv)),
                           as.numeric(as.character(kegg.enrich.sig$SH2_pv)),
                           as.numeric(as.character(kegg.enrich.sig$SH5_pv)))
rownames(kegg.enrich.sig.pv) = rownames(kegg.enrich.sig)
colnames(kegg.enrich.sig.pv) = c("KEGG_M9_adj.P.Val","KEGG_SH2_adj.P.Val","KEGG_SH5_adj.P.Val")
#write.csv(kegg.enrich.sig.pv, "../diff.exp.gene/DEG_SM/DEG_KEGG/kegg.enrich.sig.pv_transcriptomics.csv")

kegg.enrich.sig.pv.log10 = -log10(kegg.enrich.sig.pv)
kegg.enrich.sig.pv.log10[kegg.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!

pheatmap(kegg.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 70, cellheight = 27,
         show_rownames = T, show_colnames = T, 
         fontsize=12, legend=TRUE,
         headerplot = "KEEG_Transcriptomics_Analysis",
         main = "KEGG Pathway Enrichment Analysis \n(adj.P.Val) - Transcriptomics"
)

##################################################################################################
###### Gene Set Enrichment Analysis - Category 3 (GSEA-3) with adj.P.Val
##################################################################################################

#######24/07/19

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
library(tidyr)

####################################################################
### 2. Set Working Directory , Set Parameters 
####################################################################

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/filt_txt_all_2/")

thr.fc=2 #thresold_FoldChange
thr.pv =0.05 #thresold_pValue1
swati.color = c(brewer.pal(8,"Dark2")) #TA
swati.color2 = c(brewer.pal(9, "Greens")) #SM
swati.color.Blues = c(brewer.pal(8,"Blues"))

#source("http://www.bioconductor.org/biocLite.R")
#BiocManager::install("piano", dependencies=TRUE)          #this worked
#Package piano version 1.22.0 Index

geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = T)
myGsc.cat3 = loadGSC(cbind(as.character(geneCategories$gene),
                           as.character(geneCategories$category3)))

DEG_all_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)
rownames(DEG_all_M9) = DEG_all_M9$X
DEG_all_M9 = DEG_all_M9[,-1]

DEG_all_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH2_minus_spoVG_upp.csv", header = T)
rownames(DEG_all_SH2) = DEG_all_SH2$X
DEG_all_SH2 = DEG_all_SH2[,-1]

DEG_all_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH5_minus_spoVG_upp.csv", header = T)
rownames(DEG_all_SH5) = DEG_all_SH5$X
DEG_all_SH5 = DEG_all_SH5[,-1]

DEG_all = list(DEG_all_M9, DEG_all_SH2, DEG_all_SH5)
names(DEG_all) = c("M9","SH2","SH5")

gsea.enrich.cat3.all = data.frame(Name = as.character())
for (i in 1:length(DEG_all)){
  my.de = DEG_all[[i]]
  my.pval = my.de$padj
  my.fc = my.de$log2FoldChange
  length(my.pval)
  length(my.fc)
  names(my.pval) = rownames(my.de)
  names(my.fc) = rownames(my.de)
  my.fc = my.fc[!is.na(my.pval)]
  my.pval = my.pval[!is.na(my.pval)]
  length(my.pval)
  length(my.fc)
  my.gsaRes.cat3 = runGSA(geneLevelStats=my.pval,directions=my.fc,
                          gsc=myGsc.cat3,geneSetStat="reporter",adjMethod="BH")
  my.gsea.cat3 = GSAsummaryTable(my.gsaRes.cat3)
  my.gsea.cat3 = my.gsea.cat3[!my.gsea.cat3$Name == "",]  #removes the empty row
  #my.gsea.cat3 = my.gsea.cat3[order(my.gsea.cat3$Name),]  #sorts alphabetically
  colnames(my.gsea.cat3) = paste0(names(DEG_all)[i],"_",colnames(my.gsea.cat3)) #keep the first column as name 
  colnames(my.gsea.cat3)[1] = "Name"
  dim(my.gsea.cat3)
  gsea.enrich.cat3.all= merge(gsea.enrich.cat3.all,my.gsea.cat3, by="Name", all=T)
}
dim(gsea.enrich.cat3.all)
rownames(gsea.enrich.cat3.all) = gsea.enrich.cat3.all[,"Name"]
gsea.enrich.cat3.all = gsea.enrich.cat3.all[,-1]
#write.csv(gsea.enrich.cat3.all, "../diff.exp.gene/DEG_SM/DEG_GSEA-3/Final/gsea.enrich.cat3.all_transcriptomics.csv")

gsea.cat3.enrich.ND =  gsea.enrich.cat3.all[,grep("non-dir.", colnames(gsea.enrich.cat3.all))]
colnames(gsea.cat3.enrich.ND) = c("M9_Stat(non-dir.)",
                                  "M9_p(non-dir.)",
                                  "M9_padj(non-dir.)",
                                  "SH2_Stat(non-dir.)",
                                  "SH2_p(non-dir.)",
                                  "SH2_padj(non-dir.)",
                                  "SH5_Stat(non-dir.)",
                                  "SH5_p(non-dir.)",
                                  "SH5_padj(non-dir.)")

colnames(gsea.cat3.enrich.ND) = c("M9_Stat_ND",
                                  "M9_p_ND",
                                  "M9_padj_ND",
                                  "SH2_Stat_ND",
                                  "SH2_p_ND",
                                  "SH2_padj_ND",
                                  "SH5_Stat_ND",
                                  "SH5_p_ND",
                                  "SH5_padj_ND")

gsea.cat3.enrich = matrix(as.numeric(unlist(gsea.cat3.enrich.ND)), nrow = nrow(gsea.cat3.enrich.ND))
rownames(gsea.cat3.enrich) = rownames(gsea.cat3.enrich.ND)
colnames(gsea.cat3.enrich) = colnames(gsea.cat3.enrich.ND)
gsea.cat3.enrich.sig = gsea.cat3.enrich[(as.numeric(as.character(gsea.cat3.enrich[,"M9_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat3.enrich[,"SH2_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat3.enrich[,"SH5_padj_ND"]))<0.05),]
gsea.cat3.enrich.sig = na.omit(gsea.cat3.enrich.sig)
#write.csv(gsea.cat3.enrich.sig, "../diff.exp.gene/DEG_SM/DEG_GSEA-3/Final/gsea.cat3.enrich.sig_transcriptomics.csv")

gsea.cat3.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat3.enrich.sig[,"M9_padj_ND"])),
                                as.numeric(as.character(gsea.cat3.enrich.sig[,"SH2_padj_ND"])),
                                as.numeric(as.character(gsea.cat3.enrich.sig[,"SH5_padj_ND"])))
rownames(gsea.cat3.enrich.sig.pv) = rownames(gsea.cat3.enrich.sig)
colnames(gsea.cat3.enrich.sig.pv) = c("GSEA_Cat3_M9_padj","GSEA_Cat3_SH2_padj","GSEA_Cat3_SH5_padj")
gsea.cat3.enrich.sig.pv.log10 = -log10(gsea.cat3.enrich.sig.pv)
gsea.cat3.enrich.sig.pv.log10[gsea.cat3.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!

annotation.cat3.cat1 = data.frame(cat1 = unlist(lapply(rownames(gsea.cat3.enrich.sig.pv.log10), function(x)
   geneCategories$category1[geneCategories$category3==x][1])))
rownames(annotation.cat3.cat1) = rownames(gsea.cat3.enrich.sig.pv.log10)

#write.csv(annotation.cat3.cat1, "../diff.exp.gene/DEG_SM/DEG_GSEA-3/Final/annotation.cat3.cat1_transcriptomics.csv")

pheatmap(gsea.cat3.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 25, cellheight = 10,
         show_rownames = T, show_colnames = T, 
         fontsize=8, legend=TRUE,
         headerplot = "GSEA_Cat3_Transcriptomics_Analysis",
         main = "GSEA SubtiWiki Category 3 \n(padj) - Transcriptomics",
         annotation_row = annotation.cat3.cat1
)

##################################################################################################
###### Gene Set Enrichment Analysis - Category 1 (GSEA-1) with adj.P.Val
##################################################################################################

#######24/07/19

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
library(tidyr)

####################################################################
### 2. Set Working Directory , Set Parameters 
####################################################################

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/filt_txt_all_2/")

thr.fc=2 #thresold_FoldChange
thr.pv =0.05 #thresold_pValue1
swati.color = c(brewer.pal(8,"Dark2")) #TA
swati.color2 = c(brewer.pal(9, "Greens")) #SM
swati.color.Blues = c(brewer.pal(8,"Blues"))

#source("http://www.bioconductor.org/biocLite.R")
#BiocManager::install("piano", dependencies=TRUE)          #this worked
#Package piano version 1.22.0 Index

geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = T)
myGsc.cat1 = loadGSC(cbind(as.character(geneCategories$gene),
                           as.character(geneCategories$category1)))

DEG_all_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)
rownames(DEG_all_M9) = DEG_all_M9$X
DEG_all_M9 = DEG_all_M9[,-1]

DEG_all_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH2_minus_spoVG_upp.csv", header = T)
rownames(DEG_all_SH2) = DEG_all_SH2$X
DEG_all_SH2 = DEG_all_SH2[,-1]

DEG_all_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH5_minus_spoVG_upp.csv", header = T)
rownames(DEG_all_SH5) = DEG_all_SH5$X
DEG_all_SH5 = DEG_all_SH5[,-1]

DEG_all = list(DEG_all_M9, DEG_all_SH2, DEG_all_SH5)
names(DEG_all) = c("M9","SH2","SH5")

gsea.enrich.cat1.all = data.frame(Name = as.character())
for (i in 1:length(DEG_all)){
  my.de = DEG_all[[i]]
  my.pval = my.de$padj
  my.fc = my.de$log2FoldChange
  length(my.pval)
  length(my.fc)
  names(my.pval) = rownames(my.de)
  names(my.fc) = rownames(my.de)
  my.fc = my.fc[!is.na(my.pval)]
  my.pval = my.pval[!is.na(my.pval)]
  length(my.pval)
  length(my.fc)
  my.gsaRes.cat1 = runGSA(geneLevelStats=my.pval,directions=my.fc,
                          gsc=myGsc.cat1,geneSetStat="reporter",adjMethod="BH")
  my.gsea.cat1 = GSAsummaryTable(my.gsaRes.cat1)
  my.gsea.cat1 = my.gsea.cat1[!my.gsea.cat1$Name == "",]  #removes the empty row
  #my.gsea.cat1 = my.gsea.cat1[order(my.gsea.cat1$Name),]  #sorts alphabetically
  colnames(my.gsea.cat1) = paste0(names(DEG_all)[i],"_",colnames(my.gsea.cat1)) #keep the first column as name 
  colnames(my.gsea.cat1)[1] = "Name"
  dim(my.gsea.cat1)
  gsea.enrich.cat1.all= merge(gsea.enrich.cat1.all,my.gsea.cat1, by="Name", all=T)
}
dim(gsea.enrich.cat1.all)
rownames(gsea.enrich.cat1.all) = gsea.enrich.cat1.all[,"Name"]
gsea.enrich.cat1.all = gsea.enrich.cat1.all[,-1]
#write.csv(gsea.enrich.cat1.all, "../diff.exp.gene/DEG_SM/DEG_GSEA-1/Final/gsea.enrich.cat1.all_transcriptomics.csv")

gsea.cat1.enrich.ND =  gsea.enrich.cat1.all[,grep("non-dir.", colnames(gsea.enrich.cat1.all))]
colnames(gsea.cat1.enrich.ND) = c("M9_Stat(non-dir.)",
                                  "M9_p(non-dir.)",
                                  "M9_padj(non-dir.)",
                                  "SH2_Stat(non-dir.)",
                                  "SH2_p(non-dir.)",
                                  "SH2_padj(non-dir.)",
                                  "SH5_Stat(non-dir.)",
                                  "SH5_p(non-dir.)",
                                  "SH5_padj(non-dir.)")

colnames(gsea.cat1.enrich.ND) = c("M9_Stat_ND",
                                  "M9_p_ND",
                                  "M9_padj_ND",
                                  "SH2_Stat_ND",
                                  "SH2_p_ND",
                                  "SH2_padj_ND",
                                  "SH5_Stat_ND",
                                  "SH5_p_ND",
                                  "SH5_padj_ND")

gsea.cat1.enrich = matrix(as.numeric(unlist(gsea.cat1.enrich.ND)), nrow = nrow(gsea.cat1.enrich.ND))
rownames(gsea.cat1.enrich) = rownames(gsea.cat1.enrich.ND)
colnames(gsea.cat1.enrich) = colnames(gsea.cat1.enrich.ND)
gsea.cat1.enrich.sig = gsea.cat1.enrich[(as.numeric(as.character(gsea.cat1.enrich[,"M9_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat1.enrich[,"SH2_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat1.enrich[,"SH5_padj_ND"]))<0.05),]
gsea.cat1.enrich.sig = na.omit(gsea.cat1.enrich.sig)
#write.csv(gsea.cat1.enrich.sig, "../diff.exp.gene/DEG_SM/DEG_GSEA-1/Final/gsea.cat1.enrich.sig_transcriptomics.csv")

gsea.cat1.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat1.enrich.sig[,"M9_padj_ND"])),
                                as.numeric(as.character(gsea.cat1.enrich.sig[,"SH2_padj_ND"])),
                                as.numeric(as.character(gsea.cat1.enrich.sig[,"SH5_padj_ND"])))
rownames(gsea.cat1.enrich.sig.pv) = rownames(gsea.cat1.enrich.sig)
colnames(gsea.cat1.enrich.sig.pv) = c("GSEA_cat1_M9_padj","GSEA_cat1_SH2_padj","GSEA_cat1_SH5_padj")
gsea.cat1.enrich.sig.pv.log10 = -log10(gsea.cat1.enrich.sig.pv)
gsea.cat1.enrich.sig.pv.log10[gsea.cat1.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!

pheatmap(gsea.cat1.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 25, cellheight = 10,
         show_rownames = T, show_colnames = T, 
         fontsize=8, legend=TRUE,
         headerplot = "GSEA_cat1_Transcriptomics_Analysis",
         main = "GSEA SubtiWiki Category 1 \n(padj) - Transcriptomics"
)



##################################################################################################
###### Gene Set Enrichment Analysis - Category 2 (GSEA-2) with adj.P.Val
##################################################################################################

#######25/07/19

geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = T)
myGsc.cat2 = loadGSC(cbind(as.character(geneCategories$gene),
                           as.character(geneCategories$category2)))
# DEG_all_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_M9) = DEG_all_M9$X
# DEG_all_M9 = DEG_all_M9[,-1]
# DEG_all_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH2_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_SH2) = DEG_all_SH2$X
# DEG_all_SH2 = DEG_all_SH2[,-1]# 
# DEG_all_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH5_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_SH5) = DEG_all_SH5$X
# DEG_all_SH5 = DEG_all_SH5[,-1]
DEG_all = list(DEG_all_M9, DEG_all_SH2, DEG_all_SH5)
names(DEG_all) = c("M9","SH2","SH5")

gsea.enrich.cat2.all = data.frame(Name = as.character())
for (i in 1:length(DEG_all)){
  my.de = DEG_all[[i]]
  my.pval = my.de$padj
  my.fc = my.de$log2FoldChange
  length(my.pval)
  length(my.fc)
  names(my.pval) = rownames(my.de)
  names(my.fc) = rownames(my.de)
  my.fc = my.fc[!is.na(my.pval)]
  my.pval = my.pval[!is.na(my.pval)]
  length(my.pval)
  length(my.fc)
  my.gsaRes.cat2 = runGSA(geneLevelStats=my.pval,directions=my.fc,
                          gsc=myGsc.cat2,geneSetStat="reporter",adjMethod="BH")
  my.gsea.cat2 = GSAsummaryTable(my.gsaRes.cat2)
  my.gsea.cat2 = my.gsea.cat2[!my.gsea.cat2$Name == "",]  #removes the empty row
  #my.gsea.cat2 = my.gsea.cat2[order(my.gsea.cat2$Name),]  #sorts alphabetically
  colnames(my.gsea.cat2) = paste0(names(DEG_all)[i],"_",colnames(my.gsea.cat2)) #keep the first column as name 
  colnames(my.gsea.cat2)[1] = "Name"
  dim(my.gsea.cat2)
  gsea.enrich.cat2.all= merge(gsea.enrich.cat2.all,my.gsea.cat2, by="Name", all=T)
}
dim(gsea.enrich.cat2.all)
rownames(gsea.enrich.cat2.all) = gsea.enrich.cat2.all[,"Name"]
gsea.enrich.cat2.all = gsea.enrich.cat2.all[,-1]
#write.csv(gsea.enrich.cat2.all, "../diff.exp.gene/DEG_SM/DEG_GSEA-2/Final/gsea.enrich.cat2.all_transcriptomics.csv")

gsea.cat2.enrich.ND =  gsea.enrich.cat2.all[,grep("non-dir.", colnames(gsea.enrich.cat2.all))]
colnames(gsea.cat2.enrich.ND) = c("M9_Stat(non-dir.)",
                                  "M9_p(non-dir.)",
                                  "M9_padj(non-dir.)",
                                  "SH2_Stat(non-dir.)",
                                  "SH2_p(non-dir.)",
                                  "SH2_padj(non-dir.)",
                                  "SH5_Stat(non-dir.)",
                                  "SH5_p(non-dir.)",
                                  "SH5_padj(non-dir.)")

colnames(gsea.cat2.enrich.ND) = c("M9_Stat_ND",
                                  "M9_p_ND",
                                  "M9_padj_ND",
                                  "SH2_Stat_ND",
                                  "SH2_p_ND",
                                  "SH2_padj_ND",
                                  "SH5_Stat_ND",
                                  "SH5_p_ND",
                                  "SH5_padj_ND")

gsea.cat2.enrich = matrix(as.numeric(unlist(gsea.cat2.enrich.ND)), nrow = nrow(gsea.cat2.enrich.ND))
rownames(gsea.cat2.enrich) = rownames(gsea.cat2.enrich.ND)
colnames(gsea.cat2.enrich) = colnames(gsea.cat2.enrich.ND)
gsea.cat2.enrich.sig = gsea.cat2.enrich[(as.numeric(as.character(gsea.cat2.enrich[,"M9_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat2.enrich[,"SH2_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat2.enrich[,"SH5_padj_ND"]))<0.05),]
gsea.cat2.enrich.sig = na.omit(gsea.cat2.enrich.sig)
#write.csv(gsea.cat2.enrich.sig, "../diff.exp.gene/DEG_SM/DEG_GSEA-2/Final/gsea.cat2.enrich.sig_transcriptomics.csv")

gsea.cat2.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat2.enrich.sig[,"M9_padj_ND"])),
                                as.numeric(as.character(gsea.cat2.enrich.sig[,"SH2_padj_ND"])),
                                as.numeric(as.character(gsea.cat2.enrich.sig[,"SH5_padj_ND"])))
rownames(gsea.cat2.enrich.sig.pv) = rownames(gsea.cat2.enrich.sig)
colnames(gsea.cat2.enrich.sig.pv) = c("GSEA_cat2_M9_padj","GSEA_cat2_SH2_padj","GSEA_cat2_SH5_padj")
gsea.cat2.enrich.sig.pv.log10 = -log10(gsea.cat2.enrich.sig.pv)
gsea.cat2.enrich.sig.pv.log10[gsea.cat2.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!

annotation.cat2.cat1 = data.frame(cat1 = unlist(lapply(rownames(gsea.cat2.enrich.sig.pv.log10), function(x)
  geneCategories$category1[geneCategories$category3==x][1])))
rownames(annotation.cat2.cat1) = rownames(gsea.cat2.enrich.sig.pv.log10)
#write.csv(annotation.cat2.cat1, "../diff.exp.gene/DEG_SM/DEG_GSEA-2/Final/annotation.cat2.cat1_transcriptomics.csv")

pheatmap(gsea.cat2.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 25, cellheight = 10,
         show_rownames = T, show_colnames = T, 
         fontsize=8, legend=TRUE,
         headerplot = "GSEA_cat2_Transcriptomics_Analysis",
         main = "GSEA SubtiWiki Category 2 \n(padj) - Transcriptomics"
)

##################################################################################################
###### Gene Set Enrichment Analysis - Category 4 (GSEA-4) with adj.P.Val
##################################################################################################

#######25/07/19

geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = T)
myGsc.cat4 = loadGSC(cbind(as.character(geneCategories$gene),
                           as.character(geneCategories$category4)))
# DEG_all_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_M9) = DEG_all_M9$X
# DEG_all_M9 = DEG_all_M9[,-1]
# DEG_all_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH2_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_SH2) = DEG_all_SH2$X
# DEG_all_SH2 = DEG_all_SH2[,-1]# 
# DEG_all_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH5_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_SH5) = DEG_all_SH5$X
# DEG_all_SH5 = DEG_all_SH5[,-1]
DEG_all = list(DEG_all_M9, DEG_all_SH2, DEG_all_SH5)
names(DEG_all) = c("M9","SH2","SH5")

gsea.enrich.cat4.all = data.frame(Name = as.character())
for (i in 1:length(DEG_all)){
  my.de = DEG_all[[i]]
  my.pval = my.de$padj
  my.fc = my.de$log2FoldChange
  length(my.pval)
  length(my.fc)
  names(my.pval) = rownames(my.de)
  names(my.fc) = rownames(my.de)
  my.fc = my.fc[!is.na(my.pval)]
  my.pval = my.pval[!is.na(my.pval)]
  length(my.pval)
  length(my.fc)
  my.gsaRes.cat4 = runGSA(geneLevelStats=my.pval,directions=my.fc,
                          gsc=myGsc.cat4,geneSetStat="reporter",adjMethod="BH")
  my.gsea.cat4 = GSAsummaryTable(my.gsaRes.cat4)
  my.gsea.cat4 = my.gsea.cat4[!my.gsea.cat4$Name == "",]  #removes the empty row
  #my.gsea.cat4 = my.gsea.cat4[order(my.gsea.cat4$Name),]  #sorts alphabetically
  colnames(my.gsea.cat4) = paste0(names(DEG_all)[i],"_",colnames(my.gsea.cat4)) #keep the first column as name 
  colnames(my.gsea.cat4)[1] = "Name"
  dim(my.gsea.cat4)
  gsea.enrich.cat4.all= merge(gsea.enrich.cat4.all,my.gsea.cat4, by="Name", all=T)
}
dim(gsea.enrich.cat4.all)
rownames(gsea.enrich.cat4.all) = gsea.enrich.cat4.all[,"Name"]
gsea.enrich.cat4.all = gsea.enrich.cat4.all[,-1]
#write.csv(gsea.enrich.cat4.all, "../diff.exp.gene/DEG_SM/DEG_GSEA-4/Final/gsea.enrich.cat4.all_transcriptomics.csv")

gsea.cat4.enrich.ND =  gsea.enrich.cat4.all[,grep("non-dir.", colnames(gsea.enrich.cat4.all))]
colnames(gsea.cat4.enrich.ND) = c("M9_Stat(non-dir.)",
                                  "M9_p(non-dir.)",
                                  "M9_padj(non-dir.)",
                                  "SH2_Stat(non-dir.)",
                                  "SH2_p(non-dir.)",
                                  "SH2_padj(non-dir.)",
                                  "SH5_Stat(non-dir.)",
                                  "SH5_p(non-dir.)",
                                  "SH5_padj(non-dir.)")

colnames(gsea.cat4.enrich.ND) = c("M9_Stat_ND",
                                  "M9_p_ND",
                                  "M9_padj_ND",
                                  "SH2_Stat_ND",
                                  "SH2_p_ND",
                                  "SH2_padj_ND",
                                  "SH5_Stat_ND",
                                  "SH5_p_ND",
                                  "SH5_padj_ND")

gsea.cat4.enrich = matrix(as.numeric(unlist(gsea.cat4.enrich.ND)), nrow = nrow(gsea.cat4.enrich.ND))
rownames(gsea.cat4.enrich) = rownames(gsea.cat4.enrich.ND)
colnames(gsea.cat4.enrich) = colnames(gsea.cat4.enrich.ND)
gsea.cat4.enrich.sig = gsea.cat4.enrich[(as.numeric(as.character(gsea.cat4.enrich[,"M9_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat4.enrich[,"SH2_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat4.enrich[,"SH5_padj_ND"]))<0.05),]
gsea.cat4.enrich.sig = na.omit(gsea.cat4.enrich.sig)
#write.csv(gsea.cat4.enrich.sig, "../diff.exp.gene/DEG_SM/DEG_GSEA-4/Final/gsea.cat4.enrich.sig_transcriptomics.csv")

gsea.cat4.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat4.enrich.sig[,"M9_padj_ND"])),
                                as.numeric(as.character(gsea.cat4.enrich.sig[,"SH2_padj_ND"])),
                                as.numeric(as.character(gsea.cat4.enrich.sig[,"SH5_padj_ND"])))
rownames(gsea.cat4.enrich.sig.pv) = rownames(gsea.cat4.enrich.sig)
colnames(gsea.cat4.enrich.sig.pv) = c("GSEA_cat4_M9_padj","GSEA_cat4_SH2_padj","GSEA_cat4_SH5_padj")
gsea.cat4.enrich.sig.pv.log10 = -log10(gsea.cat4.enrich.sig.pv)
gsea.cat4.enrich.sig.pv.log10[gsea.cat4.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!

annotation.cat4.cat1 = data.frame(cat1 = unlist(lapply(rownames(gsea.cat4.enrich.sig.pv.log10), function(x)
  geneCategories$category1[geneCategories$category3==x][1])))
rownames(annotation.cat4.cat1) = rownames(gsea.cat4.enrich.sig.pv.log10)
#write.csv(annotation.cat4.cat1, "../diff.exp.gene/DEG_SM/DEG_GSEA-4/Final/annotation.cat4.cat1_transcriptomics.csv")

pheatmap(gsea.cat4.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 25, cellheight = 10,
         show_rownames = T, show_colnames = T, 
         fontsize=8, legend=TRUE,
         headerplot = "GSEA_cat4_Transcriptomics_Analysis",
         main = "GSEA SubtiWiki Category 4 \n(padj) - Transcriptomics",
         annotation_row = annotation.cat4.cat1
)


