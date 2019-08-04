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
library(edgeR)

####################################################################
### 2. Set Working Directory , Set Parameters 
####################################################################

setwd("~/OneDrive - University of Warwick/WORK/RESULTS/TRANSCRIPTOMICS/intersection_toRNAdo_all")

thr.fc=2 #thresold_FoldChange
thr.pv =0.05 #thresold_pValue1
swati.color = c(brewer.pal(8,"Dark2")) #TA
swati.color2 = c(brewer.pal(9, "Greens")) #SM
swati.color.Blues = c(brewer.pal(8,"Blues"))

#To install this package, start R (version "3.5") and enter:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR", version = "3.8")



#The DGEList data class
#x <- read.delim("TableOfCounts.txt",row.names="Symbol")
x <- read.csv("./filt_txt_all_2/counts.all_ED.csv", header = T)
rownames(x) = x$X
group <- factor(c(1,1,1,
                  2,2,2,
                  3,3,3,
                  4,4,4,
                  5,5,5,
                  6,6,6,
                  7,7,7,
                  8,8,8,
                  9,9,9))
y <- DGEList(counts=x, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)


#To perform quasi-likelihood F-tests:


fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)


#To perform likelihood ratio tests:
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

#Filtering
y$samples


#We filter out lowly expressed genes using the following commands
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]






