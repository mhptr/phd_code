
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

#Input data

#Reading SampleTable
sampleTable = read.csv("./sampleTable_SM.csv", header = T)
sampleTable_WT_SH5_vs_SH2 = sampleTable[c(10:12, 16:18),]


#Reading Combined Read Counts 
counts.all = read.csv("./counts.all_SM_S1583.csv", header = T) #change it 
rownames(counts.all) = counts.all$X
counts.all = counts.all[,-1]
colnames(counts.all) = gsub(".txt","",colnames(counts.all))
colnames(counts.all) = sampleTable$sampleName
counts.all_WT_SH5_vs_SH2 = counts.all[,c(10:12, 16:18)]

################################################################
#5. DESeq2 - Count matrix input
################################################################

cts.SH5_vs_SH2 = counts.all_WT_SH5_vs_SH2
##  Note: 
#In order to benefit from the default settings of the package, 
#you should put the variable of interest at the end of the formula and make sure the control level is the first level.
colData.SH5_vs_SH2 = sampleTable_WT_SH5_vs_SH2   #Match order with column order of cts.M9


#It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. 
rownames(colData.SH5_vs_SH2) = colData.SH5_vs_SH2$sampleName

all(rownames(colData.SH5_vs_SH2) %in% colnames(cts.SH5_vs_SH2)) # should be TRUE
all(rownames(colData.SH5_vs_SH2) == colnames(cts.SH5_vs_SH2)) # should be FALSE  #is TRUE
cts.SH5_vs_SH2 <- cts.SH5_vs_SH2[, rownames(colData.SH5_vs_SH2)]
all(rownames(colData.SH5_vs_SH2) == colnames(cts.SH5_vs_SH2)) # should be TRUE

library("DESeq2")
dds.SH5_vs_SH2 <- DESeqDataSetFromMatrix(countData = cts.SH5_vs_SH2,
                                 colData = colData.SH5_vs_SH2,
                                 design = ~media)  #~1 is the factore suggesting the genotype
dds.SH5_vs_SH2


# class: DESeqDataSet 
# dim: 6495 6 
# metadata(1): version
# assays(1): counts
# rownames(6495): new_4215473_4215670 new_1_148 ... new_4215437_4215505_c S1583
# rowData names(0):
#   colnames(6): S_SH2_W1 S_SH2_W2 ... S_SH5_W2 S_SH5_W3
# colData names(4): sampleName fileName genotype media



################################################################
#5. DESeq2 - Pre-filtering
################################################################

keep.SH5_vs_SH2 = rowSums(counts(dds.SH5_vs_SH2)) >= 20
dds.SH5_vs_SH2 = dds.SH5_vs_SH2[keep.SH5_vs_SH2,]

#### Define factor levels

# Method 1 - Using Factor
#dds.M9$genotype <- factor(dds.M9$genotype, levels = c("wt","ko"))

### OR ###

# Method 2 - Using relevel
dds.SH5_vs_SH2$genotype = relevel(dds.SH5_vs_SH2$media, ref = "SH2")


################################################################
#5. DESeq2 - Differential expression analysis
################################################################

dds.SH5_vs_SH2 <- DESeq(dds.SH5_vs_SH2)      # DESeq is DONE at this step!! 
res.SH5_vs_SH2 <- results(dds.SH5_vs_SH2)
res.SH5_vs_SH2                    ## WORKS! ko vs wt = i.e., condition treated vs untreated

#log2 fold change (MLE): genotype ko vs wt 
#Wald test p-value: genotype ko vs wt 
#DataFrame with 6094 rows and 6 columns


# log2 fold change (MLE): media SH5 vs SH2 
# Wald test p-value: media SH5 vs SH2 
# DataFrame with 6091 rows and 6 columns
# baseMean      log2FoldChange             lfcSE               stat
# <numeric>           <numeric>         <numeric>          <numeric>
#   new_4215473_4215670   2.63399947090185    3.99446641579004  1.69669283645403   2.35426609340686
# new_1_148             4.34689592603584    2.72783781924969  1.41348593738288   1.92986555232404
# new_24_297_c          15.6399107216104 -0.0782309907210783 0.699777257029826 -0.111794131539979
# new_150_409           11.7372250055052   -0.48429691298007 0.728995286486544 -0.664334765886047
# new_299_1067_c         112.90901996405   0.997838544284783 0.325500270852303   3.06555365275735
# ...                                ...                 ...               ...                ...
# BSU41050              76.5967849569196  -0.818414773277247 0.355238571094754  -2.30384547138308
# new_4215104_4215254_c 88.1776029516256   -0.59987832406736 0.357870536901296  -1.67624395475957
# BSU41060              99.5540651925613  -0.614397700926338 0.390098769310725  -1.57497984936464
# new_4215437_4215505_c  12.262284032787   0.975247853604324 0.725189510506607   1.34481792617633
# S1583                 4.58068175012327    2.19434163041609  1.22552841587873   1.79052692861692



################################################################
# Log fold change shrinkage for visualization and ranking
################################################################

resultsNames(dds.SH5_vs_SH2)
resLFC.SH5_vs_SH2 <- lfcShrink(dds.SH5_vs_SH2, coef="media_SH5_vs_SH2", type="apeglm") #apeglm method for effect size shrinkage
resLFC.SH5_vs_SH2


# log2 fold change (MAP): media SH5 vs SH2 
# Wald test p-value: media SH5 vs SH2 
# DataFrame with 6091 rows and 5 columns
# baseMean      log2FoldChange             lfcSE             pvalue
# <numeric>           <numeric>         <numeric>          <numeric>
#   new_4215473_4215670   2.63399947090185     2.3831196788006  2.76249083903094 0.0185593201155135
# new_1_148             4.34689592603584    1.47418234070042  1.31661882048699 0.0536234986346045
# new_24_297_c          15.6399107216104 -0.0518874420909429 0.574445301639903  0.910986640707937
# new_150_409           11.7372250055052  -0.324656356340431 0.607461020210071  0.506476074503922
# new_299_1067_c         112.90901996405    0.93159860487379 0.319767953493664 0.0021726739225278
# ...                                ...                 ...               ...                ...
# BSU41050              76.5967849569196  -0.746015429573274 0.346174309240035 0.0212313195643226
# new_4215104_4215254_c 88.1776029516256  -0.540274678814787 0.345071637322423  0.093690417777404
# BSU41060              99.5540651925613  -0.543433098507143 0.373167539943496  0.115261095857352
# new_4215437_4215505_c  12.262284032787   0.693785214139072 0.644563176802308  0.178684039799958
# S1583                 4.58068175012327    1.24979364995238  1.09247541723121 0.0733692421783611
# padj
# <numeric>
#   new_4215473_4215670   0.0329961526046681
# new_1_148             0.0854580665053312
# new_24_297_c            0.93383029763582
# new_150_409            0.588160092076392
# new_299_1067_c        0.0046191123427982
# ...                                  ...
# BSU41050              0.0372572651876372
# new_4215104_4215254_c  0.141010213660037
# BSU41060               0.168277884675726
# new_4215437_4215505_c  0.245735941842751
# S1583                  0.113395598606546
# 





################################################################
# p-values and adjusted p-values
################################################################


resOrdered.SH5_vs_SH2 <- res.SH5_vs_SH2[order(res.SH5_vs_SH2$padj),]
summary(res.SH5_vs_SH2)
# out of 6091 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2119, 35%
# LFC < 0 (down)     : 1771, 29%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
sum(res.SH5_vs_SH2$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?
#[1] 3890
res05 <- results(dds.SH5_vs_SH2, alpha=0.05) ## adjusted p-value < 0.05
summary(res05)
# out of 6091 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1970, 32%
# LFC < 0 (down)     : 1608, 26%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
sum(res05$padj < 0.05, na.rm=TRUE)
#[1] 3578



################################################################
# Exploring and exporting results - MA-plot
################################################################

##Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res.SH5_vs_SH2, ylim=c(-2,2)) #Points will be colored red if the adjusted p value is less than 0.1
##It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
plotMA(resLFC.SH5_vs_SH2, ylim=c(-2,2))
#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
#idx.M9 <- identify(res.M9$baseMean, res.M9$log2FoldChange)
#rownames(res.M9)[idx.M9]  # Doesn't end



################################################################
# Plot counts
################################################################

plotCounts(dds.SH5_vs_SH2, gene=which.min(res.SH5_vs_SH2$padj), intgroup="genotype")
d.M9 <- plotCounts(dds.SH5_vs_SH2, gene=which.min(res.SH5_vs_SH2$padj), intgroup="genotype", 
                   returnData=TRUE)
library("ggplot2")
ggplot(d.SH5_vs_SH2, aes(x=genotype, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))



################################################################
# More information on results columns
################################################################

mcols(res.SH5_vs_SH2)$description

# [1] "mean of normalized counts for all samples" "log2 fold change (MLE): media SH5 vs SH2" 
# [3] "standard error: media SH5 vs SH2"          "Wald statistic: media SH5 vs SH2"         
# [5] "Wald test p-value: media SH5 vs SH2"       "BH adjusted p-values"  


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

resSig.SH5_vs_SH2 = subset(resOrdered.SH5_vs_SH2)
#write.csv(resSig.SH5_vs_SH2, "../diff.exp.gene/DEG_SM/Diff_all_sm/resSig.SH5_vs_SH2.csv")

resSig.SH5_vs_SH2.05 <- subset(resOrdered.SH5_vs_SH2, padj < 0.05)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
#write.csv(resSig.SH5_vs_SH2.05, "../diff.exp.gene/DEG_SM/Diff_SIG_sm/resSig.SH5_vs_SH2.05.csv")
#write.csv(as.data.frame(resSig.M9.05), 
#          file="DEG_ko_vs_wt_p0.05_M9_sm.csv") #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.SH5_vs_SH2.01 <- subset(resOrdered.SH5_vs_SH2, padj < 0.01)  #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig.SH5_vs_SH2.01
#write.csv(resSig.SH5_vs_SH2.01, "../diff.exp.gene/DEG_SM/Diff_SIG_sm/resSig.SH5_vs_SH2.01.csv")

# log2 fold change (MLE): media SH5 vs SH2 
# Wald test p-value: media SH5 vs SH2 
# DataFrame with 3075 rows and 6 columns
# baseMean    log2FoldChange             lfcSE              stat
# <numeric>         <numeric>         <numeric>         <numeric>
#   BSU29570              55313.1943015603  10.0336573798977 0.273819063947076  36.6433849976093
# BSU08660              27651.3327200329  9.16519862039338 0.301266516567704  30.4222278825124
# BSU09780              4211.55972623408  6.97714105810241 0.242974174584459  28.7155664590071
# BSU16880              5400.01264552632  8.81902229789076 0.309321696537814  28.5108429075639
# BSU09230              6058.71864608124  8.97553190956399  0.31764439112116  28.2565414672801
# ...                                ...               ...               ...               ...
# BSU24660              142.320761430511 0.893457124339648 0.318254087828763  2.80737045809879
# BSU14290              65.1858906956902 -1.03096648020442 0.367450235135086 -2.80573090346617
# BSU_misc_RNA_49       201.145939832828  1.03276047606431 0.368136712772258   2.8053721354958
# new_506502_506670_c   12.4536678935624  2.43028667544768 0.866555287062776   2.8045373581243
# new_1262876_1263611_c 37.4176998690318 -1.16619022325901 0.415818212169366  -2.8045674507013
# pvalue                  padj
# <numeric>             <numeric>
#   BSU29570               5.8323363118723e-294 3.55247604756142e-290
# BSU08660              2.79199601024911e-203 8.50302384921366e-200
# BSU09780              2.43898875641538e-181 4.95196017177535e-178
# BSU16880              8.59643622193649e-179 1.30902232569538e-175
# BSU09230              1.18288755657779e-175 1.44099362142307e-172
# 



# ##########################################################
# ## VENN
# ##########################################################
# 
# 
# resSig.M9.05 = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_M9_minus_spoVG_upp.csv", header = T)
# rownames(resSig.M9.05) = resSig.M9.05$X
# resSig.SH2.05 = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_SH2_minus_spoVG_upp.csv", header = T)
# rownames(resSig.SH2.05) = resSig.SH2.05$X
# resSig.SH5.05 = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_SH5_minus_spoVG_upp.csv", header = T)
# rownames(resSig.SH5.05) = resSig.SH5.05$X
# 
# venn.plot = venn.diagram(
#   list("M9_Transcriptomic" = rownames(resSig.M9.05),
#        "SH2_Transcriptomic" = rownames(resSig.SH2.05),
#        "SH5_Transcriptomic" = rownames(resSig.SH5.05)),
#   filename = NULL,
#   col="transparent",
#   fill=swati.color[1:3], height = 3000, width = 3000,resolution =500, imagetype = "tiff",
#   main = "Venn Diagram - All Conditions \nTranscriptomic Analysis",main.pos  = c(0.5, 1.05), main.fontface = "bold",
#   main.fontfamily = "serif", main.col = "Dark Blue", main.cex = 1.5
# )
# grid.newpage()
# grid.draw(venn.plot)
# 
# 
# 
# #Getting intersects of overlaps
# #https://stackoverflow.com/questions/23559371/how-to-get-the-list-of-items-in-venn-diagram-in-r
# 
# library(reshape)
# library(R.utils)
# 
# ## data
# M9_Transcriptomic = data.frame("M9_Transcriptomic" = rownames(resSig.M9.05), "M9_Transcriptomic" = 1)
# SH2_Transcriptomic = data.frame("SH2_Transcriptomic" = rownames(resSig.SH2.05), "SH2_Transcriptomic" = 1)
# SH5_Transcriptomic = data.frame("SH5_Transcriptomic" = rownames(resSig.SH5.05), "SH5_Transcriptomic" = 1)
# 
# 
# ## a merged data frame.
# Transcripts <- list(M9_Transcriptomic = M9_Transcriptomic, 
#                     SH2_Transcriptomic = SH2_Transcriptomic,
#                     SH5_Transcriptomic = SH5_Transcriptomic)
# merged_Transcripts <- merge_recurse(Transcripts)
# 
# 
# 
# 
# ##################################################################################################
# ###### Gene Set Enrichment Analysis - Category 3 (GSEA-3) with padj
# ##################################################################################################
# 
# 
# #source("http://www.bioconductor.org/biocLite.R")
# #BiocManager::install("piano", dependencies=TRUE)          #this worked
# DEG_ko_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_vs_wt_all_M9_sm.csv", header = T)
# dim(DEG_ko_M9)
# DEG_ko_M9 = DEG_ko_M9[!apply(DEG_ko_M9, 1, function(x) any(is.na(x))), ]
# dim(DEG_ko_M9)
# #write.csv(DEG_ko_M9, "../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_M9_sm.csv")
# DEG_ko_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_vs_wt_all_SH2_sm.csv", header = T)
# dim(DEG_ko_SH2)
# DEG_ko_SH2 = DEG_ko_SH2[!apply(DEG_ko_SH2, 1, function(x) any(is.na(x))), ]
# dim(DEG_ko_SH2)
# DEG_ko_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_vs_wt_all_SH5_sm.csv", header = T)
# dim(DEG_ko_SH5)
# DEG_ko_SH5 = DEG_ko_SH5[!apply(DEG_ko_SH5, 1, function(x) any(is.na(x))), ]
# dim(DEG_ko_SH5)
# 
# 
# 
# DEG_all = merge(DEG_ko_M9, DEG_ko_SH2, by.x="X", by.y="X", all.x = TRUE)
# DEG_all = merge(DEG_all, DEG_ko_SH5, by.x="X", by.y="X", all.x = TRUE)
# dim(DEG_all)
# DEG_all = DEG_all[!apply(DEG_all, 1, function(x) any(is.na(x))), ]
# dim(DEG_all)
# 
# 
# colnames(DEG_all) = c("X", "baseMean_M9", "log2FoldChange_M9", "lfcSE_M9", "stat_M9", "pvalue_M9", "padj_M9",
#                       "baseMean_SH2", "log2FoldChange_SH2", "lfcSE_SH2", "stat_SH2", "pvalue_SH2", "padj_SH2",
#                       "baseMean_SH5", "log2FoldChange_SH5", "lfcSE_SH5", "stat_SH5", "pvalue_SH5", "padj_SH5")
# rownames(DEG_all) = DEG_all$X
# DEG_all = DEG_all[,-1]
# 
# 
# DEG_ko_M9_clean = DEG_all[,1:6]
# DEG_ko_SH2_clean = DEG_all[,7:12]
# DEG_ko_SH5_clean = DEG_all[,13:18]
# 
# 
# geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = TRUE)
# myGsc.cat3 = loadGSC(cbind(as.character(geneCategories$gene),
#                            as.character(geneCategories$category3)))
# DEG_all_cond = list(DEG_ko_M9_clean, DEG_ko_SH2_clean, DEG_ko_SH5_clean)
# names(DEG_all_cond) = c("M9","SH2","SH5")
# 
# gsea.enrich.cat3.all = c()
# for (i in 1:length(DEG_all_cond)){
#   my.de = DEG_all_cond[[i]]
#   my.pval = my.de$padj
#   my.fc = my.de$log2FoldChange
#   names(my.pval) = rownames(my.de)
#   names(my.fc) = rownames(my.de)
#   my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc
#   my.fc = my.fc[!is.na(my.fc)]    #1707 entries 
#   my.gsaRes.cat3 = runGSA(geneLevelStats=my.pval,directions=my.fc,
#                           gsc=myGsc.cat3,geneSetStat="reporter",adjMethod="BH")
#   my.gsea.cat3 = GSAsummaryTable(my.gsaRes.cat3)
#   colnames(my.gsea.cat3) = paste0(names(DEG_all_cond)[i],"_",colnames(my.gsea.cat3))
#   my.gsea.cat3 = as.matrix(my.gsea.cat3)
#   gsea.enrich.cat3.all=cbind(gsea.enrich.cat3.all,my.gsea.cat3)
# }
# 
# #write.csv(gsea.enrich.cat3.all, "../diff.exp.gene/DEG_SM/SGC_sm/gsea.enrich.cat3.all_sm.csv")
# 
# dim(gsea.enrich.cat3.all)
# #[1] 121  57
# 
# #gsea.go.sig = gsea.go[gsea.go$`p (non-dir.)` < 0.05,]
# rownames(gsea.enrich.cat3.all) = gsea.enrich.cat3.all[,"M9_Name"]
# #write.csv(gsea.enrich.cat3.all, "../diff.exp.gene/DEG_SM/SGC_sm/gsea.enrich.cat3.all_rowname_sm.csv")
# 
# #gsea.enrich.pv =  gsea.enrich.all[,grep("LB_p(non-dir.)", colnames(gsea.enrich.all))]
# gsea.cat3.enrich.ND =  gsea.enrich.cat3.all[,grep("(non-dir.)", colnames(gsea.enrich.cat3.all))]
# colnames(gsea.cat3.enrich.ND) = c("M9_Stat(non-dir.)",
#                                   "M9_p(non-dir.)",
#                                   "M9_padj(non-dir.)",
#                                   "SH2_Stat(non-dir.)",
#                                   "SH2_p(non-dir.)",
#                                   "SH2_padj(non-dir.)",
#                                   "SH5_Stat(non-dir.)",
#                                   "SH5_p(non-dir.)",
#                                   "SH5_padj(non-dir.)")
# 
# 
# colnames(gsea.cat3.enrich.ND) = c("M9_Stat_ND",
#                                   "M9_p_ND",
#                                   "M9_padj_ND",
#                                   "SH2_Stat_ND",
#                                   "SH2_p_ND",
#                                   "SH2_padj_ND",
#                                   "SH5_Stat_ND",
#                                   "SH5_p_ND",
#                                   "SH5_padj_ND")
# 
# gsea.cat3.enrich = matrix(as.numeric(unlist(gsea.cat3.enrich.ND)), nrow = nrow(gsea.cat3.enrich.ND))
# rownames(gsea.cat3.enrich) = rownames(gsea.cat3.enrich.ND)
# colnames(gsea.cat3.enrich) = colnames(gsea.cat3.enrich.ND)
# 
# gsea.cat3.enrich.sig = gsea.cat3.enrich[(as.numeric(as.character(gsea.cat3.enrich[,"M9_padj_ND"]))<0.05|
#                                            as.numeric(as.character(gsea.cat3.enrich[,"SH2_padj_ND"]))<0.05|
#                                            as.numeric(as.character(gsea.cat3.enrich[,"SH5_padj_ND"]))<0.05),]
# 
# #write.csv(gsea.cat3.enrich.sig, "../diff.exp.gene/DEG_SM/SGC_sm/gsea.cat3.enrich.sig_sm.csv")
# 
# gsea.cat3.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat3.enrich.sig[,"M9_padj_ND"])),
#                                 as.numeric(as.character(gsea.cat3.enrich.sig[,"SH2_padj_ND"])),
#                                 as.numeric(as.character(gsea.cat3.enrich.sig[,"SH5_padj_ND"])))
# rownames(gsea.cat3.enrich.sig.pv) = rownames(gsea.cat3.enrich.sig)
# colnames(gsea.cat3.enrich.sig.pv) = c("GSEA_Cat3_M9_padj","GSEA_Cat3_SH2_padj","GSEA_Cat3_SH5_padj")
# gsea.cat3.enrich.sig.pv.log10 = -log10(gsea.cat3.enrich.sig.pv)
# gsea.cat3.enrich.sig.pv.log10[gsea.cat3.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!
# 
# #write.csv(gsea.cat3.enrich.sig.pv,
# #          "../diff.exp.gene/DEG_SM/SGC_sm/gsea.cat3.enrich.sig.pv_sm.csv")
# 
# ####      !!!   ANNOTATION need FIXING     !!! 
# 
# #annotation.cat3.cat1 = data.frame(cat1 = unlist(lapply(rownames(gsea.cat3.enrich.sig.pv.log10), function(x)
# #  geneCategories$category1[geneCategories$category3==x][1])))
# #rownames(annotation.cat3.cat1) = rownames(gsea.cat3.enrich.sig.pv.log10)
# #write.csv(annotation.cat3.cat1, "../diff.exp.gene/DEG_SM/SGC_sm/annotation.cat3.cat1.csv")
# 
# pheatmap(gsea.cat3.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
#          cluster_rows=TRUE, cluster_cols=FALSE,
#          cellwidth = 35, cellheight = 25,
#          show_rownames = T, show_colnames = T, 
#          fontsize=8, legend=TRUE,
#          headerplot = "GSEA_Cat3_Transcriptomics_Analysis",
#          main = "GSEA SubtiWiki Category 3 \n(padj) - Transcriptomics"
#          #annotation_row = annotation.cat3.cat1
# )
# 
# 
# 
# ##################################################################################################
# ###### Gene Set Enrichment Analysis - Category 1 (GSEA-1) with padj
# ##################################################################################################
# 
# 
# #source("http://www.bioconductor.org/biocLite.R")
# #BiocManager::install("piano", dependencies=TRUE)          #this worked
# DEG_ko_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_vs_wt_all_M9_sm.csv", header = T)
# dim(DEG_ko_M9)
# DEG_ko_M9 = DEG_ko_M9[!apply(DEG_ko_M9, 1, function(x) any(is.na(x))), ]
# dim(DEG_ko_M9)
# #write.csv(DEG_ko_M9, "../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_M9_sm.csv")
# DEG_ko_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_vs_wt_all_SH2_sm.csv", header = T)
# dim(DEG_ko_SH2)
# DEG_ko_SH2 = DEG_ko_SH2[!apply(DEG_ko_SH2, 1, function(x) any(is.na(x))), ]
# dim(DEG_ko_SH2)
# DEG_ko_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/DEG_ko_vs_wt_all_SH5_sm.csv", header = T)
# dim(DEG_ko_SH5)
# DEG_ko_SH5 = DEG_ko_SH5[!apply(DEG_ko_SH5, 1, function(x) any(is.na(x))), ]
# dim(DEG_ko_SH5)
# 
# 
# DEG_all = merge(DEG_ko_M9, DEG_ko_SH2, by.x="X", by.y="X", all.x = TRUE)
# DEG_all = merge(DEG_all, DEG_ko_SH5, by.x="X", by.y="X", all.x = TRUE)
# dim(DEG_all)
# DEG_all = DEG_all[!apply(DEG_all, 1, function(x) any(is.na(x))), ]
# dim(DEG_all)
# 
# 
# colnames(DEG_all) = c("X", "baseMean_M9", "log2FoldChange_M9", "lfcSE_M9", "stat_M9", "pvalue_M9", "padj_M9",
#                       "baseMean_SH2", "log2FoldChange_SH2", "lfcSE_SH2", "stat_SH2", "pvalue_SH2", "padj_SH2",
#                       "baseMean_SH5", "log2FoldChange_SH5", "lfcSE_SH5", "stat_SH5", "pvalue_SH5", "padj_SH5")
# rownames(DEG_all) = DEG_all$X
# DEG_all = DEG_all[,-1]
# 
# 
# DEG_ko_M9_clean = DEG_all[,1:6]
# DEG_ko_SH2_clean = DEG_all[,7:12]
# DEG_ko_SH5_clean = DEG_all[,13:18]
# 
# 
# geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = TRUE)
# myGsc.cat1 = loadGSC(cbind(as.character(geneCategories$gene),
#                            as.character(geneCategories$category1)))
# DEG_all_cond = list(DEG_ko_M9_clean, DEG_ko_SH2_clean, DEG_ko_SH5_clean)
# names(DEG_all_cond) = c("M9","SH2","SH5")
# 
# gsea.enrich.cat1.all = c()
# for (i in 1:length(DEG_all_cond)){
#   my.de = DEG_all_cond[[i]]
#   my.pval = my.de$padj
#   my.fc = my.de$log2FoldChange
#   names(my.pval) = rownames(my.de)
#   names(my.fc) = rownames(my.de)
#   my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc
#   my.fc = my.fc[!is.na(my.fc)]    #1707 entries 
#   my.gsaRes.cat1 = runGSA(geneLevelStats=my.pval,directions=my.fc,
#                           gsc=myGsc.cat1,geneSetStat="reporter",adjMethod="BH")
#   my.gsea.cat1 = GSAsummaryTable(my.gsaRes.cat1)
#   colnames(my.gsea.cat1) = paste0(names(DEG_all_cond)[i],"_",colnames(my.gsea.cat1))
#   my.gsea.cat1 = as.matrix(my.gsea.cat1)
#   gsea.enrich.cat1.all=cbind(gsea.enrich.cat1.all,my.gsea.cat1)
# }
# 
# #write.csv(gsea.enrich.cat1.all, "../diff.exp.gene/DEG_SM/SGC_sm/gsea.enrich.cat1.all_sm.csv")
# 
# dim(gsea.enrich.cat1.all)
# #[1]  6 57
# 
# #gsea.go.sig = gsea.go[gsea.go$`p (non-dir.)` < 0.05,]
# rownames(gsea.enrich.cat1.all) = gsea.enrich.cat1.all[,"M9_Name"]
# #write.csv(gsea.enrich.cat1.all, "../diff.exp.gene/DEG_SM/SGC_sm/gsea.enrich.cat1.all_rowname_sm.csv")
# 
# #gsea.enrich.pv =  gsea.enrich.all[,grep("LB_p(non-dir.)", colnames(gsea.enrich.all))]
# gsea.cat1.enrich.ND =  gsea.enrich.cat1.all[,grep("(non-dir.)", colnames(gsea.enrich.cat1.all))]
# colnames(gsea.cat1.enrich.ND) = c("M9_Stat(non-dir.)",
#                                   "M9_p(non-dir.)",
#                                   "M9_padj(non-dir.)",
#                                   "SH2_Stat(non-dir.)",
#                                   "SH2_p(non-dir.)",
#                                   "SH2_padj(non-dir.)",
#                                   "SH5_Stat(non-dir.)",
#                                   "SH5_p(non-dir.)",
#                                   "SH5_padj(non-dir.)")
# 
# 
# colnames(gsea.cat1.enrich.ND) = c("M9_Stat_ND",
#                                   "M9_p_ND",
#                                   "M9_padj_ND",
#                                   "SH2_Stat_ND",
#                                   "SH2_p_ND",
#                                   "SH2_padj_ND",
#                                   "SH5_Stat_ND",
#                                   "SH5_p_ND",
#                                   "SH5_padj_ND")
# 
# gsea.cat1.enrich = matrix(as.numeric(unlist(gsea.cat1.enrich.ND)), nrow = nrow(gsea.cat1.enrich.ND))
# rownames(gsea.cat1.enrich) = rownames(gsea.cat1.enrich.ND)
# colnames(gsea.cat1.enrich) = colnames(gsea.cat1.enrich.ND)
# 
# gsea.cat1.enrich.sig = gsea.cat1.enrich[(as.numeric(as.character(gsea.cat1.enrich[,"M9_padj_ND"]))<0.05|
#                                            as.numeric(as.character(gsea.cat1.enrich[,"SH2_padj_ND"]))<0.05|
#                                            as.numeric(as.character(gsea.cat1.enrich[,"SH5_padj_ND"]))<0.05),]
# 
# #write.csv(gsea.cat1.enrich.sig, "../diff.exp.gene/DEG_SM/SGC_sm/gsea.cat1.enrich.sig_sm.csv")
# 
# gsea.cat1.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat1.enrich.sig[,"M9_padj_ND"])),
#                                 as.numeric(as.character(gsea.cat1.enrich.sig[,"SH2_padj_ND"])),
#                                 as.numeric(as.character(gsea.cat1.enrich.sig[,"SH5_padj_ND"])))
# rownames(gsea.cat1.enrich.sig.pv) = rownames(gsea.cat1.enrich.sig)
# colnames(gsea.cat1.enrich.sig.pv) = c("GSEA_Cat1_M9_padj","GSEA_Cat1_SH2_padj","GSEA_Cat1_SH5_padj")
# gsea.cat1.enrich.sig.pv.log10 = -log10(gsea.cat1.enrich.sig.pv)
# gsea.cat1.enrich.sig.pv.log10[gsea.cat1.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.1 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!
# 
# #write.csv(gsea.cat1.enrich.sig.pv,
# #          "../diff.exp.gene/DEG_SM/SGC_sm/gsea.cat1.enrich.sig.pv_sm.csv")
# 
# 
# pheatmap(gsea.cat1.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
#          cluster_rows=TRUE, cluster_cols=FALSE,
#          cellwidth = 50, cellheight = 30,
#          show_rownames = T, show_colnames = T, 
#          fontsize=8, legend=TRUE,
#          headerplot = "GSEA_Cat1_Transcriptomics_Analysis",
#          main = "GSEA SubtiWiki Category 1 \n(padj) - Transcriptomics"
# )
# 
# ###################################  ALL GOOD SO FAR  ################################### 
# 
# 
# #14/07/19
# 
# #################################################################
# ###                   Volcano Plots
# ##################################################################
# 
# setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/filt_txt_all_2/")
# library(ggplot2)
# # volcano_trans = function(DEG, tag) {
# #   DEG_File_Location = paste("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/", DEG, sep = "")
# #   DEG_File = read.csv(DEG_File_Location, header = T)
# #   rownames(DEG_File) = DEG_File$X
# #   DEG_File = DEG_File[,-1]
# #   
# #   adj.P.Val =0.05
# #   logFC=0.05
# #   g = ggplot(data=DEG_File,
# #     aes(title, x=logFC, y=-log10(adj.P.Val), colour=adj.P.Val<0.05)) +   
# #     #colour=-log10(P.Value) +
# #     ggtitle("Volcano Plot - Transcriptomics Analysis") +      #CHANGE THE NAME HERE 
# #     theme(
# #       plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
# #       axis.title.x = element_text(color="black", size=12, face="bold"),
# #       axis.title.y = element_text(color="black", size=12, face="bold") 
# #     ) + 
# #     geom_hline(yintercept=-log10(adj.P.Val), linetype="dashed", color = "red") + #thr.fc is 2
# #         geom_point(alpha=0.4, size=1.75) +
# #     ylim(c(0,23))+  xlim(c(-12,12))+   #LB, M9, SH2
# #     xlab("log2 Fold Change (Del_spoVG vs WT)") + ylab("-log10 adj.P.Val")
# #   print(g)
# # }
# # volcano_trans("DEG_ko_vs_wt_p0.05_M9_minus_spoVG_upp.csv", "M9")
# 
# 
# DEG = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH2_minus_spoVG_upp.csv", header = T)
# rownames(DEG) = DEG$X
# DEG = DEG[,-1]
# thr.adj.pv =0.05
# g = ggplot(data=DEG, 
#     aes(title, x=log2FoldChange, y=-log10(padj), colour=padj<0.05)) +
#   ggtitle("Volcano Plot - Transcriptomics Analysis") +      #CHANGE THE NAME HERE 
#   theme(
#     plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
#     axis.title.x = element_text(color="black", size=12, face="bold"),
#     axis.title.y = element_text(color="black", size=12, face="bold") 
#   ) + 
#   geom_hline(yintercept=-log10(thr.adj.pv), linetype="dashed", color = "red") +
#   #geom_vline(xintercept=-log2(thr.fc), linetype="dashed", color = "red") +  #thr.fc is 2
#   #geom_vline(xintercept=log2(thr.fc), linetype="dashed", color = "red") +
#   geom_point(alpha=0.4, size=1.75) +
#   #ylim(c(0,23))+  xlim(c(-12,12))+   #LB, M9, SH2
#   ylim(c(0,100))+  xlim(c(-0.25,0.25))+   #SH5
#   xlab("log2 Fold Change (Del_spoVG vs WT)") + ylab("-log10 adj.P.Val")
# print(g)
# 
# 
# 
# #https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html#introduction
# # if (!requireNamespace('BiocManager', quietly = TRUE))
# #   install.packages('BiocManager')
# # BiocManager::install('EnhancedVolcano')
# #devtools::install_github('kevinblighe/EnhancedVolcano')
# library(EnhancedVolcano)
# 
# EnhancedVolcano(DEG,
#                 lab = rownames(DEG),
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 xlim = c(-8, 8),
#                 title = 'Volcano Plot - Transcriptomics Analysis',
#                 pCutoff = 10e-16,
#                 FCcutoff = 1.5,
#                 transcriptPointSize = 1.5,
#                 transcriptLabSize = 3.0)
# 
# 
# EnhancedVolcano(DEG,
#                 lab = rownames(DEG),
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 xlim = c(-6, 6),
#                 ylim = c(0,77),
#                 title = 'Volcano Plot - Transcriptomics Analysis',
#                 pCutoff = 10e-16,
#                 FCcutoff = 1.5,
#                 transcriptPointSize = 1.5,
#                 transcriptLabSize = 3.0,
#                 col=c('black', 'orange', 'green', 'red3'),
#                 colAlpha = 1)
# 
# 
# 
# EnhancedVolcano(DEG,
#                 lab = rownames(DEG),
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 xlim = c(-6, 6),
#                 pCutoff = 10e-12,
#                 FCcutoff = 1.5,
#                 cutoffLineType = 'twodash',
#                 cutoffLineWidth = 0.8,
#                 transcriptPointSize = 3.0,
#                 transcriptLabSize = 4.0,
#                 colAlpha = 1,
#                 legend=c('NS','Log (base 2) fold-change','P value',
#                          'P value & Log (base 2) fold-change'),
#                 legendPosition = 'right',
#                 legendLabSize = 16,
#                 legendIconSize = 5.0)
# 
# 
# 
# 
# 
# DEG = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_M9_minus_spoVG_upp.csv", header = T)
# rownames(DEG) = DEG$X
# DEG = DEG[,-1]
# 
# 
# EnhancedVolcano(DEG,
#                 lab = rownames(DEG),
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 selectLab = c('new_3788276_3788397_c', 'BSU17180', 'BSU00500'),
#                 xlim = c(-5.5,8),
#                 ylim = c(0, 25),
#                 xlab = bquote(~Log[2]~ 'fold change'),
#                 pCutoff = 10e-14,
#                 FCcutoff = 2.0,
#                 transcriptPointSize = 3.0,
#                 transcriptLabSize = 5.0,
#                 transcriptLabCol = 'black',
#                 transcriptLabFace = 'bold',
#                 boxedlabels = TRUE,
#                 colAlpha = 4/5,
#                 legend=c('NS','Log (base 2) fold-change','padj',
#                          'padj & Log (base 2) fold-change'),
#                 legendPosition = 'right',
#                 legendLabSize = 14,
#                 legendIconSize = 4.0,
#                 drawConnectors = TRUE,
#                 widthConnectors = 1.0,
#                 colConnectors = 'black')
# 
# 
# 
# 
# 
# ###################
# 
# #15/07/19
# 
# library(ggplot2)
# library(scales)
# library(limma)
# # # set up an example dataset (see ?lmFit)
# # sd <- 0.3*sqrt(4/rchisq(100,df=4))
# # y <- matrix(rnorm(100*6,sd=sd),100,6)
# # rownames(y) <- paste("Gene",1:100)
# # y[1:2,4:6] <- y[1:2,4:6] + 2
# # design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
# # options(digits=3)
# 
# # Ordinary fit
# 
# dim(DEG$log2FoldChange)
# 
# fit <- lmFit(DEG$log2FoldChange, DEG$padj)
# fit <- eBayes(fit)
# tt <- topTable(fit,coef=2, n = Inf)
# 
# # transformation function for reverse log1p axis
# revlog_trans <- function(base = exp(1)) {
#   trans <- function(x) -log1p(x)
#   inv <- function(x) expm1(-x)
#   scales::trans_new("revlog1p", trans, inv, domain = c(0, Inf))
# }
# 
# ggplot(tt, aes(x = logFC, y = P.Value)) +
#   scale_fill_gradient(low = "lightgray", high = "navy") +
#   scale_color_gradient(low = "lightgray", high = "navy") +
#   scale_y_continuous(trans = revlog_trans(), expand = c(0.005, 0.005)) +
#   expand_limits(y = c(0, 1)) +
#   stat_density_2d(aes(fill = ..level..), geom = "polygon",
#                   show.legend = FALSE) +
#   geom_point(data = subset(tt, P.Value < 0.05), 
#              color = "red", alpha = 0.5) +
#   geom_vline(xintercept = 0) +
#   geom_hline(yintercept = 0) +
#   geom_hline(yintercept = 0.05, linetype = "dashed") +
#   geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
#   theme_linedraw() + 
#   theme(panel.grid = element_blank()) +
#   xlab("Fold change (log2)") +
#   ylab("P-Value") +
#   annotate("text", x = min(tt$logFC), y = 1,
#            label = "Nominally significant",
#            color = "black", hjust = 0) +
#   annotate("point", x = min(tt$logFC) - 0.05, y = 1, color = "red")
# 
# 
# ######################################################################################
# ### 5. Installing KEGGREST  #No need to run - Already downloaded csv exists  
# ######################################################################################
# 
# 
# ## try http:// if https:// URLs are not supported
# #source("https://bioconductor.org/biocLite.R")
# #pkgs <- rownames(installed.packages())
# #biocLite(pkgs, type="source")                     #something worked 
# 
# 
# ## try http:// if https:// URLs are not supported
# #source("https://bioconductor.org/biocLite.R")
# #biocLite("KEGGREST")                              #worked 
# 
# #library(KEGGREST)
# #bsu.geneName = keggList("bsu")
# #bsu.geneName.df = data.frame(geneName = bsu.geneName)
# #bsu.geneName.df$geneID = names(bsu.geneName)
# 
# #bsu.uniprot = keggConv("bsu", "uniprot")
# #bsu.uniprot.df = data.frame(geneID = bsu.uniprot)
# #bsu.uniprot.df$uniProt = names(bsu.uniprot)
# 
# #bsu.pathways.gene = keggLink("pathway", "bsu")
# #bsu.pathways.gene.df = data.frame(keggPathwayID = bsu.pathways.gene)
# #bsu.pathways.gene.df$geneID = names(bsu.pathways.gene)
# 
# #bsu.keggPathName = keggList("pathway", "bsu")
# #bsu.keggPathName.df = data.frame(keggPathName = bsu.keggPathName)
# #bsu.keggPathName.df$keggPathwayID = names(bsu.keggPathName)
# #bsu.keggPathName.df$keggPathName = gsub(" - Bacillus subtilis subsp. subtilis 168","",bsu.keggPathName.df$keggPathName)
# 
# #bsu.kegg.pathway = merge(merge(bsu.geneName.df,bsu.uniprot.df,by="geneID"),bsu.pathways.gene.df,by="geneID")
# #bsu.kegg.pathway$keggPathName = unlist(lapply(bsu.kegg.pathway$keggPathwayID, function(x)
# #  bsu.keggPathName.df$keggPathName[bsu.keggPathName.df$keggPathwayID==x]))
# #bsu.kegg.pathway$geneID = gsub("bsu:","",bsu.kegg.pathway$geneID)
# #bsu.kegg.pathway$uniProt = gsub("up:","",bsu.kegg.pathway$uniProt)
# #bsu.kegg.pathway$keggPathwayID = gsub("path:","",bsu.kegg.pathway$keggPathwayID)
# 
# #write.csv(bsu.kegg.pathway,"./data/bsu.kegg.pathway.csv",row.names = F)
# 
# 
# 
# 
# ######################################################################################
# ### 5. KEGG Pathway Enrichment
# ######################################################################################
# 
# library(pheatmap)
# 
# DEG_all_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_M9) = DEG_all_M9$X
# DEG_all_M9 = DEG_all_M9[,-1]
# 
# DEG_all_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH2_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_SH2) = DEG_all_SH2$X
# DEG_all_SH2 = DEG_all_SH2[,-1]
# 
# DEG_all_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH5_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_SH5) = DEG_all_SH5$X
# DEG_all_SH5 = DEG_all_SH5[,-1]
# 
# bsu.kegg.pathway = read.csv("../bsu.kegg.pathway.csv", header = T)
# DEG_all = list(DEG_all_M9,DEG_all_SH2,DEG_all_SH5)
# names(DEG_all) = c("M9","SH2","SH5")
# kegg.enrich.all = c()                              #make dataframe 
# for(i in 1:length(DEG_all)){
#   my.de = DEG_all[[i]]
#   my.de.genes = rownames(my.de)[my.de$padj<0.05]
#   m = length(my.de.genes)                             #m = length of Sig Diff Ex pr
#   n = dim(my.de)[1] - m                               #n = length of all diff. protein - length of Sig Diff Ex pr
#   my.kegg.enrich = c()   
#   for(j in 1:length(unique(bsu.kegg.pathway$keggPathName))) {
#     my.path = unique(bsu.kegg.pathway$keggPathName)[j]
#     my.path.gene = unique(bsu.kegg.pathway$geneID[bsu.kegg.pathway$keggPathName == my.path])
#     k = length(my.path.gene)             # k = Total no of genes in pathway
#     x = sum(my.de.genes%in%my.path.gene) # x = no of differentially expressed genes in a particular pathway
#     my.de.path.gene = paste(my.de.genes[my.de.genes%in%my.path.gene],collapse = "|")
#     pv = phyper(x-1,m,n,k,lower.tail = F) 
#     my.kegg.enrich = rbind(my.kegg.enrich,cbind(pv,k,x,my.de.path.gene))
#   }
#   colnames(my.kegg.enrich) = paste0(names(DEG_all)[i],"_",colnames(my.kegg.enrich))
#   kegg.enrich.all = as.data.frame(cbind(kegg.enrich.all,my.kegg.enrich)) 
# }
# rownames(kegg.enrich.all) = unique(bsu.kegg.pathway$keggPathName)
# #kegg.enrich.all = kegg.enrich.all[!apply(kegg.enrich.all, 1, function(x) any(is.na(x))),]
# 
# kegg.enrich.all = kegg.enrich.all[as.numeric(as.character(kegg.enrich.all$M9_k))>2,]
# kegg.enrich.sig = kegg.enrich.all[(as.numeric(as.character(kegg.enrich.all$M9_pv))<0.05|
#                                      as.numeric(as.character(kegg.enrich.all$SH2_pv))<0.05|
#                                      as.numeric(as.character(kegg.enrich.all$SH5_pv))<0.05),]
# #write.csv(kegg.enrich.sig, "../diff.exp.gene/DEG_SM/DEG_KEGG/kegg.enrich.sig_transcriptomics.csv")
# 
# kegg.enrich.sig.pv = cbind(as.numeric(as.character(kegg.enrich.sig$M9_pv)),
#                            as.numeric(as.character(kegg.enrich.sig$SH2_pv)),
#                            as.numeric(as.character(kegg.enrich.sig$SH5_pv)))
# rownames(kegg.enrich.sig.pv) = rownames(kegg.enrich.sig)
# colnames(kegg.enrich.sig.pv) = c("KEGG_M9_adj.P.Val","KEGG_SH2_adj.P.Val","KEGG_SH5_adj.P.Val")
# #write.csv(kegg.enrich.sig.pv, "../diff.exp.gene/DEG_SM/DEG_KEGG/kegg.enrich.sig.pv_transcriptomics.csv")
# 
# kegg.enrich.sig.pv.log10 = -log10(kegg.enrich.sig.pv)
# kegg.enrich.sig.pv.log10[kegg.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!
# 
# pheatmap(kegg.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
#          cluster_rows=TRUE, cluster_cols=FALSE,
#          cellwidth = 70, cellheight = 27,
#          show_rownames = T, show_colnames = T, 
#          fontsize=12, legend=TRUE,
#          headerplot = "KEEG_Transcriptomics_Analysis",
#          main = "KEGG Pathway Enrichment Analysis \n(adj.P.Val) - Transcriptomics"
# )
# 
# 
# 
# ##################################################################################################
# ###### Gene Set Enrichment Analysis - Category 3 (GSEA-3) with adj.P.Val
# ##################################################################################################
# 
# 
# #source("http://www.bioconductor.org/biocLite.R")
# #BiocManager::install("piano", dependencies=TRUE)          #this worked
# 
# library(piano)
# #Package piano version 1.22.0 Index
# geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = T)
# myGsc.cat3 = loadGSC(cbind(as.character(geneCategories$gene),
#                            as.character(geneCategories$category3)))
# 
# DEG_all_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_M9) = DEG_all_M9$X
# DEG_all_M9 = DEG_all_M9[,-1]
# 
# DEG_all_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH2_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_SH2) = DEG_all_SH2$X
# DEG_all_SH2 = DEG_all_SH2[,-1]
# 
# DEG_all_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_SH5_minus_spoVG_upp.csv", header = T)
# rownames(DEG_all_SH5) = DEG_all_SH5$X
# DEG_all_SH5 = DEG_all_SH5[,-1]
# 
# DEG_all = list(DEG_all_M9, DEG_all_SH2, DEG_all_SH5)
# names(DEG_all) = c("M9","SH2","SH5")
# 
# # gsea.enrich.cat3.all = c()
# # for (i in 1:length(DEG_all)){
# #   my.de = DEG_all[[i]]
# #   my.pval = my.de$padj
# #   my.fc = my.de$log2FoldChange
# #   names(my.pval) = rownames(my.de)
# #   names(my.fc) = rownames(my.de)
# #   my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc #m9 pval = NA BSU11750
# #   my.fc = my.fc[!is.na(my.fc)]     
# #   my.gsaRes.cat3 = runGSA(geneLevelStats=my.pval,directions=my.fc,
# #                           gsc=myGsc.cat3,geneSetStat="reporter",adjMethod="BH")
# #   my.gsea.cat3 = GSAsummaryTable(my.gsaRes.cat3)
# #   colnames(my.gsea.cat3) = paste0(names(DEG_all)[i],"_",colnames(my.gsea.cat3))
# #   my.gsea.cat3 = as.matrix(my.gsea.cat3)
# #   gsea.enrich.cat3.all=cbind(gsea.enrich.cat3.all,my.gsea.cat3)
# # }
# 
# gsea.enrich.cat3.all = c()
# for (i in 1:length(DEG_all)){
#   my.de = DEG_all[[i]]
#   my.pval = my.de$padj
#   my.fc = my.de$log2FoldChange
#   length(my.pval)
#   length(my.fc)
#   names(my.pval) = rownames(my.de)
#   names(my.fc) = rownames(my.de)
#   my.pval = my.pval[!is.na(my.pval)]
#   my.fc = my.pval[!is.na(my.pval)]
#   length(my.pval)
#   length(my.fc)
#   #m9 pval & fc = 6092, m9 pval - na = 5383, m9 fc - na = 6092
#   #SH2 pval & fc = 5776, SH2 pval - na = 5769, SH2 fc - na = 5776
#   #SH5 pval & fc = 6184, SH5 pval - na = 4864, SH5 fc - na = 6184
#   # my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc #m9 pval = NA BSU11750
#   my.fc = my.fc[!is.na(my.fc)]
#   my.gsaRes.cat3 = runGSA(geneLevelStats=my.pval,directions=my.fc,
#                           gsc=myGsc.cat3,geneSetStat="reporter",adjMethod="BH")
#   # M9 = 3774 genes and 121 gene sets
#   # SH2 = 4006 genes and 122 gene sets
#   # SH5 = 3503 genes and 122 gene sets
#   my.gsea.cat3 = GSAsummaryTable(my.gsaRes.cat3)
#   colnames(my.gsea.cat3) = paste0(names(DEG_all)[i],"_",colnames(my.gsea.cat3))
#   my.gsea.cat3 = as.matrix(my.gsea.cat3)
#   dim(my.gsea.cat3)
#   # M9 = [1] 121  16
#   # SH2 = [1] 122  16
#   # SH5 = [1] 122  16
#   write.csv(my.gsea.cat3, "../diff.exp.gene/DEG_SM/DEG_GSEA-3/my.gsea.cat3.csv", row.names = FALSE)
#   gsea.enrich.cat3.all=cbind(gsea.enrich.cat3.all,my.gsea.cat3)
# }
# dim(gsea.enrich.all)
# 
# my.gsea.cat3_M9 = read.csv("../diff.exp.gene/DEG_SM/DEG_GSEA-3/my.gsea.cat3_M9.csv", header = T)
# my.gsea.cat3_SH2 = read.csv("../diff.exp.gene/DEG_SM/DEG_GSEA-3/my.gsea.cat3_SH2.csv", header = T)
# my.gsea.cat3_SH5 = read.csv("../diff.exp.gene/DEG_SM/DEG_GSEA-3/my.gsea.cat3_SH5.csv", header = T)
# gsea.enrich.cat3.all = as.data.frame(cbind(my.gsea.cat3_M9, my.gsea.cat3_SH2, my.gsea.cat3_SH5))
# 
# #gsea.go.sig = gsea.go[gsea.go$`p (non-dir.)` < 0.05,]
# rownames(gsea.enrich.cat3.all) = gsea.enrich.cat3.all[,"M9_Name"]
# #write.csv(gsea.enrich.cat3.all, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/gsea.enrich.cat3.all.csv")
# 
# #gsea.enrich.pv =  gsea.enrich.all[,grep("LB_p(non-dir.)", colnames(gsea.enrich.all))]
# gsea.cat3.enrich.ND =  gsea.enrich.cat3.all[,grep("non-dir.", colnames(gsea.enrich.cat3.all))]
# colnames(gsea.cat3.enrich.ND) = c("M9_Stat(non-dir.)",
#                                   "M9_p(non-dir.)",
#                                   "M9_padj(non-dir.)",
#                                   "SH2_Stat(non-dir.)",
#                                   "SH2_p(non-dir.)",
#                                   "SH2_padj(non-dir.)",
#                                   "SH5_Stat(non-dir.)",
#                                   "SH5_p(non-dir.)",
#                                   "SH5_padj(non-dir.)")
# 
# colnames(gsea.cat3.enrich.ND) = c("M9_Stat_ND",
#                                   "M9_p_ND",
#                                   "M9_padj_ND",
#                                   "SH2_Stat_ND",
#                                   "SH2_p_ND",
#                                   "SH2_padj_ND",
#                                   "SH5_Stat_ND",
#                                   "SH5_p_ND",
#                                   "SH5_padj_ND")
# 
# gsea.cat3.enrich = matrix(as.numeric(unlist(gsea.cat3.enrich.ND)), nrow = nrow(gsea.cat3.enrich.ND))
# rownames(gsea.cat3.enrich) = rownames(gsea.cat3.enrich.ND)
# colnames(gsea.cat3.enrich) = colnames(gsea.cat3.enrich.ND)
# gsea.cat3.enrich.sig = gsea.cat3.enrich[(as.numeric(as.character(gsea.cat3.enrich[,"M9_padj_ND"]))<0.05|
#                                            as.numeric(as.character(gsea.cat3.enrich[,"SH2_padj_ND"]))<0.05|
#                                            as.numeric(as.character(gsea.cat3.enrich[,"SH5_padj_ND"]))<0.05),]
# #write.csv(gsea.cat3.enrich.sig, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/gsea.cat3.enrich.sig.csv")
# 
# gsea.cat3.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat3.enrich.sig[,"M9_padj_ND"])),
#                                 as.numeric(as.character(gsea.cat3.enrich.sig[,"SH2_padj_ND"])),
#                                 as.numeric(as.character(gsea.cat3.enrich.sig[,"SH5_padj_ND"])))
# rownames(gsea.cat3.enrich.sig.pv) = rownames(gsea.cat3.enrich.sig)
# colnames(gsea.cat3.enrich.sig.pv) = c("GSEA_Cat3_M9_padj","GSEA_Cat3_SH2_padj","GSEA_Cat3_SH5_padj")
# gsea.cat3.enrich.sig.pv.log10 = -log10(gsea.cat3.enrich.sig.pv)
# gsea.cat3.enrich.sig.pv.log10[gsea.cat3.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!
# 
# annotation.cat3.cat1 = data.frame(cat1 = unlist(lapply(rownames(gsea.cat3.enrich.sig.pv.log10), function(x)
#   geneCategories$category1[geneCategories$category3==x][1])))
# rownames(annotation.cat3.cat1) = rownames(gsea.cat3.enrich.sig.pv.log10)
# #write.csv(annotation.cat3.cat1, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/annotation.cat3.cat1.csv")
# 
# pheatmap(gsea.cat3.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
#          cluster_rows=TRUE, cluster_cols=FALSE,
#          cellwidth = 25, cellheight = 10,
#          show_rownames = T, show_colnames = T, 
#          fontsize=8, legend=TRUE,
#          headerplot = "GSEA_Cat3_Transcriptomics_Analysis",
#          main = "GSEA SubtiWiki Category 3 \n(padj) - Transcriptomics",
#          annotation_row = annotation.cat3.cat1
# )
# 
# 
# #######################################
# 
# #17/07/19
# # 
# # library(piano)  
# # 
# # 
# # DEG_M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_M9_minus_spoVG_upp.csv", header = T)
# # rownames(DEG_M9) = DEG_M9$X
# # DEG_M9 = DEG_M9[,-1]
# # 
# # DEG_SH2 = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_SH2_minus_spoVG_upp.csv", header = T)
# # rownames(DEG_SH2) = DEG_SH2$X
# # DEG_SH2 = DEG_SH2[,-1]
# # 
# # DEG_SH5 = read.csv("../diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_SH5_minus_spoVG_upp.csv", header = T)
# # rownames(DEG_SH5) = DEG_SH5$X
# # DEG_SH5 = DEG_SH5[,-1]
# # 
# # geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = T)
# # myGsc.cat3 = loadGSC(cbind(as.character(geneCategories$gene),
# #                            as.character(geneCategories$category3)))
# # # DEG_all = list(DEG_SH2, DEG_SH5)
# # # names(DEG_all) = c("SH2","SH5")
# # # 
# # # gsea.enrich.cat3.all = c()
# # # for (i in 1:length(DEG_all)){
# # #   my.de = DEG_all[[i]]
# # #   my.pval = my.de$padj
# # #   my.fc = my.de$log2FoldChange
# # #   names(my.pval) = rownames(my.de)
# # #   names(my.fc) = rownames(my.de)
# # #   my.pval = my.pval[!is.na(my.pval)] 
# # #   my.fc = my.pval[!is.na(my.pval)]
# # #   #m9 pval & fc = 6092, m9 pval - na = 5383, m9 fc - na = 6092
# # #   #SH2 pval & fc = 5776, SH2 pval - na = 5769, SH2 fc - na = 5776
# # #   #SH5 pval & fc = 6184, SH5 pval - na = 4864, SH5 fc - na = 6184
# # #   # my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc #m9 pval = NA BSU11750
# # #   my.fc = my.fc[!is.na(my.fc)]     
# # #   my.gsaRes.cat3 = runGSA(geneLevelStats=my.pval,directions=my.fc,
# # #                           gsc=myGsc.cat3,geneSetStat="reporter",adjMethod="BH")
# # #   my.gsea.cat3 = GSAsummaryTable(my.gsaRes.cat3)
# # #   colnames(my.gsea.cat3) = paste0(names(DEG_all)[i],"_",colnames(my.gsea.cat3))
# # #   my.gsea.cat3 = as.matrix(my.gsea.cat3)
# # #   gsea.enrich.cat3.all=cbind(gsea.enrich.cat3.all,my.gsea.cat3)
# # # }
# # 
# # DEG_all = list(DEG_M9, DEG_SH2, DEG_SH5)
# # names(DEG_all) = c("M9", "SH2","SH5")
# # 
# # gsea.enrich.cat3.all = c()
# # for (i in 1:length(DEG_all)){
# #   my.de = DEG_all[[i]]
# #   my.pval = my.de$padj
# #   my.fc = my.de$log2FoldChange
# #   names(my.pval) = rownames(my.de)
# #   names(my.fc) = rownames(my.de)
# #   my.pval = my.pval[!is.na(my.pval)] 
# #   my.fc = my.pval[!is.na(my.pval)]
# #   # my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc #m9 pval = NA BSU11750
# #   my.fc = my.fc[!is.na(my.fc)]     
# #   my.gsaRes.cat3 = runGSA(geneLevelStats=my.pval,directions=my.fc,
# #                           gsc=myGsc.cat3,geneSetStat="reporter",adjMethod="BH")
# #   my.gsea.cat3 = GSAsummaryTable(my.gsaRes.cat3)
# #   colnames(my.gsea.cat3) = paste0(names(DEG_all)[i],"_",colnames(my.gsea.cat3))
# #   my.gsea.cat3 = as.matrix(my.gsea.cat3)
# #   gsea.enrich.cat3.all=cbind(gsea.enrich.cat3.all,my.gsea.cat3)
# # }
# # 
# # #dim(gsea.enrich.all)
# # #[1] 111  76
# # 
# # #gsea.go.sig = gsea.go[gsea.go$`p (non-dir.)` < 0.05,]
# # rownames(gsea.enrich.cat3.all) = gsea.enrich.cat3.all[,"SH2_Name"]
# # 
# # #write.csv(gsea.enrich.cat3.all, "../diff.exp.gene/DEG_SM/DEG_GSEA-3/gsea.enrich.cat3.all_transcriptomics_M9.csv")
# # 
# # #gsea.enrich.pv =  gsea.enrich.all[,grep("LB_p(non-dir.)", colnames(gsea.enrich.all))]
# # gsea.cat3.enrich.ND =  gsea.enrich.cat3.all[,grep("(non-dir.)", colnames(gsea.enrich.cat3.all))]
# # colnames(gsea.cat3.enrich.ND) = c("M9_Stat(non-dir.)",
# #                                   "M9_p(non-dir.)",
# #                                   "M9_padj(non-dir.)",
# #                                   "SH2_Stat(non-dir.)",
# #                                   "SH2_p(non-dir.)",
# #                                   "SH2_padj(non-dir.)",
# #                                   "SH5_Stat(non-dir.)",
# #                                   "SH5_p(non-dir.)",
# #                                   "SH5_padj(non-dir.)")
# # 
# # 
# # colnames(gsea.cat3.enrich.ND) = c("LB_Stat_ND",
# #                                   "LB_p_ND",
# #                                   "LB_padj_ND",
# #                                   "M9_Stat_ND",
# #                                   "M9_p_ND",
# #                                   "M9_padj_ND",
# #                                   "SH2_Stat_ND",
# #                                   "SH2_p_ND",
# #                                   "SH2_padj_ND",
# #                                   "SH5_Stat_ND",
# #                                   "SH5_p_ND",
# #                                   "SH5_padj_ND")
# # 
# # gsea.cat3.enrich = matrix(as.numeric(unlist(gsea.cat3.enrich.ND)), nrow = nrow(gsea.cat3.enrich.ND))
# # rownames(gsea.cat3.enrich) = rownames(gsea.cat3.enrich.ND)
# # colnames(gsea.cat3.enrich) = colnames(gsea.cat3.enrich.ND)
# # 
# # gsea.cat3.enrich.sig = gsea.cat3.enrich[(as.numeric(as.character(gsea.cat3.enrich[,"LB_padj_ND"]))<0.05|
# #                                            as.numeric(as.character(gsea.cat3.enrich[,"M9_padj_ND"]))<0.05|
# #                                            as.numeric(as.character(gsea.cat3.enrich[,"SH2_padj_ND"]))<0.05|
# #                                            as.numeric(as.character(gsea.cat3.enrich[,"SH5_padj_ND"]))<0.05),]
# # 
# # 
# # #write.csv(gsea.cat3.enrich.sig, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/gsea.cat3.enrich.sig.csv")
# # 
# # 
# # gsea.cat3.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat3.enrich.sig[,"LB_padj_ND"])),
# #                                 as.numeric(as.character(gsea.cat3.enrich.sig[,"M9_padj_ND"])),
# #                                 as.numeric(as.character(gsea.cat3.enrich.sig[,"SH2_padj_ND"])),
# #                                 as.numeric(as.character(gsea.cat3.enrich.sig[,"SH5_padj_ND"])))
# # 
# # 
# # rownames(gsea.cat3.enrich.sig.pv) = rownames(gsea.cat3.enrich.sig)
# # colnames(gsea.cat3.enrich.sig.pv) = c("GSEA_Cat3_LB_padj","GSEA_Cat3_M9_padj","GSEA_Cat3_SH2_padj","GSEA_Cat3_SH5_padj")
# # gsea.cat3.enrich.sig.pv.log10 = -log10(gsea.cat3.enrich.sig.pv)
# # gsea.cat3.enrich.sig.pv.log10[gsea.cat3.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!
# # 
# # annotation.cat3.cat1 = data.frame(cat1 = unlist(lapply(rownames(gsea.cat3.enrich.sig.pv.log10), function(x)
# #   geneCategories$category1[geneCategories$category3==x][1])))
# # rownames(annotation.cat3.cat1) = rownames(gsea.cat3.enrich.sig.pv.log10)
# # 
# # #write.csv(annotation.cat3.cat1, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/annotation.cat3.cat1.csv")
# # 
# # pheatmap(gsea.cat3.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
# #          cluster_rows=TRUE, cluster_cols=FALSE,
# #          cellwidth = 25, cellheight = 10,
# #          show_rownames = T, show_colnames = T, 
# #          fontsize=8, legend=TRUE,
# #          headerplot = "GSEA_Cat3_Transcriptomics_Analysis",
# #          main = "GSEA SubtiWiki Category 3 \n(padj) - Transcriptomics",
# #          annotation_row = annotation.cat3.cat1
# # )
# # 
# # 
# # ######## cat1
# # 
# # 
# # library(piano)  
# # geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = T)
# # myGsc.cat1 = loadGSC(cbind(as.character(geneCategories$gene),
# #                            as.character(geneCategories$category1)))
# # DEG_all = list(DEG_M9, DEG_SH2, DEG_SH5)
# # names(DEG_all) = c("M9", "SH2","SH5")
# # 
# # gsea.enrich.cat1.all = c()
# # for (i in 1:length(DEG_all)){
# #   my.de = DEG_all[[i]]
# #   my.pval = my.de$padj
# #   my.fc = my.de$log2FoldChange
# #   names(my.pval) = rownames(my.de)
# #   names(my.fc) = rownames(my.de)
# #   my.pval = my.pval[!is.na(my.pval)] 
# #   my.fc = my.pval[!is.na(my.pval)]
# #   #m9 pval & fc = 6092, m9 pval - na = 5383, m9 fc - na = 6092
# #   #SH2 pval & fc = 5776, SH2 pval - na = 5769, SH2 fc - na = 5776
# #   #SH5 pval & fc = 6184, SH5 pval - na = 4864, SH5 fc - na = 6184
# #   # my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc #m9 pval = NA BSU11750
# #   my.fc = my.fc[!is.na(my.fc)]     
# #   my.gsaRes.cat1 = runGSA(geneLevelStats=my.pval,directions=my.fc,
# #                           gsc=myGsc.cat1,geneSetStat="reporter",adjMethod="BH")
# #   my.gsea.cat1 = GSAsummaryTable(my.gsaRes.cat1)
# #   colnames(my.gsea.cat1) = paste0(names(DEG_all)[i],"_",colnames(my.gsea.cat1))
# #   my.gsea.cat1 = as.matrix(my.gsea.cat1)
# #   gsea.enrich.cat1.all=cbind(gsea.enrich.cat1.all,my.gsea.cat1)
# # }
# # 
# # #dim(gsea.enrich.cat1.all)
# # #[1]  6 48
# # 
# # #gsea.go.sig = gsea.go[gsea.go$`p (non-dir.)` < 0.05,]
# # rownames(gsea.enrich.cat1.all) = gsea.enrich.cat1.all[,"M9_Name"]
# # 
# # #write.csv(gsea.enrich.cat1.all, "../diff.exp.gene/DEG_SM/DEG_GSEA-1/gsea.enrich.cat1.all_transcriptomics_M9.csv")
# # 
# # #gsea.enrich.pv =  gsea.enrich.all[,grep("LB_p(non-dir.)", colnames(gsea.enrich.all))]
# # gsea.cat1.enrich.ND =  gsea.enrich.cat1.all[,grep("(non-dir.)", colnames(gsea.enrich.cat1.all))]
# # colnames(gsea.cat1.enrich.ND) = c("M9_Stat(non-dir.)",
# #                                   "M9_p(non-dir.)",
# #                                   "M9_padj(non-dir.)",
# #                                   "SH2_Stat(non-dir.)",
# #                                   "SH2_p(non-dir.)",
# #                                   "SH2_padj(non-dir.)",
# #                                   "SH5_Stat(non-dir.)",
# #                                   "SH5_p(non-dir.)",
# #                                   "SH5_padj(non-dir.)")
# # 
# # 
# # colnames(gsea.cat1.enrich.ND) = c("M9_Stat_ND",
# #                                   "M9_p_ND",
# #                                   "M9_padj_ND",
# #                                   "SH2_Stat_ND",
# #                                   "SH2_p_ND",
# #                                   "SH2_padj_ND",
# #                                   "SH5_Stat_ND",
# #                                   "SH5_p_ND",
# #                                   "SH5_padj_ND")
# # 
# # gsea.cat1.enrich = matrix(as.numeric(unlist(gsea.cat1.enrich.ND)), nrow = nrow(gsea.cat1.enrich.ND))
# # rownames(gsea.cat1.enrich) = rownames(gsea.cat1.enrich.ND)
# # colnames(gsea.cat1.enrich) = colnames(gsea.cat1.enrich.ND)
# # 
# # gsea.cat1.enrich.sig = gsea.cat1.enrich[(as.numeric(as.character(gsea.cat1.enrich[,"M9_padj_ND"]))<0.05|
# #                                            as.numeric(as.character(gsea.cat1.enrich[,"SH2_padj_ND"]))<0.05|
# #                                            as.numeric(as.character(gsea.cat1.enrich[,"SH5_padj_ND"]))<0.05),]
# # 
# # 
# # #write.csv(gsea.cat1.enrich.sig, "../diff.exp.gene/DEG_SM/DEG_GSEA-1/gsea.cat1.enrich.sig.csv")
# # 
# # 
# # gsea.cat1.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat1.enrich.sig[,"M9_padj_ND"])),
# #                                 as.numeric(as.character(gsea.cat1.enrich.sig[,"SH2_padj_ND"])),
# #                                 as.numeric(as.character(gsea.cat1.enrich.sig[,"SH5_padj_ND"])))
# # 
# # 
# # rownames(gsea.cat1.enrich.sig.pv) = rownames(gsea.cat1.enrich.sig)
# # colnames(gsea.cat1.enrich.sig.pv) = c("GSEA_Cat1_M9_padj","GSEA_Cat1_SH2_padj","GSEA_Cat1_SH5_padj")
# # gsea.cat1.enrich.sig.pv.log10 = -log10(gsea.cat1.enrich.sig.pv)
# # gsea.cat1.enrich.sig.pv.log10[gsea.cat1.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!
# # 
# # pdf("../diff.exp.gene/DEG_SM/Figures/GSEA-1/gsea.cat1.enrich.sig.pv.log10_transcriptomics.pdf",width=10, height=10)
# # pheatmap(gsea.cat1.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
# #          cluster_rows=TRUE, cluster_cols=FALSE,
# #          cellwidth = 25, cellheight = 10,
# #          show_rownames = T, show_colnames = T, 
# #          fontsize=8, legend=TRUE,
# #          headerplot = "GSEA_Cat1_Transcriptomics_Analysis",
# #          main = "GSEA SubtiWiki Category 1 \n(padj) - Transcriptomics"
# # )
# # dev.off()
# 
# 
# 
# #######################
# 
# #21/07/19
# 
# #crude heatMap 
# 
# pheatmap(
#   counts.all_WT_SH5_vs_SH2,
#   scale="row",
#   cluster_rows=F, cluster_cols=T,
#   show_rownames = F, show_colnames = T, 
#   fontsize=10, legend=TRUE,
#   main = "HeatMap \nTotal Read Counts - SH5_vs_SH2\nTranscriptomic"
# )
# 
# 
# ####################################################################
# ### PCA Analysis  
# ####################################################################
# 
# 
# ######################################
# ### PCA - SH5_vs_SH2
# ######################################
# 
# pca = prcomp(t(counts.all_WT_SH5_vs_SH2), center = TRUE, scale = FALSE)
# percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
# pca2 = cbind(as.data.frame(pca$x), sampleTable_WT_SH5_vs_SH2)   #check the order of ko and WT
# 
# #write.csv(pca2, "./Combined /Filtered /CSV/20181217_PCA_Combined.csv")
# 
# #pdf("../Figures/LB/PCA_Plot_LB.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
# ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype,shape=pca2$media)) + geom_point(size=3) + 
#   xlab(paste0("PC1: ",percent.var[1],"% variance")) +
#   ylab(paste0("PC2: ",percent.var[2], "% variance")) +
#   labs(colour = "Genotype",shape="Media") +
#   ggtitle("PCA Plot - SH5_vs_SH2 \nTranscriptomic Analysis - Raw Data") + #CHANGE THE NAME HERE 
#   theme(
#     plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
#     axis.title.x = element_text(color="black", size=12, face="bold"),
#     axis.title.y = element_text(color="black", size=12, face="bold")
#   ) + 
#   #geom_text(aes(PC1, PC2, colour=conditions),label=pca2$genotype) +    #adds labels to each data point
#   coord_fixed()
# #dev.off()                  #turn this OFF if just want to see the picture in the Plots
# 
# 
