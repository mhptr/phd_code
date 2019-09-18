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

####################################################################
### 2. Set Working Directory 
####################################################################

setwd("~/OneDrive - University of Warwick/WORK/RESULTS/PROTEOMICS/FINAL Result")


a####################################################################
### 3. Set Parameters  
####################################################################

thr.fc=2 #thresold_FoldChange
thr.pv =0.05 #thresold_pValue1
swati.color = c(brewer.pal(8,"Dark2")) #TA
swati.color2 = c(brewer.pal(9, "Greens")) #SM
swati.color.Blues = c(brewer.pal(8,"Blues"))


####################################################################
### 4. Load Data/ .txt file 
####################################################################


data.pr.all = read.csv("./Perseus Output Files /combined/20181217_combined_BSU_clean_20190107.csv", header = T)
dim(data.pr.all)

sampleName=colnames(data.pr.all)[1:24]

sample.details = data.frame(cbind(sampleName=sampleName,
                                  media=unlist(lapply(sampleName,function(x)
                                    unlist(strsplit(as.character(x), "_"))[1]
                                  )),
                                  genotype=gsub("[1|2|3]", "", unlist(lapply(sampleName,function(x)  #global subtitute
                                    unlist(strsplit(as.character(x), "_"))[2]
                                  )))))

#experimentalist=unlist(lapply(sampleName,function(x)
#  unlist(strsplit(as.character(x), "_"))[1]
#))
#experimentalist[experimentalist %in% experimentalist[-grep("SH", experimentalist)]] = "Emma"
#experimentalist[experimentalist %in% experimentalist[grep("SH", experimentalist)]] = "Swati"

#sample.details$experimentalist = experimentalist

####################################################################
### 5. Annotate Data  
####################################################################


data.pr = data.pr.all[,1:24] #[Row,Column]
annotation.all = data.pr.all[,25:dim(data.pr.all)[2]] #[2] means the columns

#par(mfrow=c(4,4))
#for (i in 1:4) {
#  for(j in 1:4) {
#    plot(data.pr[,i], data.pr[,j])
#  }
#}


#Correlation Plots 
plot(data.pr[,1], data.pr[,2])

dim(data.pr) 

data.pr.v2 = data.pr[!apply(data.pr, 1, function(x) any(is.na(x))), ]
#pheatmap(data.pr.v2, scale = "row")




####################################################################
### 6. Set Conditions  
####################################################################

## We must make sure that the order in the condition is same as
## the column in the data.pr.all


#conditions = factor(c("SH2_D", "SH2_D", "SH2_D", "SH2_W", "SH2_W", "SH2_W", "SH5_D", "SH5_D", "SH5_D", "SH5_W", "SH5_W", "SH5_W", 
#                      "LB_D", "LB_D", "LB_D", "LB_W", "LB_W", "LB_W", "M9_D", "M9_D", "M9_D", "M9_W", "M9_W", "M9_W"))  





####################################################################
### 7. Imputate the Missing Values (LFQ Intensities) 
####################################################################


# minimun expression value across protein is replacement for NA
replaced.exp = min(data.pr[!is.na(data.pr)])  #derives the minimum value from the matrix data.pr
data.pr.imputed = c() # this is imputed expression value
for(i in 1:dim(data.pr)[1]) {
  my.pr.wt.LB = data.pr[i,paste0(sample.details$media,"_",sample.details$genotype)  == "LB_W"]
  my.pr.ko.LB = data.pr[i,paste0(sample.details$media,"_",sample.details$genotype) == "LB_D"]
  my.pr.wt.M9 = data.pr[i,paste0(sample.details$media,"_",sample.details$genotype) == "M9_W"]
  my.pr.ko.M9 = data.pr[i,paste0(sample.details$media,"_",sample.details$genotype) == "M9_D"]
  my.pr.wt.SH2 = data.pr[i,paste0(sample.details$media,"_",sample.details$genotype) == "SH2_W"]
  my.pr.ko.SH2 = data.pr[i,paste0(sample.details$media,"_",sample.details$genotype) == "SH2_D"]
  my.pr.wt.SH5 = data.pr[i,paste0(sample.details$media,"_",sample.details$genotype) == "SH5_W"]
  my.pr.ko.SH5 = data.pr[i,paste0(sample.details$media,"_",sample.details$genotype) == "SH5_D"]
  # imputing WT out of 3 replicates, if more than 2 NA then replae by min 
  # else replace by the avg of the two expression values
  if(sum(is.na(my.pr.wt.LB))>=2) {
    my.pr.wt.LB[is.na(my.pr.wt.LB)] = replaced.exp
  } else {
    my.pr.wt.LB[is.na(my.pr.wt.LB)] = mean(my.pr.wt.LB[!is.na(my.pr.wt.LB)])
  }
  if(sum(is.na(my.pr.wt.M9))>=2) {
    my.pr.wt.M9[is.na(my.pr.wt.M9)] = replaced.exp
  } else {
    my.pr.wt.M9[is.na(my.pr.wt.M9)] = mean(my.pr.wt.M9[!is.na(my.pr.wt.M9)])
  }
  if(sum(is.na(my.pr.wt.SH2))>=2) {
    my.pr.wt.SH2[is.na(my.pr.wt.SH2)] = replaced.exp
  } else {
    my.pr.wt.SH2[is.na(my.pr.wt.SH2)] = mean(my.pr.wt.SH2[!is.na(my.pr.wt.SH2)])
  }
  if(sum(is.na(my.pr.wt.SH5))>=2) {
    my.pr.wt.SH5[is.na(my.pr.wt.SH5)] = replaced.exp
  } else {
    my.pr.wt.SH5[is.na(my.pr.wt.SH5)] = mean(my.pr.wt.SH5[!is.na(my.pr.wt.SH5)])
  }
  # imputing Delta_spoVG  # imputing Delta_spoVG  # imputing Delta_spoVG  # imputing Delta_spoVG
  if(sum(is.na(my.pr.ko.LB))>=2) {
    my.pr.ko.LB[is.na(my.pr.ko.LB)] = replaced.exp
  } else {
    my.pr.ko.LB[is.na(my.pr.ko.LB)] = mean(my.pr.ko.LB[!is.na(my.pr.ko.LB)])
  }  
  if(sum(is.na(my.pr.ko.M9))>=2) {
    my.pr.ko.M9[is.na(my.pr.ko.M9)] = replaced.exp
  } else {
    my.pr.ko.M9[is.na(my.pr.ko.M9)] = mean(my.pr.ko.M9[!is.na(my.pr.ko.M9)])
  }     
  if(sum(is.na(my.pr.ko.SH2))>=2) {
    my.pr.ko.SH2[is.na(my.pr.ko.SH2)] = replaced.exp
  } else {
    my.pr.ko.SH2[is.na(my.pr.ko.SH2)] = mean(my.pr.ko.SH2[!is.na(my.pr.ko.SH2)])
  }     
  if(sum(is.na(my.pr.ko.SH5))>=2) {
    my.pr.ko.SH5[is.na(my.pr.ko.SH5)] = replaced.exp
  } else {
    my.pr.ko.SH5[is.na(my.pr.ko.SH5)] = mean(my.pr.ko.SH5[!is.na(my.pr.ko.SH5)])
  }     
  
  data.pr.imputed = rbind(data.pr.imputed,
                          cbind(my.pr.ko.LB,my.pr.wt.LB, my.pr.ko.M9,my.pr.wt.M9,my.pr.ko.SH2,my.pr.wt.SH2,my.pr.ko.SH5,my.pr.wt.SH5))    
  #check and alter the order of knockout(ko) or the wild-type(wt) 
  #based on the data.pr.all
}
dim(data.pr.imputed)

##pheatmap(data.pr.imputed, scale = "row")
#pheatmap(data.pr.imputed, scale="row",
#         cluster_rows=F, cluster_cols=T,
#         show_rownames = F, show_colnames = T, 
#         fontsize=10, legend=TRUE,
#         main = "HeatMap \nTotal Identified Proteins (Imputated) \nProteomics"
#)


sampleName = colnames(data.pr.imputed)
sample.details.imputed = data.frame(cbind(sampleName=sampleName,
                                          media=unlist(lapply(sampleName,function(x)
                                            unlist(strsplit(as.character(x), "_"))[1]
                                          )),
                                          genotype=gsub("[1|2|3]", "", unlist(lapply(sampleName,function(x)  #global subtitute
                                            unlist(strsplit(as.character(x), "_"))[2]
                                          )))))

experimentalist=unlist(lapply(sampleName,function(x)
  unlist(strsplit(as.character(x), "_"))[1]
))

sample.details.imputed$media_genot = paste0(sample.details.imputed$media,'.',
                                            sample.details.imputed$genotype)
#experimentalist[experimentalist %in% experimentalist[-grep("SH", experimentalist)]] = "Emma"
#experimentalist[experimentalist %in% experimentalist[grep("SH", experimentalist)]] = "Swati"
#sample.details.imputed$experimentalist = experimentalist



####################################################################
### 8. Assign Row Names 
####################################################################

rownames(data.pr.imputed) = annotation.all$T..ENSG   #Proteins with 1 BSU number only




####################################################################
### 9. PCA Analysis  
####################################################################


######################################
### 9.1 PCA - All
######################################

pca = prcomp(t(data.pr.imputed), center = TRUE, scale = FALSE)
percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
pca2 = cbind(as.data.frame(pca$x), sample.details.imputed)   #check the order of ko and WT

#write.csv(pca2, "./Combined /Filtered /CSV/20181217_PCA_Combined.csv")

#pdf("../Figures/LB/PCA_Plot_LB.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
#ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype,shape=pca2$media)) + geom_point(size=3) + 
#  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
#  ylab(paste0("PC2: ",percent.var[2], "% variance")) +
#  labs(colour = "Genotype",shape="Media") +
#  ggtitle("PCA Plot - All Conditions \nProteomic Analysis") + #CHANGE THE NAME HERE 
#  theme(
#    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
#    axis.title.x = element_text(color="black", size=12, face="bold"),
#    axis.title.y = element_text(color="black", size=12, face="bold")
#  ) + 
#  #geom_text(aes(PC1, PC2, colour=conditions),label=pca2$genotype) +    #adds labels to each data point
#  coord_fixed()
#dev.off()                  #turn this OFF if just want to see the picture in the Plots




################################################################
#9.2 PCA - Conditionwise
################################################################


data.pr.imputed.LB = data.pr.imputed[,1:6]
data.pr.imputed.M9 = data.pr.imputed[,7:12]
data.pr.imputed.SH2 = data.pr.imputed[,13:18]
data.pr.imputed.SH5 = data.pr.imputed[,19:24]

########

#pca = prcomp(t(data.pr.imputed.LB), center = TRUE, scale = FALSE)
#pca = prcomp(t(data.pr.imputed.M9), center = TRUE, scale = FALSE)
#pca = prcomp(t(data.pr.imputed.SH2), center = TRUE, scale = FALSE)
#pca = prcomp(t(data.pr.imputed.SH5), center = TRUE, scale = FALSE)
#percent.var = round(100*pca$sdev^2/sum(pca$sdev^2))
#pca2 = cbind(as.data.frame(pca$x), sample.details.imputed)   #check the order of ko and WT

#pdf("../Figures/LB/PCA_Plot_LB.pdf",width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
#ggplot(pca2, aes(PC1, PC2, colour=pca2$genotype)) + geom_point(size=3) + 
#  xlab(paste0("PC1: ",percent.var[1],"% variance")) +
#  ylab(paste0("PC2: ",percent.var[2], "% variance")) +
#  labs(colour = "Genotype",shape="Media") +
  #ggtitle("PCA Plot - LB \nProteomic Analysis") + #CHANGE THE NAME HERE 
  #ggtitle("PCA Plot - M9 \nProteomic Analysis") + 
  #ggtitle("PCA Plot - SH2 \nProteomic Analysis") +  
  #ggtitle("PCA Plot - SH5 \nProteomic Analysis") + 
#  theme(
#    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
#    axis.title.x = element_text(color="black", size=12, face="bold"),
#    axis.title.y = element_text(color="black", size=12, face="bold")
#  ) + 
#  #geom_text(aes(PC1, PC2, colour=conditions),label=pca2$genotype) +    #adds labels to each data point
#  coord_fixed()
#dev.off()                  #turn this OFF if just want to see the picture in the Plots






####################################################################
### 10. Set up the Data Frame   
####################################################################


#### 10.1. Design
design <- model.matrix(~0+sample.details.imputed$media_genot)


#colnames(design) <- levels(conditions)
colnames(design) = gsub("sample.details.imputed\\$media_genot","",colnames(design))  #"\\$" maked the dollar sign as non-special character
rownames(design) <- colnames(data.pr.imputed)


#### 10.2. lmfit

#genotype = unique(conditions)

####################################################################
### 5. Set REFERENCE and Making CONTRAST 
####################################################################

####################################################################
### 5.1 Del vs WT 
####################################################################


mc = makeContrasts(
  de_KO_vs_WT_in_LB = LB.D-LB.W,
  de_KO_vs_WT_in_M9 = M9.D-M9.W,
  de_KO_vs_WT_in_SH2 = SH2.D-SH2.W,
  de_KO_vs_WT_in_SH5 = SH5.D-SH5.W,
  #diff = (((LB.D-SH2.D)-M9.D)-SH5.D) - (((LB.W-SH2.W)-M9.W)-SH5.W),
  levels = design
)

### MAKING CONTRAST
lm.fit <- lmFit(data.pr.imputed, design)
c.fit = contrasts.fit(lm.fit, mc)
eb = eBayes(c.fit)


####################################################################
### 5. DIFFRENTIAL EXPRESSION (Delta vs WT)
####################################################################
#### 5. differentiually expressed proteins Knockout referene to WT
diffexp.pr.all.LB = topTable(eb, adjust="BH",coef = 1,number = dim(data.pr.imputed)[1])
diffexp.pr.all.M9 = topTable(eb, adjust="BH",coef = 2,number = dim(data.pr.imputed)[1])
diffexp.pr.all.SH2 = topTable(eb, adjust="BH",coef = 3,number = dim(data.pr.imputed)[1])
diffexp.pr.all.SH5 = topTable(eb, adjust="BH",coef = 4,number = dim(data.pr.imputed)[1])

#write.csv(diffexp.pr.all.LB,"./Analysis - R-Script /Data/combined/output/diff.exp./20181217_diffexp.pr.all_LB.csv")
#write.csv(diffexp.pr.all.M9,"./Analysis - R-Script /Data/combined/output/diff.exp./20181217_diffexp.pr.all_M9.csv")
#write.csv(diffexp.pr.all.SH2,"./Analysis - R-Script /Data/combined/output/diff.exp./20181217_diffexp.pr.all_SH2.csv")
#write.csv(diffexp.pr.all.SH5,"./Analysis - R-Script /Data/combined/output/diff.exp./20181217_diffexp.pr.all_SH5.csv")

#diffexp.pr.all.D_vs_W = topTable(eb, adjust="BH",coef = 5,number = dim(data.pr.imputed)[1])
#sum(diffexp.pr.all.D_vs_W$adj.P.Val<0.05)


##################
##  SIGNIFICANT
##################

diffexp.pr.all.LB_SIG = diffexp.pr.all.LB[(diffexp.pr.all.LB$adj.P.Val <0.05),]
diffexp.pr.all.M9_SIG = diffexp.pr.all.M9[(diffexp.pr.all.M9$adj.P.Val <0.05),]
diffexp.pr.all.SH2_SIG = diffexp.pr.all.SH2[(diffexp.pr.all.SH2$adj.P.Val <0.05),]
diffexp.pr.all.SH5_SIG = diffexp.pr.all.SH5[(diffexp.pr.all.SH5$adj.P.Val <0.05),]

#write.csv(diffexp.pr.all.LB_SIG,"./Analysis - R-Script /Data/combined/output/diff.exp.SIG/20181217_diffexp.pr.all_LB_SIG.csv")
#write.csv(diffexp.pr.all.M9_SIG,"./Analysis - R-Script /Data/combined/output/diff.exp.SIG/20181217_diffexp.pr.all_M9_SIG.csv")
#write.csv(diffexp.pr.all.SH2_SIG,"./Analysis - R-Script /Data/combined/output/diff.exp.SIG/20181217_diffexp.pr.all_SH2_SIG.csv")
#write.csv(diffexp.pr.all.SH5_SIG,"./Analysis - R-Script /Data/combined/output/diff.exp.SIG/20181217_diffexp.pr.all_SH5_SIG.csv")



##########################################################
##                     VENN
##########################################################

diffexp.pr.all.LB = read.csv("./Analysis - R-Script /Data/combined/output/diff.exp./diff.exp._minus_spoVG_upp/20190514_diffexp.pr.all_LB_minus_spoVG_upp.csv", header = T)
rownames(diffexp.pr.all.LB) = diffexp.pr.all.LB$X
diffexp.pr.all.LB = diffexp.pr.all.LB[,-1]

diffexp.pr.all.M9 = read.csv("./Analysis - R-Script /Data/combined/output/diff.exp./diff.exp._minus_spoVG_upp/20190514_diffexp.pr.all_M9_minus_spoVG_upp.csv", header = T)
rownames(diffexp.pr.all.M9) = diffexp.pr.all.M9$X
diffexp.pr.all.M9 = diffexp.pr.all.M9[,-1]

diffexp.pr.all.SH2 = read.csv("./Analysis - R-Script /Data/combined/output/diff.exp./diff.exp._minus_spoVG_upp/20190514_diffexp.pr.all_SH2_minus_spoVG_upp.csv", header = T)
rownames(diffexp.pr.all.SH2) = diffexp.pr.all.SH2$X
diffexp.pr.all.SH2 = diffexp.pr.all.SH2[,-1]

diffexp.pr.all.SH5 = read.csv("./Analysis - R-Script /Data/combined/output/diff.exp./diff.exp._minus_spoVG_upp/20190514_diffexp.pr.all_SH5_minus_spoVG_upp.csv", header = T)
rownames(diffexp.pr.all.SH5) = diffexp.pr.all.SH5$X
diffexp.pr.all.SH5 = diffexp.pr.all.SH5[,-1]


rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05]
rownames(diffexp.pr.all.M9)[diffexp.pr.all.M9$adj.P.Val<0.05]
rownames(diffexp.pr.all.SH2)[diffexp.pr.all.SH2$adj.P.Val<0.05]
rownames(diffexp.pr.all.SH5)[diffexp.pr.all.SH5$adj.P.Val<0.05]
#rownames(diffexp.pr.all.D_vs_W)[diffexp.pr.all.D_vs_W$adj.P.Val<0.05]


venn.plot = venn.diagram(
  list("LB" = rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05],
       "M9" = rownames(diffexp.pr.all.M9)[diffexp.pr.all.M9$adj.P.Val<0.05],
       "SH2" = rownames(diffexp.pr.all.SH2)[diffexp.pr.all.SH2$adj.P.Val<0.05],
       "SH5" = rownames(diffexp.pr.all.SH5)[diffexp.pr.all.SH5$adj.P.Val<0.05]),
  filename = NULL,
  col="transparent",
  fill=swati.color[1:4], height = 3000, width = 3000,resolution =500, imagetype = "tiff",
  main = "Venn Diagram - All Conditions \nProteomic Analysis",main.pos  = c(0.5, 1.05), main.fontface = "bold",
  main.fontfamily = "Helvetica", main.col = "Dark Blue", 
  main.cex = 1.5, cat.cex = 2, cat.fontfamily = "Helvetica", 
  cex = 2.3, fontfamily = "Helvetica",fontface = "bold"
)
grid.newpage()
grid.draw(venn.plot)


#################################################################
###                   Heat_Maps_Diff.SIG_exp_proteins
##################################################################

pheatmap(data.pr.imputed, scale="row",
         cluster_rows=F, cluster_cols=T,
         show_rownames = F, show_colnames = T, 
         fontsize=10, legend=TRUE,
         main = "HeatMap \nTotal Identified Proteins (Imputated) \nProteomics")

#data.pr.imputed_LB = data.pr.imputed[,1:6]
#write.csv(data.pr.imputed_LB, "./Analysis - R-Script /Data/combined/output/imputated/data.pr.imputed_LB.csv")

data.pr.imputed_LB = read.csv("./Analysis - R-Script /Data/combined/output/imputated/data.pr.imputed_LB.csv", header = T)
diffexp.pr.all.LB_SIG_annot = read.csv("./Analysis - R-Script /Data/combined/output/diffexp.pr.all.LB_SIG_Annotated.csv", header = T)
data.pr.imputed_LB_diff.exp.SIG = merge(diffexp.pr.all.LB_SIG, data.pr.imputed_LB, by.x="X",by.y="X",all.x=TRUE)
rownames(data.pr.imputed_LB_diff.exp.SIG) = diffexp.pr.all.LB_SIG_annot$name
data.pr.imputed_LB_diff.exp.SIG = data.pr.imputed_LB_diff.exp.SIG[,-1]
#write.csv(data.pr.imputed_LB_diff.exp.SIG, "./Analysis - R-Script /Data/combined/output/data.pr.imputed_LB_diff.exp.SIG.csv")
data.pr.imputed_LB_diff.exp.SIG_top10 = read.csv("./Analysis - R-Script /Data/combined/output/imputated_diff.exp_SIG/data.pr.imputed_LB_diff.exp.SIG_top10.csv", header = T)
rownames(data.pr.imputed_LB_diff.exp.SIG_top10) = data.pr.imputed_LB_diff.exp.SIG_top10$X
#data.pr.imputed_LB_diff.exp.SIG_Down = data.pr.imputed_LB_diff.exp.SIG_Down[,-1]

data.pr.imputed_LB_diff.exp.SIG_top10_sample = data.pr.imputed_LB_diff.exp.SIG_top10[,8:13]
rownames(data.pr.imputed_LB_diff.exp.SIG_top10_sample) = rownames(data.pr.imputed_LB_diff.exp.SIG_top10_sample)
colnames(data.pr.imputed_LB_diff.exp.SIG_top10_sample) = colnames(data.pr.imputed_LB_diff.exp.SIG_top10_sample)


pheatmap(data.pr.imputed_LB_diff.exp.SIG_top10_sample,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 25, cellheight = 10,
         show_rownames = T, show_colnames = T, 
         fontsize=8, legend=TRUE,
         headerplot = "LB_diff_exp_SIG_Down",
         main = "LB_diff_exp_SIG_Down"
         #annotation_row = annotation.cat3.cat1
)





gsea.cat3.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat3.enrich.sig[,"LB_padj_ND"])),
                                as.numeric(as.character(gsea.cat3.enrich.sig[,"M9_padj_ND"])),
                                as.numeric(as.character(gsea.cat3.enrich.sig[,"SH2_padj_ND"])),
                                as.numeric(as.character(gsea.cat3.enrich.sig[,"SH5_padj_ND"])))


rownames(gsea.cat3.enrich.sig.pv) = rownames(gsea.cat3.enrich.sig)
colnames(gsea.cat3.enrich.sig.pv) = c("GSEA_Cat3_LB_adj.P.Val","GSEA_Cat3_M9_adj.P.Val","GSEA_Cat3_SH2_adj.P.Val","GSEA_Cat3_SH5_adj.P.Val")
gsea.cat3.enrich.sig.pv.log10 = -log10(gsea.cat3.enrich.sig.pv)
gsea.cat3.enrich.sig.pv.log10[gsea.cat3.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!

annotation.cat3.cat1 = data.frame(cat1 = unlist(lapply(rownames(gsea.cat3.enrich.sig.pv.log10), function(x)
  geneCategories$category1[geneCategories$category3==x][1])))
rownames(annotation.cat3.cat1) = rownames(gsea.cat3.enrich.sig.pv.log10)

#write.csv(annotation.cat3.cat1, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/annotation.cat3.cat1.csv")

pheatmap(gsea.cat3.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 25, cellheight = 10,
         show_rownames = T, show_colnames = T, 
         fontsize=8, legend=TRUE,
         headerplot = "GSEA_Cat3_Proteomics_Analysis",
         main = "GSEA SubtiWiki Category 3 \n(adj.P.Val) - Proteomics",
         annotation_row = annotation.cat3.cat1
)




#################################################################
###                   Volcano Plots
##################################################################


thr.adj.pv =0.05

g = ggplot(#data=diffexp.pr.all.LB, 
           #data=diffexp.pr.all.M9, 
           #data=diffexp.pr.all.SH2, 
           #data=diffexp.pr.all.SH5,
           aes(title, x=logFC, y=-log10(adj.P.Val), colour=adj.P.Val<0.05)) +   #colour=-log10(P.Value)

  #ggtitle("Volcano Plot - LB \nProteomic Analysis") +      #CHANGE THE NAME HERE 
  #ggtitle("Volcano Plot - M9 \nProteomic Analysis") +
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


bsu.kegg.pathway = read.csv("./Analysis - R-Script /Data/bsu.kegg.pathway.csv")
diffexp.pr.all.conds = list(diffexp.pr.all.LB,diffexp.pr.all.M9,diffexp.pr.all.SH2,diffexp.pr.all.SH5)
names(diffexp.pr.all.conds) = c("LB","M9","SH2","SH5")
kegg.enrich.all = c()                              #make dataframe 
for(i in 1:length(diffexp.pr.all.conds)){
  my.de = diffexp.pr.all.conds[[i]]
  #my.de.genes = rownames(my.de)[my.de$P.Value<0.05] 
  my.de.genes = rownames(my.de)[my.de$adj.P.Val<0.05]
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
  colnames(my.kegg.enrich) = paste0(names(diffexp.pr.all.conds)[i],"_",colnames(my.kegg.enrich))
  kegg.enrich.all = as.data.frame(cbind(kegg.enrich.all,my.kegg.enrich)) 
}
rownames(kegg.enrich.all) = unique(bsu.kegg.pathway$keggPathName)
#kegg.enrich.all = kegg.enrich.all[!apply(kegg.enrich.all, 1, function(x) any(is.na(x))),]

kegg.enrich.all = kegg.enrich.all[as.numeric(as.character(kegg.enrich.all$LB_k))>2,]
kegg.enrich.sig = kegg.enrich.all[(as.numeric(as.character(kegg.enrich.all$LB_pv))<0.05|
                                     as.numeric(as.character(kegg.enrich.all$M9_pv))<0.05|
                                     as.numeric(as.character(kegg.enrich.all$SH2_pv))<0.05|
                                     as.numeric(as.character(kegg.enrich.all$SH5_pv))<0.05),]
#write.csv(kegg.enrich.sig, "./Analysis - R-Script /Data/combined/output/KEGG/kegg.enrich.sig.csv")


kegg.enrich.sig.pv = cbind(as.numeric(as.character(kegg.enrich.sig$LB_pv)),
                           as.numeric(as.character(kegg.enrich.sig$M9_pv)),
                           as.numeric(as.character(kegg.enrich.sig$SH2_pv)),
                           as.numeric(as.character(kegg.enrich.sig$SH5_pv)))


rownames(kegg.enrich.sig.pv) = rownames(kegg.enrich.sig)
colnames(kegg.enrich.sig.pv) = c("KEGG_LB_adj.P.Val","KEGG_M9_adj.P.Val","KEGG_SH2_adj.P.Val","KEGG_SH5_adj.P.Val")

#write.csv(kegg.enrich.sig.pv, "./Analysis - R-Script /Data/combined/output/KEGG/kegg.enrich.sig.pv.csv")

kegg.enrich.sig.pv.log10 = -log10(kegg.enrich.sig.pv)
kegg.enrich.sig.pv.log10[kegg.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!

pheatmap(kegg.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 70, cellheight = 27,
         show_rownames = T, show_colnames = T, 
         fontsize=12, legend=TRUE,
         headerplot = "KEEG_Proteomics_Analysis",
         main = "KEGG Pathway Enrichment Analysis \n(adj.P.Val) - Proteomics"
         )



##################################################################################################
###### Gene Set Enrichment Analysis - Category 3 (GSEA-3) with adj.P.Val
##################################################################################################


#source("http://www.bioconductor.org/biocLite.R")
#BiocManager::install("piano", dependencies=TRUE)          #this worked
geneCategories = read.csv("./Analysis - R-Script /SubtiWiki Exports /geneCategories.csv")
myGsc.cat3 = loadGSC(cbind(as.character(geneCategories$gene),
                          as.character(geneCategories$category3)))
diffexp.pr.all.conds = list(diffexp.pr.all.LB,diffexp.pr.all.M9,diffexp.pr.all.SH2,diffexp.pr.all.SH5)
names(diffexp.pr.all.conds) = c("LB","M9","SH2","SH5")

gsea.enrich.cat3.all = c()
for (i in 1:length(diffexp.pr.all.conds)){
  my.de = diffexp.pr.all.conds[[i]]
  my.pval = my.de$adj.P.Val
  my.fc = my.de$logFC
  names(my.pval) = rownames(my.de)
  names(my.fc) = rownames(my.de)
  my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc
  my.fc = my.fc[!is.na(my.fc)]    #1707 entries 
  my.gsaRes.cat3 = runGSA(geneLevelStats=my.pval,directions=my.fc,
                        gsc=myGsc.cat3,geneSetStat="reporter",adjMethod="BH")
  my.gsea.cat3 = GSAsummaryTable(my.gsaRes.cat3)
  colnames(my.gsea.cat3) = paste0(names(diffexp.pr.all.conds)[i],"_",colnames(my.gsea.cat3))
  my.gsea.cat3 = as.matrix(my.gsea.cat3)
  gsea.enrich.cat3.all=cbind(gsea.enrich.cat3.all,my.gsea.cat3)
}

#dim(gsea.enrich.all)
#[1] 111  76

#gsea.go.sig = gsea.go[gsea.go$`p (non-dir.)` < 0.05,]
rownames(gsea.enrich.cat3.all) = gsea.enrich.cat3.all[,"LB_Name"]

#write.csv(gsea.enrich.cat3.all, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/gsea.enrich.cat3.all.csv")

#gsea.enrich.pv =  gsea.enrich.all[,grep("LB_p(non-dir.)", colnames(gsea.enrich.all))]
gsea.cat3.enrich.ND =  gsea.enrich.cat3.all[,grep("(non-dir.)", colnames(gsea.enrich.cat3.all))]
colnames(gsea.cat3.enrich.ND) = c("LB_Stat(non-dir.)",
                                  "LB_p(non-dir.)",
                                  "LB_padj(non-dir.)",
                                  "M9_Stat(non-dir.)",
                                  "M9_p(non-dir.)",
                                  "M9_padj(non-dir.)",
                                  "SH2_Stat(non-dir.)",
                                  "SH2_p(non-dir.)",
                                  "SH2_padj(non-dir.)",
                                  "SH5_Stat(non-dir.)",
                                  "SH5_p(non-dir.)",
                                  "SH5_padj(non-dir.)")


colnames(gsea.cat3.enrich.ND) = c("LB_Stat_ND",
                                  "LB_p_ND",
                                  "LB_padj_ND",
                                  "M9_Stat_ND",
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

gsea.cat3.enrich.sig = gsea.cat3.enrich[(as.numeric(as.character(gsea.cat3.enrich[,"LB_padj_ND"]))<0.05|
                                     as.numeric(as.character(gsea.cat3.enrich[,"M9_padj_ND"]))<0.05|
                                     as.numeric(as.character(gsea.cat3.enrich[,"SH2_padj_ND"]))<0.05|
                                     as.numeric(as.character(gsea.cat3.enrich[,"SH5_padj_ND"]))<0.05),]


#write.csv(gsea.cat3.enrich.sig, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/gsea.cat3.enrich.sig.csv")


gsea.cat3.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat3.enrich.sig[,"LB_padj_ND"])),
                           as.numeric(as.character(gsea.cat3.enrich.sig[,"M9_padj_ND"])),
                           as.numeric(as.character(gsea.cat3.enrich.sig[,"SH2_padj_ND"])),
                           as.numeric(as.character(gsea.cat3.enrich.sig[,"SH5_padj_ND"])))


rownames(gsea.cat3.enrich.sig.pv) = rownames(gsea.cat3.enrich.sig)
colnames(gsea.cat3.enrich.sig.pv) = c("GSEA_Cat3_LB_adj.P.Val","GSEA_Cat3_M9_adj.P.Val","GSEA_Cat3_SH2_adj.P.Val","GSEA_Cat3_SH5_adj.P.Val")
gsea.cat3.enrich.sig.pv.log10 = -log10(gsea.cat3.enrich.sig.pv)
gsea.cat3.enrich.sig.pv.log10[gsea.cat3.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!

annotation.cat3.cat1 = data.frame(cat1 = unlist(lapply(rownames(gsea.cat3.enrich.sig.pv.log10), function(x)
  geneCategories$category1[geneCategories$category3==x][1])))
rownames(annotation.cat3.cat1) = rownames(gsea.cat3.enrich.sig.pv.log10)

#write.csv(annotation.cat3.cat1, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/annotation.cat3.cat1.csv")

pheatmap(gsea.cat3.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 25, cellheight = 10,
         show_rownames = T, show_colnames = T, 
         fontsize=8, legend=TRUE,
         headerplot = "GSEA_Cat3_Proteomics_Analysis",
         main = "GSEA SubtiWiki Category 3 \n(adj.P.Val) - Proteomics",
         annotation_row = annotation.cat3.cat1
)




##################################################################################################
###### Gene Set Enrichment Analysis - Category 1 (GSEA-1) - adj.P.Val
##################################################################################################


#source("http://www.bioconductor.org/biocLite.R")
#BiocManager::install("piano", dependencies=TRUE)          #this worked
geneCategories = read.csv("./Analysis - R-Script /SubtiWiki Exports /geneCategories.csv")
myGsc.cat1 = loadGSC(cbind(as.character(geneCategories$gene),
                           as.character(geneCategories$category1)))
diffexp.pr.all.conds = list(diffexp.pr.all.LB,diffexp.pr.all.M9,diffexp.pr.all.SH2,diffexp.pr.all.SH5)
names(diffexp.pr.all.conds) = c("LB","M9","SH2","SH5")

gsea.enrich.cat1.all = c()
for (i in 1:length(diffexp.pr.all.conds)){
  my.de = diffexp.pr.all.conds[[i]]
  my.pval = my.de$adj.P.Val
  my.fc = my.de$logFC
  names(my.pval) = rownames(my.de)
  names(my.fc) = rownames(my.de)
  my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc
  my.fc = my.fc[!is.na(my.fc)]    #1707 entries 
  my.gsaRes.cat1 = runGSA(geneLevelStats=my.pval,directions=my.fc,
                          gsc=myGsc.cat1,geneSetStat="reporter",adjMethod="BH")
  my.gsea.cat1 = GSAsummaryTable(my.gsaRes.cat1)
  colnames(my.gsea.cat1) = paste0(names(diffexp.pr.all.conds)[i],"_",colnames(my.gsea.cat1))
  my.gsea.cat1 = as.matrix(my.gsea.cat1)
  gsea.enrich.cat1.all=cbind(gsea.enrich.cat1.all,my.gsea.cat1)
}

#dim(gsea.enrich.all)
#[1] 111  76

#gsea.go.sig = gsea.go[gsea.go$`p (non-dir.)` < 0.05,]
#write.csv(gsea.go, "./LB/gsea.go_LB.csv")
#write.csv(gsea.go.sig, "./LB/gsea.go.Sig_LB.csv")

rownames(gsea.enrich.cat1.all) = gsea.enrich.cat1.all[,"LB_Name"]

#write.csv(gsea.enrich.cat1.all, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/gsea.enrich.cat1.all.csv")

#gsea.enrich.pv =  gsea.enrich.all[,grep("LB_p(non-dir.)", colnames(gsea.enrich.all))]
gsea.cat1.enrich.ND =  gsea.enrich.cat1.all[,grep("(non-dir.)", colnames(gsea.enrich.cat1.all))]
colnames(gsea.cat1.enrich.ND) = c("LB_Stat(non-dir.)",
                                  "LB_p(non-dir.)",
                                  "LB_padj(non-dir.)",
                                  "M9_Stat(non-dir.)",
                                  "M9_p(non-dir.)",
                                  "M9_padj(non-dir.)",
                                  "SH2_Stat(non-dir.)",
                                  "SH2_p(non-dir.)",
                                  "SH2_padj(non-dir.)",
                                  "SH5_Stat(non-dir.)",
                                  "SH5_p(non-dir.)",
                                  "SH5_padj(non-dir.)")


colnames(gsea.cat1.enrich.ND) = c("LB_Stat_ND",
                                  "LB_p_ND",
                                  "LB_padj_ND",
                                  "M9_Stat_ND",
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

gsea.cat1.enrich.sig = gsea.cat1.enrich[(as.numeric(as.character(gsea.cat1.enrich[,"LB_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat1.enrich[,"M9_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat1.enrich[,"SH2_padj_ND"]))<0.05|
                                           as.numeric(as.character(gsea.cat1.enrich[,"SH5_padj_ND"]))<0.05),]

gsea.cat1.enrich.sig.pv = cbind(as.numeric(as.character(gsea.cat1.enrich.sig[,"LB_padj_ND"])),
                                as.numeric(as.character(gsea.cat1.enrich.sig[,"M9_padj_ND"])),
                                as.numeric(as.character(gsea.cat1.enrich.sig[,"SH2_padj_ND"])),
                                as.numeric(as.character(gsea.cat1.enrich.sig[,"SH5_padj_ND"])))


rownames(gsea.cat1.enrich.sig.pv) = rownames(gsea.cat1.enrich.sig)
colnames(gsea.cat1.enrich.sig.pv) = c("GSEA_Cat1_LB_adj.P.Val","GSEA_Cat1_M9_adj.P.Val","GSEA_Cat1_SH2_adj.P.Val","GSEA_Cat1_SH5_adj.P.Val")
gsea.cat1.enrich.sig.pv.log10 = -log10(gsea.cat1.enrich.sig.pv)
gsea.cat1.enrich.sig.pv.log10[gsea.cat1.enrich.sig.pv.log10 < -log10(0.05)] = 0 #all the values less 1.3 (that is >0.05) are white, hence all colours denote SIG pVal ONLY!

pheatmap(gsea.cat1.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10), 
         cluster_rows=TRUE, cluster_cols=FALSE,
         cellwidth = 60, cellheight = 60,
         show_rownames = T, show_colnames = T, 
         fontsize=11, legend=TRUE,
         main = "GSEA SubtiWiki Category 1 \n(adj.P.Val)- Proteomics"
)




################################################
######REGULATORS Interaction Maps
################################################


library(igraph)
regulations = read.csv("./Analysis - R-Script /SubtiWiki Exports /regulations.csv")
regulations2 = regulations[!as.character(regulations$regulator.locus.tag)==as.character(regulations$locus.tag),]  #removes auotoregulat


###############################
##REGULATORS - Universe
###############################


diffexp.pr.all.LB = read.csv("./Analysis - R-Script /Data/combined/output/diff.exp./diff.exp._minus_spoVG_upp/20190514_diffexp.pr.all_LB_minus_spoVG_upp.csv", header = T)
rownames(diffexp.pr.all.LB) = diffexp.pr.all.LB$X
diffexp.pr.all.LB = diffexp.pr.all.LB[,-1]
  

regulations2.U = regulations2[regulations2$regulator.locus.tag%in%rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05] | 
                                 regulations2$locus.tag%in%rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05], ]


#regulations2.U = regulations2 
g.U = graph.data.frame(cbind(as.character(regulations2.U$regulator.locus.tag),
                              as.character(regulations2.U$locus.tag)), directed=F)  #no direction decided at this step
g.U = igraph::simplify(g.U)  #sort of finds the unique interactions
nodes.U  = rownames(as.matrix(V(g.U))) #V is the vertex which is the node
edge.list.U = as.data.frame(get.edgelist(g.U))
degree.U = as.matrix(igraph::degree(g.U),ncol=1) #print the number of edges per vertex
degree.U = as.data.frame(degree.U)
colnames(degree.U) = "Degree"

#plot(g.U,vertex.size=3,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0,
#     mark.border="grey",
#     #vertex.color=dim(degree.U),
#     pch=19,vertex.label=NA)    #check pch 





##########################
##REGULATORS - LB - Non-directional
##########################


regulations2.LB = regulations2[regulations2$regulator.locus.tag%in%rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05] | 
  regulations2$locus.tag%in%rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05], ]

g.LB = graph.data.frame(cbind(as.character(regulations2.LB$regulator.locus.tag),
                            as.character(regulations2.LB$locus.tag)), directed=T)  #no direction decided at this step

#add gene names with the bsu numbers

Gene = read.csv("./Analysis - R-Script /SubtiWiki Exports /Gene.csv", header =T)
rownames(Gene) = rownames(Gene$locus)

g.LB = igraph::simplify(g.LB)  #sort of finds the unique interactions

nodes.LB  = rownames(as.matrix(V(g.LB))) #V is the vertex which is the node
edge.list.LB = as.data.frame(get.edgelist(g.LB))
degree.LB = as.matrix(igraph::degree(g.LB),ncol=1) #print the number of edges per vertex
degree.LB = as.data.frame(degree.LB)
colnames(degree.LB) = "Degree"


#Self exploratory 18/Jan/19

#degree.LB = degree.LB[!apply(degree.LB, function(x) is.na(degree.LB), ...)]
#degree.LB[(as.numeric(degree.LB$Degree)<500), ] = "blue"
#degree.LB[(as.numeric(degree.LB$Degree)>500), ] = "red"
#kegg.enrich.sig.pv.log10[kegg.enrich.sig.pv.log10 < -log10(0.05)] = 0
#kegg.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10)


#Self exploratory 17/Jan/19
plot(g.LB,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0.4, 
     mark.border="red",pch=19,vertex.label=degree.LB$Degree, #rownames(degree.LB[(degree.LB$Degree)>10 ,]) = gives BSU numbers in the name 
     xlim=c(-0.8,0.8), ylim = c(-0.5,1.5),
     #vertex.color=colorRampPalette(c("white", "red")) (10),
     vertex.color=degree.U$Degree,        #gives multiple colours 
     #vertex.color=dim(degree.U),  #gives two colours : orange and blue
     vertex.label.color="black", vertex.label.dist=1.5,
     #vertex.shape=diffexp.pr.all.LB_SIG$P.Value,
     vertex.shape="circle",  #how to set two shapes for two kinds of expression data
     vertex.size=degree.LB$Degree, scale=degree.LB$Degree)    #check pch #[degree.LB$Degree<300 = 800,]
#how to make the size of each each in a certain margin length eg all value >20 =20
title("Regulatory Networks in LB - Proteomics",cex.main=1,col.main="black")

#biggest dot is BSU25200(30) = sigA row 12
#second biggest  is BSU25490 (9) = hrcA row 11
#third BSU00370 (8) = abrB row 47



#Done with Tauqeer
# plot(g.LB,vertex.size=3,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0,
#       mark.border="grey",pch=19,vertex.label=NA)    #check pch


#Self exploratory 17/Sep/19
plot(g.LB,vertex.frame.color=,vertex.label.dist=0.5, edge.arrow.size=0.4, 
     mark.border="red",pch=19,vertex.label=degree.LB$Degree, #rownames(degree.LB[(degree.LB$Degree)>10 ,]) = gives BSU numbers in the name 
     xlim=c(-0.8,0.8), ylim = c(-0.5,1.5),
     #vertex.color=colorRampPalette(c("white", "red")) (10),
     vertex.color=degree.U$Degree,        #gives multiple colours 
     #vertex.color=dim(degree.U),  #gives two colours : orange and blue
     vertex.label.color="black", vertex.label.dist=1.5,
     #vertex.shape=diffexp.pr.all.LB_SIG$P.Value,
     vertex.shape="circle",  #how to set two shapes for two kinds of expression data
     vertex.size=(2.3)*(degree.LB$Degree), 
     scale=degree.LB$Degree)    #check pch #[degree.LB$Degree<300 = 800,]
#how to make the size of each each in a certain margin length eg all value >20 =20
title("Regulatory Networks in LB - Proteomics",cex.main=1,col.main="black")







































