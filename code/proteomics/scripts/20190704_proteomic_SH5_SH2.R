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
####################################################################
### 5. Annotate Data  
####################################################################


data.pr = data.pr.all[,1:24] #[Row,Column]
annotation.all = data.pr.all[,25:dim(data.pr.all)[2]] #[2] means the columns
#Correlation Plots 
plot(data.pr[,1], data.pr[,2])

dim(data.pr) 

data.pr.v2 = data.pr[!apply(data.pr, 1, function(x) any(is.na(x))), ]


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

####################################################################
### 8. Assign Row Names 
####################################################################

rownames(data.pr.imputed) = annotation.all$T..ENSG   #Proteins with 1 BSU number only


####################################################################
### 10. Set up the Data Frame   
####################################################################


#### 10.1. Design
design <- model.matrix(~0+sample.details.imputed$media_genot)


#colnames(design) <- levels(conditions)
colnames(design) = gsub("sample.details.imputed\\$media_genot","",colnames(design))  #"\\$" maked the dollar sign as non-special character
rownames(design) <- colnames(data.pr.imputed)

####################################################################
### 5. Set REFERENCE and Making CONTRAST 
####################################################################

####################################################################
### 5.1 Del vs WT 
####################################################################


# mc = makeContrasts(
#   de_KO_vs_WT_in_LB = LB.D-LB.W,
#   de_KO_vs_WT_in_M9 = M9.D-M9.W,
#   de_KO_vs_WT_in_SH2 = SH2.D-SH2.W,
#   de_KO_vs_WT_in_SH5 = SH5.D-SH5.W,
#   #diff = (((LB.D-SH2.D)-M9.D)-SH5.D) - (((LB.W-SH2.W)-M9.W)-SH5.W),
#   levels = design
# )

mc = makeContrasts(
  de_WT_vs_WT_in_SH5 = SH5.W-SH2.W,
  levels = design
)


### MAKING CONTRAST
lm.fit <- lmFit(data.pr.imputed, design)
c.fit = contrasts.fit(lm.fit, mc)
eb = eBayes(c.fit)


####################################################################
### 5. DIFFRENTIAL EXPRESSION (WT SH5 vs WT SH2)
####################################################################
#### 5. differentiually expressed proteins WT SH5 referene to WT SH2
diffexp.pr.all.SH5_vs_SH2 = topTable(eb, adjust="BH",coef = 1,number = dim(data.pr.imputed)[1])

##################
##  SIGNIFICANT
##################

diffexp.pr_SIG_SH5_vs_SH2 = diffexp.pr.all.SH5_vs_SH2[(diffexp.pr.all.SH5_vs_SH2$adj.P.Val <0.05),]
#write.csv(diffexp.pr_SIG_SH5_vs_SH2, "./Analysis - R-Script /Data/combined/output/diff.exp.SIG/diffexp.pr_SIG_SH5_vs_SH2.csv")



