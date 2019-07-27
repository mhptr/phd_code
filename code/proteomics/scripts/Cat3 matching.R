

#install.packages("tidyverse")

library(tidyverse)

setwd("~/OneDrive - University of Warwick/WORK/RESULTS/PROTEOMICS/FINAL Result")

geneCat = read.csv("./Analysis - R-Script /SubtiWiki Exports /geneCategories.csv", header = T)
cat3 = geneCat[,c(2,3,6)]
cat3 = unique(cat3)
#View(data.frame(unique(cat3.only$category3)))

Heat_shock_proteins = cat3 %>%
  filter(category3 == "Heat shock proteins") 


General_stress_proteins = cat3 %>% 
filter(category3 == "General stress proteins (controlled by SigB)")
dim(General_stress_proteins)

Miscellaneous_metabolic_pathways = cat3 %>%
  filter(category3 == "Miscellaneous metabolic pathways")
dim(Miscellaneous_metabolic_pathways)

Biosynthesis_of_antibacterial_compounds = cat3 %>%
  filter(category3 == "Biosynthesis of antibacterial compounds")
dim(Biosynthesis_of_antibacterial_compounds)

Sporulation_proteins = cat3 %>%
  filter(category3 == "Sporulation proteins")
dim(Sporulation_proteins)


Translation = cat3 %>%
  filter(category3 == "Translation")
dim(Translation)



#write.csv(Sporulation_proteins, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/Sporulation_proteins.csv")
#write.csv(Translation, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/Translation.csv")
#write.csv(Heat_shock_proteins, "./Analysis - R-Script /Data/combined/output/GSEA_SGC/Heat_shock_proteins.csv")



cat3.enriched = c ("General stress proteins (controlled by SigB)",
                   "Miscellaneous metabolic pathways",
                   "Biosynthesis of antibacterial compounds",
                   "Sporulation proteins",
                   "Translation",
                   "PBSX prophag",
                   "Utilization of specific carbon sources",
                   "Heat shock proteins")
cat3.paths = c()
for (i in 1:length(cat3.enriched)){
  my.de = cat3.enriched[[i]]
  my.path.gene = cat3 %>%
    filter(category3 == cat3.enriched )
  colnames(my.path.gene) = cat3.enriched
  my.path.gene = as.matrix(my.path.gene)
  cat3.paths = rbind(cat3.paths, my.path.gene)
}

  




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











HSP$gene


LB_Sig = read.csv("./Analysis - R-Script /Data/combined/output/diff.exp.SIG/20181217_diffexp.pr.all_LB_SIG.csv", header = T)
LB_Sig_Cat3 = merge(LB_Sig, cat3.only, by.x="X", by.y="gene")
LB_Sig_Cat3 = merge(cat3.only, LB_Sig, by.x="gene", by.y="X", by.x=T)

#gsea.cat3.enrich.ND =  gsea.enrich.cat3.all[,grep("(non-dir.)", colnames(gsea.enrich.cat3.all))]



as.list(cat3.only$category3["Heat shock protein"])














