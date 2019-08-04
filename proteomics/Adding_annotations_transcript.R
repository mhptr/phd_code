#This scripts adds different annotations by merging .csv files


setwd("~/OneDrive - University of Warwick/WORK/RESULTS/TRANSCRIPTOMICS/intersection_toRNAdo_all/diff.exp.gene/")

geneWizard = read.csv("../SubtiWiki Exports /Genes export wizard.csv", header = T)
geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header=T)
geneRegulators = read.csv("../SubtiWiki Exports /regulations.csv", header=T)

################################
## ELD - LB - ko vs wt 
################################

DEG_ko_LB = read.csv("./DEG_ED/DEG_ko_vs_wt_p0.05_LB_eld.csv", header = T)
DEG_ko_LB_Annot = merge(DEG_ko_LB, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_ko_LB_Annot) = DEG_ko_LB_Annot$X
DEG_ko_LB_Annot = DEG_ko_LB_Annot[,-1]
#write.csv(diffexp.pr.all.LB_SIG_Annotated, "./combined/output/diffexp.pr.all.LB_SIG_Annotated.csv")

################################
## ELD - LB - dc vs wt 
################################

DEG_dc_LB = read.csv("./DEG_ED/DEG_dc_vs_wt_p0.05_LB_eld.csv", header = T)
DEG_dc_LB_Annot = merge(DEG_dc_LB, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_dc_LB_Annot) = DEG_dc_LB_Annot$X
DEG_dc_LB_Annot = DEG_dc_LB_Annot[,-1]
#write.csv(diffexp.pr.all.LB_SIG_Annotated, "./combined/output/diffexp.pr.all.LB_SIG_Annotated.csv")



###################################
###################################   Merging the Cat3 and Diff_Sig names s
###################################

setwd("~/OneDrive - University of Warwick/WORK/RESULTS/PROTEOMICS/FINAL Result/Analysis - R-Script /Data")


#################
### M9
#################

Heat_shock_proteins = read.csv("./combined/output/GSEA_SGC/Heat_shock_proteins.csv", header = T)
Heat_shock_proteins = Heat_shock_proteins[,-1]

diffexp.pr.all.M9_SIG_Annotated = read.csv("./combined/output/diffexp.pr.all.M9_SIG_Annotated.csv", header = T)
Heat_shock_proteins_M9 = merge(diffexp.pr.all.M9_SIG_Annotated, Heat_shock_proteins, by.x="X", by.y="gene", all.y = TRUE)
rownames(Heat_shock_proteins_M9) = Heat_shock_proteins_M9$X
colnames(Heat_shock_proteins_M9) = colnames(Heat_shock_proteins_M9)

Heat_shock_proteins_M9 = Heat_shock_proteins_M9[,cbind("name", "description",
                                         "product", "essential", "names", 
                                         "logFC", "adj.P.Val")]

#write.csv(Heat_shock_proteins_M9, "./combined/output/GSEA_SGC/Heat_shock_proteins_M9.csv")



#################
### SH2
#################

Translation = read.csv("./combined/output/GSEA_SGC/Translation.csv", header = T)
Translation = Translation[,-1]

diffexp.pr.all.SH2_SIG_Annotated = read.csv("./combined/output/diffexp.pr.all.SH2_SIG_Annotated.csv", header = T)
Translation_SH2 = merge(diffexp.pr.all.SH2_SIG_Annotated, Translation, by.x="X", by.y="gene", all.y = TRUE)
rownames(Translation_SH2) = Translation_SH2$X
colnames(Translation_SH2) = colnames(Translation_SH2)

Translation_SH2 = Translation_SH2[,cbind("name", "description",
                                                           "product", "essential", "names", 
                                                           "logFC", "adj.P.Val")]

#write.csv(Translation_SH2, "./combined/output/GSEA_SGC/Translation_SH2.csv")




#################
### SH5
#################

Sporulation_proteins = read.csv("./combined/output/GSEA_SGC/Sporulation_proteins.csv", header = T)
Sporulation_proteins = Sporulation_proteins[,-1]

diffexp.pr.all.SH5_SIG_Annotated = read.csv("./combined/output/diffexp.pr.all.SH5_SIG_Annotated.csv", header = T)
Sporulation_proteins_SH5 = merge(diffexp.pr.all.SH5_SIG_Annotated, Sporulation_proteins, by.x="X", by.y="gene", all.y = TRUE)
rownames(Sporulation_proteins_SH5) = Sporulation_proteins_SH5$X
colnames(Sporulation_proteins_SH5) = colnames(Sporulation_proteins_SH5)

Sporulation_proteins_SH5 = Sporulation_proteins_SH5[,cbind("name", "description",
                                                           "product", "essential", "names", 
                                                           "logFC", "adj.P.Val")]

#write.csv(Sporulation_proteins_SH5, "./combined/output/GSEA_SGC/Sporulation_proteins_SH5.csv")
#Sporulation_proteins_SH5 = Sporulation_proteins_SH5[!is.na(Sporulation_proteins_SH5)>2,]










