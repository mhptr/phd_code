#This scripts adds different annotations by merging .csv files


setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/diff.exp.gene/")

geneWizard = read.csv("../SubtiWiki Exports /Genes export wizard.csv", header = T)
geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header=T)
geneRegulators = read.csv("../SubtiWiki Exports /regulations.csv", header=T)

################################
## SM - M9 - ko vs wt 
################################

DEG_ko_M9 = read.csv("./DEG_SM/DEG_ko_vs_wt_p0.05_M9_sm.csv", header = T)
DEG_ko_M9_Annot = merge(DEG_ko_M9, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_ko_M9_Annot) = DEG_ko_M9_Annot$X
DEG_ko_M9_Annot = DEG_ko_M9_Annot[,-1]
#write.csv(DEG_ko_M9_Annot, "./DEG_SM/Annotated_sm/DEG_ko_M9_Annot_sm.csv")

DEG_ko_M9_Annot = read.csv("./DEG_SM/Annotated_sm/DEG_ko_M9_Annot_sm.csv", header = T)
DEG_ko_M9_Annot_Reg = merge(DEG_ko_M9_Annot, geneRegulators, by.x="X", by.y="locus.tag")
#DEG_ko_M9_Annot_Reg = merge(DEG_ko_M9_Annot, geneRegulators, by.x="X", by.y="regulator.locus.tag")
rownames(DEG_ko_M9_Annot_Reg) = DEG_ko_M9_Annot_Reg$X
DEG_ko_M9_Annot_Reg = DEG_ko_M9_Annot_Reg[,-1]
#write.csv(DEG_ko_M9_Annot_Reg, "./DEG_SM/Annotated_sm/DEG_ko_M9_Annot_Reg_sm.csv")


################################
## SM - SH2 - ko vs wt 
################################


DEG_ko_SH2 = read.csv("./DEG_SM/DEG_ko_vs_wt_p0.05_SH2_sm.csv", header = T)
DEG_ko_SH2_Annot = merge(DEG_ko_SH2, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_ko_SH2_Annot) = DEG_ko_SH2_Annot$X
DEG_ko_SH2_Annot = DEG_ko_SH2_Annot[,-1]
#write.csv(DEG_ko_SH2_Annot, "./DEG_SM/Annotated_sm/DEG_ko_SH2_Annot_sm.csv")

DEG_ko_SH2_Annot = read.csv("./DEG_SM/Annotated_sm/DEG_ko_SH2_Annot_sm.csv", header = T)
DEG_ko_SH2_Annot_Reg = merge(DEG_ko_SH2_Annot, geneRegulators, by.x="X", by.y="locus.tag")
#DEG_ko_SH2_Annot_Reg = merge(DEG_ko_SH2_Annot, geneRegulators, by.x="X", by.y="regulator.locus.tag")
rownames(DEG_ko_SH2_Annot_Reg) = DEG_ko_SH2_Annot_Reg$X
DEG_ko_SH2_Annot_Reg = DEG_ko_SH2_Annot_Reg[,-1]
#write.csv(DEG_ko_SH2_Annot_Reg, "./DEG_SM/DEG_ko_SH2_Annot_Reg_sm.csv")


################################
## SM - SH5 - ko vs wt 
################################


DEG_ko_SH5 = read.csv("./DEG_SM/DEG_ko_vs_wt_p0.05_SH5_sm.csv", header = T)
DEG_ko_SH5_Annot = merge(DEG_ko_SH5, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_ko_SH5_Annot) = DEG_ko_SH5_Annot$X
DEG_ko_SH5_Annot = DEG_ko_SH5_Annot[,-1]
#write.csv(DEG_ko_SH5_Annot, "./DEG_SM/Annotated_sm/DEG_ko_SH5_Annot_sm.csv")

DEG_ko_SH5_Annot = read.csv("./DEG_SM/Annotated_sm/DEG_ko_SH5_Annot_sm.csv", header = T)
DEG_ko_SH5_Annot_Reg = merge(DEG_ko_SH5_Annot, geneRegulators, by.x="X", by.y="locus.tag")
#DEG_ko_SH5_Annot_Reg = merge(DEG_ko_SH5_Annot, geneRegulators, by.x="X", by.y="regulator.locus.tag")
rownames(DEG_ko_SH5_Annot_Reg) = DEG_ko_SH5_Annot_Reg$X
DEG_ko_SH5_Annot_Reg = DEG_ko_SH5_Annot_Reg[,-1]
#write.csv(DEG_ko_SH5_Annot_Reg, "./DEG_SM/Annotated_sm/DEG_ko_SH5_Annot_Reg_sm.csv")



################################################################
## 
################################################################



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










