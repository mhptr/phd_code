

setwd("~/OneDrive - University of Warwick/WORK/Results/Proteomics/FINAL Result/Analysis - R-Script /Data")


#read SIG files
diffexp.pr.all.LB_SIG = read.csv("./combined/output/diff.exp.SIG/20181217_diffexp.pr.all_LB_SIG.csv", header = T)
diffexp.pr.all.M9_SIG = read.csv("./combined/output/diff.exp.SIG/20181217_diffexp.pr.all_M9_SIG.csv", header = T)
diffexp.pr.all.SH2_SIG = read.csv("./combined/output/diff.exp.SIG/20181217_diffexp.pr.all_SH2_SIG.csv", header = T)
diffexp.pr.all.SH5_SIG = read.csv("./combined/output/diff.exp.SIG/20181217_diffexp.pr.all_SH5_SIG.csv", header = T)

#read GeneExportWizard
geneWizard = read.csv("../SubtiWiki Exports /Genes export wizard.csv", header = T)

#read GeneCategories
geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header=T)

#read GeneRegulators
geneRegulators = read.csv("../SubtiWiki Exports /regulations.csv", header=T)



#Merge Files

#LB
diffexp.pr.all.LB_SIG_Annotated = merge(diffexp.pr.all.LB_SIG, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(diffexp.pr.all.LB_SIG_Annotated) = diffexp.pr.all.LB_SIG_Annotated$X
diffexp.pr.all.LB_SIG_Annotated = diffexp.pr.all.LB_SIG_Annotated[,-1]
#write.csv(diffexp.pr.all.LB_SIG_Annotated, "./combined/output/diffexp.pr.all.LB_SIG_Annotated.csv")


diffexp.pr.all.LB_SIG_Annotated = merge(diffexp.pr.all.LB_SIG, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
diffexp.LB_SIG_Annotated_Cat3 = merge(diffexp.pr.all.LB_SIG_Annotated, geneCategories, by.x="X", by.y="gene", all.x=TRUE)
rownames(diffexp.LB_SIG_Annotated_Cat3) = diffexp.LB_SIG_Annotated_Cat3$X
diffexp.LB_SIG_Annotated_Cat3 = diffexp.LB_SIG_Annotated_Cat3[,-1]



#M9
diffexp.pr.all.M9_SIG_Annotated = merge(diffexp.pr.all.M9_SIG, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(diffexp.pr.all.M9_SIG_Annotated) = diffexp.pr.all.M9_SIG_Annotated$X
diffexp.pr.all.M9_SIG_Annotated = diffexp.pr.all.M9_SIG_Annotated[,-1]
#write.csv(diffexp.pr.all.M9_SIG_Annotated, "./combined/output/diffexp.pr.all.M9_SIG_Annotated.csv")

#SH2
diffexp.pr.all.SH2_SIG_Annotated = merge(diffexp.pr.all.SH2_SIG, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(diffexp.pr.all.SH2_SIG_Annotated) = diffexp.pr.all.SH2_SIG_Annotated$X
diffexp.pr.all.SH2_SIG_Annotated = diffexp.pr.all.SH2_SIG_Annotated[,-1]
#write.csv(diffexp.pr.all.SH2_SIG_Annotated, "./combined/output/diffexp.pr.all.SH2_SIG_Annotated.csv")

#SH5
diffexp.pr.all.SH5_SIG_Annotated = merge(diffexp.pr.all.SH5_SIG, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(diffexp.pr.all.SH5_SIG_Annotated) = diffexp.pr.all.SH5_SIG_Annotated$X
diffexp.pr.all.SH5_SIG_Annotated = diffexp.pr.all.SH5_SIG_Annotated[,-1]
#write.csv(diffexp.pr.all.SH5_SIG_Annotated, "./combined/output/diffexp.pr.all.SH5_SIG_Annotated.csv")



#Diff.pr_SH5_vs_SH2

Diff.pr_SIG_SH5_vs_SH2 = read.csv("./combined/output/diff.exp.SIG/diffexp.pr_SIG_SH5_vs_SH2.csv", header = T)
Diff.pr_SIG_SH5_vs_SH2_Annotated = merge(Diff.pr_SIG_SH5_vs_SH2, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(Diff.pr_SIG_SH5_vs_SH2_Annotated) = Diff.pr_SIG_SH5_vs_SH2_Annotated$X
Diff.pr_SIG_SH5_vs_SH2_Annotated = Diff.pr_SIG_SH5_vs_SH2_Annotated[,-1]
#write.csv(Diff.pr_SIG_SH5_vs_SH2_Annotated, "./combined/output/diff.exp.SIG/Annotated/Diff.pr_SIG_SH5_vs_SH2_Annotated.csv")




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


#########################################################################################################
## Merge Annotations to common elements
#########################################################################################################


#SH2_SH5

Common_SH2_SH5 = read.csv("./combined/output/venn/Common_SH2_SH5.csv", header = T)

Common_SH2_SH5_Annot = merge(Common_SH2_SH5, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(Common_SH2_SH5_Annot) = Common_SH2_SH5_Annot$X
Common_SH2_SH5_Annot = Common_SH2_SH5_Annot[,-1]
write.csv(Common_SH2_SH5_Annot, "./combined/output/venn/Common_SH2_SH5_Annot.csv")





