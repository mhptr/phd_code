#Adding annotations to DPP 
# 03/08/2019

setwd("~/OneDrive - University of Warwick/WORK/Results/Proteomics/FINAL Result/Analysis - R-Script /Data")

#tring a new branch

# #read SIG files minus spoVG minus upp
# 
# DPP_LB = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/06062019_diffexp.pr.all_LB_SIG_minus_spoVG_upp.csv", header = T)
# DPP_M9 = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/06062019_diffexp.pr.all_M9_SIG_minus_spoVG_upp.csv", header = T)
# DPP_SH2 = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/03082019_diffexp.pr.all_SH2_SIG_minus_spoVG_upp.csv", header = T)
# DPP_SH5 = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/03082019_diffexp.pr.all_SH5_SIG_minus_spoVG_upp.csv", header = T)
# DPP_SH5_vs_SH2 = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/diffexp.pr_SIG_WT_SH5_vs_SH2_minus_spoVG_upp.csv", header = T)
# 
# #read SubtiWiki files 
# geneWizard = read.csv("../SubtiWiki Exports /Genes export wizard.csv", header = T)
# geneWizard = geneWizard[,(1:12)]
# 
# #Merge Files
# 
# DPP_LB_annot = merge(DPP_LB, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
# DPP_M9_annot = merge(DPP_M9, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
# DPP_SH2_annot = merge(DPP_SH2, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
# DPP_SH5_annot = merge(DPP_SH5, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
# DPP_SH5_vs_SH2_annot = merge(DPP_SH5_vs_SH2, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
# 
# # write.csv(DPP_LB_annot, "./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_LB_minus_spoVG_upp_annot.csv", row.names = F)
# # write.csv(DPP_M9_annot, "./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_M9_minus_spoVG_upp_annot.csv", row.names = F)
# # write.csv(DPP_SH2_annot, "./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_SH2_minus_spoVG_upp_annot.csv", row.names = F)
# # write.csv(DPP_SH5_annot, "./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_SH5_minus_spoVG_upp_annot.csv", row.names = F)
# # write.csv(DPP_SH5_vs_SH2_annot, "./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_SH5_vs_SH2_annot_minus_spoVG_upp_annot.csv", row.names = F)


#########################################################################################################
## Merge Cat2, Cat3 and Cat4 with proteomic data = 03/08/19
#########################################################################################################

geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header = T)
# cat2_coping_with_stress = geneCategories[geneCategories$category2 == "Coping with stress",]
# write.csv(cat2_coping_with_stress, "../SubtiWiki Exports /split_categories/cat_all/cat2_coping_with_stress.csv", row.names = F)

# cat3_Biofilm_formation = geneCategories[geneCategories$category3 == "Biofilm formation",]
# write.csv(cat3_Biofilm_formation, "../SubtiWiki Exports /split_categories/cat3/cat3_Biofilm_formation.csv", row.names = F)






cat2_coping_with_stress = read.csv("../SubtiWiki Exports /split_categories/cat_all/cat2_coping_with_stress.csv", header = T)
cat3_Biofilm_formation = read.csv("../SubtiWiki Exports /split_categories/cat3/cat3_Biofilm_formation.csv", header = T)



DPP_LB = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_LB_minus_spoVG_upp_annot.csv", header = T)
DPP_M9 = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_M9_minus_spoVG_upp_annot.csv", header = T)
DPP_SH2 = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_SH2_minus_spoVG_upp_annot.csv", header = T)
DPP_SH5 = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_SH5_minus_spoVG_upp_annot.csv", header = T)
DPP_SH5_vs_SH2 = read.csv("./combined/output/diff.exp.SIG/data/diff.exp.SIG_minus_spoVG_upp/Annotated/DPP_SH5_vs_SH2_annot_minus_spoVG_upp_annot.csv", header = T)

#DPP

DPP_LB_cat2_coping_with_stress = merge(DPP_LB, cat2_coping_with_stress, by.x="X", by.y="gene", all.x = TRUE)
dim(DPP_LB_cat2_coping_with_stress)
# write.csv(DPP_LB_cat2_coping_with_stress, "./combined/output/GSEA_SGC/DPP/LB/DPP_LB_cat2_coping_with_stress.csv", row.names = F)
DPP_M9_cat2_coping_with_stress = merge(DPP_M9, cat2_coping_with_stress, by.x="X", by.y="gene", all.x = TRUE)
dim(DPP_M9_cat2_coping_with_stress)
# write.csv(DPP_M9_cat2_coping_with_stress, "./combined/output/GSEA_SGC/DPP/M9/DPP_M9_cat2_coping_with_stress.csv", row.names = F)
#want to all the columns and rows with no valid values
#write.csv(DPP_M9_cat2_coping_with_stress, "./combined/output/GSEA_SGC/DPP/M9/DPP_M9_cat2_coping_with_stress.csv", row.names = F,  na="")
DPP_SH2_cat2_coping_with_stress = merge(DPP_SH2, cat2_coping_with_stress, by.x="X", by.y="gene", all.x = TRUE)
dim(DPP_SH2_cat2_coping_with_stress)
# write.csv(DPP_SH2_cat2_coping_with_stress, "./combined/output/GSEA_SGC/DPP/SH2/DPP_SH2_cat2_coping_with_stress.csv", row.names = F)
DPP_SH5_cat2_coping_with_stress = merge(DPP_SH5, cat2_coping_with_stress, by.x="X", by.y="gene", all.x = TRUE)
dim(DPP_SH5_cat2_coping_with_stress)
# write.csv(DPP_SH5_cat2_coping_with_stress, "./combined/output/GSEA_SGC/DPP/SH5/DPP_LB_cat2_coping_with_stress.csv", row.names = F)
DPP_SH5_vs_SH2_cat2_coping_with_stress = merge(DPP_SH5_vs_SH2, cat2_coping_with_stress, by.x="X", by.y="gene", all.x = TRUE)
dim(DPP_SH5_vs_SH2_cat2_coping_with_stress)
# write.csv(DPP_SH5_vs_SH2_cat2_coping_with_stress, "./combined/output/GSEA_SGC/DPP/SH5_vs_SH2/DPP_SH5_vs_SH2_cat2_coping_with_stress.csv", row.names = F)


DPP_M9_cat3_Biofilm_formation = merge(DPP_M9, cat3_Biofilm_formation, by.x="X", by.y="gene", all.x = TRUE)
dim(DPP_M9_cat3_Biofilm_formation)
# write.csv(DPP_M9_cat3_Biofilm_formation, "./combined/output/GSEA_SGC/DPP/M9/DPP_M9_cat3_Biofilm_formation.csv", row.names = F)



##########################################################################################
##########################################################################################
##########################################################################################



# Coping_with_stress_cat2 = read.csv("../SubtiWiki Exports /split_categories/cat2/Coping_with_stress.csv", header = T)
# General_stress_proteins = read.csv("../SubtiWiki Exports /split_categories/cat3/cat2_coping_with_stress/General_stress_proteins.csv", header = T)
# Cell_envelope_stress_proteins
# Acid_stress_proteins
# Heat_shock_proteins
# Cold_stress_proteins
# Coping_with_hyper_osmotic_stress
# Resistance_against_oxidative_and_electrophile_stress
# Resistance_against_other_toxic_compounds
# Resistance_against_toxic_metals
# Resistance_against_toxins
# Biosynthesis_of_antibacterial_compounds
# Toxins_antitoxins_and_immunity_against_toxins
# Biofilm_formation




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
# write.csv(Common_SH2_SH5_Annot, "./combined/output/venn/Common_SH2_SH5_Annot.csv")













