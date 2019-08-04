
#install.packages("tidyverse")

library(tidyverse)

setwd("~/OneDrive - University of Warwick/WORK/RESULTS/PROTEOMICS/FINAL Result")

geneCat = read.csv("./Analysis - R-Script /SubtiWiki Exports /geneCategories.csv", header = T)

cat2 = geneCat[,c(2,3,5)]
dim(cat2)
cat2 = unique(cat2)
dim(cat2)
write.csv(cat2, "./Analysis - R-Script /SubtiWiki Exports /split_categories/cat2.csv", row.names = F)

cat3 = geneCat[,c(2,3,6)]
cat3 = unique(cat3)
#View(data.frame(unique(cat3.only$category3)))
write.csv(cat3, "./Analysis - R-Script /SubtiWiki Exports /split_categories/cat3.csv", row.names = F)

cat4 = geneCat[,c(2,3,7)]
cat4 = unique(cat4)
write.csv(cat4, "./Analysis - R-Script /SubtiWiki Exports /split_categories/cat4.csv", row.names = F)

########

# cat2_coping_with_stress = geneCategories[geneCategories$category2 == "Coping with stress",]
#write.csv(cat2_coping_with_stress, "../SubtiWiki Exports /split_categories/cat_all/cat2_coping_with_stress.csv", row.names = F)

#Cat 2 = Coping with stress

Coping_with_stress = cat2 %>% 
  filter(category2 == "Coping with stress")
dim(Coping_with_stress)
# write.csv(Coping_with_stress, "./Analysis - R-Script /SubtiWiki Exports /split_categories/cat2/Coping_with_stress.csv", row.names = F)

General_stress_proteins = cat3 %>% 
  filter(category3 == "General stress proteins (controlled by SigB)")
dim(General_stress_proteins)
# write.csv(General_stress_proteins, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/General_stress_proteins.csv", row.names = F)

Cell_envelope_stress_proteins = cat3 %>% 
  filter(category3 == "Cell envelope stress proteins (controlled by SigM, V, W, X, Y)")
dim(Cell_envelope_stress_proteins)
#write.csv(Cell_envelope_stress_proteins, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Cell_envelope_stress_proteins.csv", row.names = F)

Acid_stress_proteins = cat3 %>%
  filter(category3 == "Acid stress proteins (controlled by YvrI-YvrHa)")
#write.csv(Acid_stress_proteins, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Acid_stress_proteins.csv", row.names = F)

Heat_shock_proteins = cat3 %>%
  filter(category3 == "Heat shock proteins")
rownames(Heat_shock_proteins) = Heat_shock_proteins$gene
#write.csv(Heat_shock_proteins, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Heat_shock_proteins.csv", row.names = F)

Cold_stress_proteins = cat3 %>%
  filter(category3 == "Cold stress proteins")
rownames(Cold_stress_proteins) = Cold_stress_proteins$gene
# write.csv(Cold_stress_proteins, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Cold_stress_proteins.csv", row.names = F)

Coping_with_hyper_osmotic_stress = cat3 %>%
  filter(category3 == "Coping with hyper-osmotic stress")
rownames(Coping_with_hyper_osmotic_stress) = Coping_with_hyper_osmotic_stress$gene
# write.csv(Coping_with_hyper_osmotic_stress, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Coping_with_hyper_osmotic_stress.csv", row.names = F)

Resistance_against_oxidative_and_electrophile_stress = cat3 %>%
  filter(category3 == "Resistance against oxidative and electrophile stress")
rownames(Resistance_against_oxidative_and_electrophile_stress) = Resistance_against_oxidative_and_electrophile_stress$gene
# write.csv(Resistance_against_oxidative_and_electrophile_stress, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Resistance_against_oxidative_and_electrophile_stress.csv", row.names = F)

Resistance_against_other_toxic_compounds = cat3 %>%
  filter(category3 == "Resistance against other toxic compounds (nitric oxide, phenolic acids, flavonoids, oxalate)")
rownames(Resistance_against_other_toxic_compounds) = Resistance_against_other_toxic_compounds$gene
# write.csv(Resistance_against_other_toxic_compounds, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Resistance_against_other_toxic_compounds.csv", row.names = F)

Resistance_against_toxic_metals = cat3 %>%
  filter(category3 == "Resistance against toxic metals")
Resistance_against_toxic_metals = rbind(Resistance_against_toxic_metals, cat3 %>%
                                          filter(category3 == "Resistance against toxic metals/ based on similarity"))
rownames(Resistance_against_toxic_metals) = Resistance_against_toxic_metals$gene
# write.csv(Resistance_against_toxic_metals, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Resistance_against_toxic_metals.csv", row.names = F)

Resistance_against_toxins = cat3 %>%
  filter(category3 == "Resistance against toxins/ antibiotics")
Resistance_against_toxins = rbind(Resistance_against_toxins, cat3 %>%
                                          filter(category3 == "Resistance against toxins/ antibiotics/ based on similarity"))
rownames(Resistance_against_toxins) = Resistance_against_toxins$gene
# write.csv(Resistance_against_toxins, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Resistance_against_toxins.csv", row.names = F)


Biosynthesis_of_antibacterial_compounds = cat3 %>%
  filter(category3 == "Biosynthesis of antibacterial compounds")
Biosynthesis_of_antibacterial_compounds  = rbind(Biosynthesis_of_antibacterial_compounds , cat3 %>%
                                    filter(category3 == "Biosynthesis of antibacterial compounds/ based on similarity"))
rownames(Biosynthesis_of_antibacterial_compounds ) = Biosynthesis_of_antibacterial_compounds $gene
# write.csv(Biosynthesis_of_antibacterial_compounds , "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Biosynthesis_of_antibacterial_compounds .csv", row.names = F)


Toxins_antitoxins_and_immunity_against_toxins = cat3 %>%
  filter(category3 == "Toxins, antitoxins and immunity against toxins")
Toxins_antitoxins_and_immunity_against_toxins = rbind(Toxins_antitoxins_and_immunity_against_toxins, cat3 %>%
                                                   filter(category3 == "Toxins, antitoxins and immunity against toxins/ based on similarity"))
rownames(Toxins_antitoxins_and_immunity_against_toxins) = Toxins_antitoxins_and_immunity_against_toxins$gene
# write.csv(Toxins_antitoxins_and_immunity_against_toxins   , "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/cat2_coping_with_stress/Toxins_antitoxins_and_immunity_against_toxins.csv", row.names = F)


#Cat 3 = Biofilm formation
Biofilm_formation   = cat3 %>%
  filter(category3 == "Biofilm formation")
Biofilm_formation = cbind(Biofilm_formation, cat4$category4)
rownames(Biofilm_formation  ) = Biofilm_formation$gene
# write.csv(Biofilm_formation, "./Analysis - R-Script /SubtiWiki Exports /split_categories/conditions_cat3/Biofilm_formation.csv", row.names = F)


# cat3_Biofilm_formation = geneCategories[geneCategories$category3 == "Biofilm formation",]
# write.csv(cat3_Biofilm_formation, "../SubtiWiki Exports /split_categories/cat3/cat3_Biofilm_formation.csv", row.names = F)






















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

Swarming = cat3 %>%
  filter(category3 == "Swarming")
dim(Swarming)

Motility_and_chemotaxis = cat3 %>%
  filter(category3 == "Motility and chemotaxis")
rownames(Motility_and_chemotaxis) = Motility_and_chemotaxis$gene
Motility_and_chemotaxis = Motility_and_chemotaxis[,-1]
dim(Motility_and_chemotaxis)

Flagellar_proteins = cat4 %>%
  filter(category4 == "Flagellar proteins")
rownames(Flagellar_proteins) = Flagellar_proteins$gene
Flagellar_proteins = Flagellar_proteins[,-1]
dim(Flagellar_proteins)

























# cat3.enriched = c ("General stress proteins (controlled by SigB)",
#                    "Miscellaneous metabolic pathways",
#                    "Biosynthesis of antibacterial compounds",
#                    "Sporulation proteins",
#                    "Translation",
#                    "PBSX prophag",
#                    "Utilization of specific carbon sources",
#                    "Heat shock proteins")
# cat3.paths = c()
# for (i in 1:length(cat3.enriched)){
#   my.de = cat3.enriched[[i]]
#   my.path.gene = cat3 %>%
#     filter(category3 == cat3.enriched )
#   colnames(my.path.gene) = cat3.enriched
#   my.path.gene = as.matrix(my.path.gene)
#   cat3.paths = rbind(cat3.paths, my.path.gene)
# }
# 
#   
# my.pval = my.de$adj.P.Val
#   my.fc = my.de$logFC
#   names(my.pval) = rownames(my.de)
#   names(my.fc) = rownames(my.de)
#   my.pval = my.pval[!is.na(my.fc)]  #changing the number of pVal based on number of fc
#   my.fc = my.fc[!is.na(my.fc)]    #1707 entries 
#   my.gsaRes.cat3 = runGSA(geneLevelStats=my.pval,directions=my.fc,
#                           gsc=myGsc.cat3,geneSetStat="reporter",adjMethod="BH")
#   my.gsea.cat3 = GSAsummaryTable(my.gsaRes.cat3)
#   colnames(my.gsea.cat3) = paste0(names(diffexp.pr.all.conds)[i],"_",colnames(my.gsea.cat3))
#   my.gsea.cat3 = as.matrix(my.gsea.cat3)
#   gsea.enrich.cat3.all=cbind(gsea.enrich.cat3.all,my.gsea.cat3)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# HSP$gene
# 
# 
# LB_Sig = read.csv("./Analysis - R-Script /Data/combined/output/diff.exp.SIG/20181217_diffexp.pr.all_LB_SIG.csv", header = T)
# LB_Sig_Cat3 = merge(LB_Sig, cat3.only, by.x="X", by.y="gene")
# LB_Sig_Cat3 = merge(cat3.only, LB_Sig, by.x="gene", by.y="X", by.x=T)
# 
# #gsea.cat3.enrich.ND =  gsea.enrich.cat3.all[,grep("(non-dir.)", colnames(gsea.enrich.cat3.all))]
# 
# 
# 
# as.list(cat3.only$category3["Heat shock protein"])














