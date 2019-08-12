#This scripts adds different annotations by merging .csv files


setwd("~/OneDrive - University of Warwick/WORK/Results/Trans_Prot_COMBINED/Data/")

geneWizard = read.csv("./SubtiWiki Exports /Genes export wizard.csv", header = T)
geneCategories = read.csv("./SubtiWiki Exports /geneCategories.csv", header=T)
geneRegulators = read.csv("./SubtiWiki Exports /regulations.csv", header=T)

################################
## Trans_Prot_COMBINED_M9
################################

Common_M9 = read.csv("./Result/M9_common_imput_2.csv", header = T)
Common_M9_Annot = merge(Common_M9, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(Common_M9_Annot) = Common_M9_Annot$X
Common_M9_Annot = Common_M9_Annot[,-1]
#write.csv(Common_M9_Annot, "./Result_Annotated/Common_M9_Annot.csv")

################################
## Trans_Prot_COMBINED_SH2
################################

Common_SH2 = read.csv("./SH2_common_imput_2.csv", header = T)
Common_SH2_Annot = merge(Common_SH2, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(Common_SH2_Annot) = Common_SH2_Annot$X
Common_SH2_Annot = Common_SH2_Annot[,-1]
#write.csv(Common_SH2_Annot, "./Result_Annotated/Common_SH2_Annot.csv")

################################
## Trans_Prot_COMBINED_SH5
################################

Common_SH5 = read.csv("./SH5_common_imput_2.csv", header = T)
Common_SH5_Annot = merge(Common_SH5, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(Common_SH5_Annot) = Common_SH5_Annot$X
Common_SH5_Annot = Common_SH5_Annot[,-1]
#write.csv(Common_SH5_Annot, "./Result_Annotated/Common_SH5_Annot.csv")

################################
## Adding annotation - 12/07/19
################################

setwd("~/OneDrive - University of Warwick/WORK/Results/Trans_Prot_COMBINED/Data/")

geneWizard = read.csv("./SubtiWiki Exports /Genes export wizard.csv", header = T)
geneCategories = read.csv("./SubtiWiki Exports /geneCategories.csv", header=T)
geneRegulators = read.csv("./SubtiWiki Exports /regulations.csv", header=T)

DEG_M9 = read.csv("./DEG_SIG_minus_spoVG_upp/DEG_ko_M9_Annot_minus_spoVG_upp.csv", header=T)
rownames(DEG_M9) = DEG_M9$X
DEG_M9 = DEG_M9[,-1]

DPP_M9 = read.csv("./DPP_SIG_minus_spoVG_upp/04062019_diffexp.pr.all_M9_SIG_minus_spoVG_upp.csv", header=T)
rownames(DEG_M9) = DEG_M9$X
DEG_M9 = DEG_M9[,-1]


#M9

intersect_TP_M9 = read.csv("./Result/intersect_TP_M9.csv", header = T)
intersect_TP_M9_Annot = merge(intersect_TP_M9, geneWizard, by.x="a3", by.y="locus", all.x = TRUE)
rownames(intersect_TP_M9_Annot) = intersect_TP_M9_Annot$a3
intersect_TP_M9_Annot = intersect_TP_M9_Annot[,-1]
# write.csv(intersect_TP_M9_Annot, "./Result/Result_Annotated/intersect_TP_M9_Annot.csv")
intersect_TP_M9_Annot = merge(intersect_TP_M9_Annot, DEG_M9)

#SH2

intersect_TP_SH2 = read.csv("./Result/intersect_TP_SH2.csv", header = T)
intersect_TP_SH2_Annot = merge(intersect_TP_SH2, geneWizard, by.x="a3", by.y="locus", all.x = TRUE)
rownames(intersect_TP_SH2_Annot) = intersect_TP_SH2_Annot$a3
intersect_TP_SH2_Annot = intersect_TP_SH2_Annot[,-1]
#write.csv(intersect_TP_SH2_Annot, "./Result/Result_Annotated/intersect_TP_SH2_Annot.csv")


#SH5

intersect_TP_SH5 = read.csv("./Result/intersect_TP_SH5.csv", header = T)
intersect_TP_SH5_Annot = merge(intersect_TP_SH5, geneWizard, by.x="a3", by.y="locus", all.x = TRUE)
rownames(intersect_TP_SH5_Annot) = intersect_TP_SH5_Annot$a3
intersect_TP_SH5_Annot = intersect_TP_SH5_Annot[,-1]
#write.csv(intersect_TP_SH5_Annot, "./Result/Result_Annotated/intersect_TP_SH5_Annot.csv")




#####################


## 16/07/19


venn_intersect_annotate = function(DEG,DPP, tag) {
  
  DegFileLocation = paste("./DEG_SIG_minus_spoVG_upp/", DEG, sep = "")
  DppFileLocation = paste("./DPP_SIG_minus_spoVG_upp/", DPP, sep = "")
  DEG_File = read.csv(DegFileLocation, header = T)
  DPP_File = read.csv(DppFileLocation, header = T)
  rownames(DEG_File) = DEG_File$X
  DEG_File = DEG_File[,-1]
  rownames(DPP_File) = DPP_File$X
  DPP_File = DPP_File[,-1]
  
  # Make Venn
  venn.plot = venn.diagram(
    list("DEG" = rownames(DEG_File),
         "DPP" = rownames(DPP_File)),
    filename = NULL,
    col="transparent",
    fill=swati.color[7:8], height = 3000, width = 3000,resolution =500, imagetype = "tiff",
    main =  paste("Venn Diagram - Common Elements", tag, sep = " "),
    main.pos  = c(0.5, 1.05), main.fontface = "bold",
    main.fontfamily = "Helvetica", main.col = "Dark Blue", 
    main.cex = 1.5, cat.cex = 2, cat.fontfamily = "Helvetica", 
    cex = 2.3, fontfamily = "Helvetica",fontface = "bold",
    cat.pos = 0, margin = 0) #cat.dist = 0.07)
  grid.newpage()
  grid.draw(venn.plot)
  
  # Make overlap
  DEG_File_rowname = rownames(DEG_File)
  DPP_File_rowname = rownames(DPP_File)
  conditions = list(DEG_File_rowname, DPP_File_rowname)
  common = calculate.overlap(conditions)
  
  #intersect DEG and DPP
  intersect_TP = common["a3"]
  write.csv(intersect_TP, paste("./Result/Test/intersect_TP_", tag, ".csv", sep = ""),row.names=FALSE)
  DEG_minus_intersect_TP = common["a1"]
  write.csv(DEG_minus_intersect_TP, paste("./Result/Test/DEG_minus_intersect_TP_", tag, ".csv", sep = ""),row.names=FALSE)
  DPP_minus_intersect_TP = common["a2"]
  write.csv(DPP_minus_intersect_TP, paste("./Result/Test/DPP_minus_intersect_TP_", tag, ".csv", sep = ""),row.names=FALSE)
  
  # Add proteomics and transcriptomics data
  
  intersect_TP = merge(intersect_TP, geneWizard, by.x="a3", by.y="locus", all.x = TRUE)
  rownames(intersect_TP_SH5_Annot) = intersect_TP_SH5_Annot$a3
  intersect_TP_SH5_Annot = intersect_TP_SH5_Annot[,-1]
}


DEG_M9 = read.csv("./DEG_SIG_minus_spoVG_upp/DEG_ko_M9_Annot_minus_spoVG_upp.csv", header = T)
rownames(DEG_M9) = DEG_M9$X
DPP_M9 = read.csv("./DPP_SIG_minus_spoVG_upp/04062019_diffexp.pr.all_M9_SIG_minus_spoVG_upp.csv", header = T)
rownames(DPP_M9) = DPP_M9$X
intersect_TP_M9 = read.csv("./Result/Result_Annotated/intersect_TP_M9_Annot.csv", header = T)
rownames(intersect_TP_M9) = intersect_TP_M9$X
intersect_TP_M9 = merge(intersect_TP_M9, DEG_M9, by.x="X", by.y="X", all.x = TRUE)
intersect_TP_M9 = merge(intersect_TP_M9, DPP_M9, by.x="X", by.y="X", all.x = TRUE)
intersect_TP_M9 = intersect_TP_M9[,cbind(13:18,30:35,1:12)]
rownames(intersect_TP_M9) = intersect_TP_M9$X
# write.csv(intersect_TP_M9, "./Result/Result_Annotated/Result_TP_Annotated/intersect_TP_M9.csv")



DEG_SH2 = read.csv("./DEG_SIG_minus_spoVG_upp/DEG_ko_SH2_Annot_minus_spoVG_upp.csv", header = T)
rownames(DEG_SH2) = DEG_SH2$X
DPP_SH2 = read.csv("./DPP_SIG_minus_spoVG_upp/04062019_diffexp.pr.all_SH2_SIG_minus_spoVG_upp.csv", header = T)
rownames(DPP_SH2) = DPP_SH2$X
intersect_TP_SH2 = read.csv("./Result/Result_Annotated/intersect_TP_SH2_Annot.csv", header = T)
rownames(intersect_TP_SH2) = intersect_TP_SH2$X
intersect_TP_SH2 = merge(intersect_TP_SH2, DEG_SH2, by.x="X", by.y="X", all.x = TRUE)
intersect_TP_SH2 = merge(intersect_TP_SH2, DPP_SH2, by.x="X", by.y="X", all.x = TRUE)
intersect_TP_SH2 = intersect_TP_SH2[,cbind(13:18,30:35,1:12)]
rownames(intersect_TP_SH2) = intersect_TP_SH2$X
# write.csv(intersect_TP_SH2, "./Result/Result_Annotated/Result_TP_Annotated/intersect_TP_SH2.csv")


DEG_SH5 = read.csv("./DEG_SIG_minus_spoVG_upp/DEG_ko_SH5_Annot_minus_spoVG_upp.csv", header = T)
rownames(DEG_SH5) = DEG_SH5$X
DPP_SH5 = read.csv("./DPP_SIG_minus_spoVG_upp/04062019_diffexp.pr.all_SH5_SIG_minus_spoVG_upp.csv", header = T)
rownames(DPP_SH5) = DPP_SH5$X

intersect_TP_SH5 = read.csv("./Result/Result_Annotated/intersect_TP_SH5_Annot.csv", header = T)
rownames(intersect_TP_SH5) = intersect_TP_SH5$X
intersect_TP_SH5 = merge(intersect_TP_SH5, DEG_SH5, by.x="X", by.y="X", all.x = TRUE, no.dups = TRUE)
intersect_TP_SH5 = merge(intersect_TP_SH5, DPP_SH5, by.x="X", by.y="X", all.x = TRUE)
intersect_TP_SH5 = intersect_TP_SH5[,cbind(13:18,30:35,1:12)]
rownames(intersect_TP_SH5) = intersect_TP_SH5$X
#write.csv(intersect_TP_SH5, "./Result/Test/Result_Annotated/Result_TP_Annotated/intersect_TP_SH5.csv")


DEG_SH5 = read.csv("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_ko_vs_wt_p0.05_SH5_minus_spoVG_upp.csv", header = T)
rownames(DEG_SH5) = DEG_SH5$X
DPP_SH5 = read.csv("./DPP_SIG_minus_spoVG_upp/04062019_diffexp.pr.all_SH5_SIG_minus_spoVG_upp.csv", header = T)
rownames(DPP_SH5) = DPP_SH5$X
intersect_TP_SH5 = read.csv("./Result/intersect_TP_SH5.csv", header = T)
rownames(intersect_TP_SH5) = intersect_TP_SH5$a3
intersect_TP_SH5 = merge(intersect_TP_SH5, DEG_SH5, by.x="a3", by.y="X", all.x = TRUE, no.dups = TRUE)
intersect_TP_SH5 = merge(intersect_TP_SH5, DPP_SH5, by.x="a3", by.y="X", all.x = TRUE)
intersect_TP_SH5 = merge(intersect_TP_SH5, geneWizard, by.x="a3", by.y="locus", all.x = TRUE)
dim(intersect_TP_SH5)
intersect_TP_SH5 = intersect_TP_SH5[,cbind(1:24)]
rownames(intersect_TP_SH5) = intersect_TP_SH5$a3
dim(intersect_TP_SH5)
intersect_TP_SH5 = intersect_TP_SH5[,-1]
dim(intersect_TP_SH5)
#write.csv(intersect_TP_SH5, "./Result/Result_Annotated/Result_TP_Annotated/intersect_TP_SH5.csv")



DEG_WT_SH5_vs_SH2 = read.csv("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/diff.exp.gene/DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_WT_SH5_vs_SH2_p0.05_minus_spoVG_upp.csv", header = T)
rownames(DEG_WT_SH5_vs_SH2) = DEG_WT_SH5_vs_SH2$X
DPP_WT_SH5_vs_SH2 = read.csv("~/OneDrive - University of Warwick/WORK/Results/Proteomics/FINAL Result/Analysis - R-Script /Data/combined/output/diff.exp.SIG/diff.exp.SIG_minus_spoVG_upp/diffexp.pr_SIG_WT_SH5_vs_SH2_minus_spoVG_upp.csv", header = T)
rownames(DPP_WT_SH5_vs_SH2) = DPP_WT_SH5_vs_SH2$X
intersect_TP_WT_SH5_vs_SH2 = read.csv("./Result/Test/intersect_TP_WT.csv", header = T)
rownames(intersect_TP_WT_SH5_vs_SH2) = intersect_TP_WT_SH5_vs_SH2$a3
intersect_TP_WT_SH5_vs_SH2 = merge(intersect_TP_WT_SH5_vs_SH2, DEG_WT_SH5_vs_SH2, by.x="a3", by.y="X", all.x = TRUE, no.dups = TRUE)
intersect_TP_WT_SH5_vs_SH2 = merge(intersect_TP_WT_SH5_vs_SH2, DPP_WT_SH5_vs_SH2, by.x="a3", by.y="X", all.x = TRUE)
intersect_TP_WT_SH5_vs_SH2 = merge(intersect_TP_WT_SH5_vs_SH2, geneWizard, by.x="a3", by.y="locus", all.x = TRUE)
dim(intersect_TP_WT_SH5_vs_SH2)
intersect_TP_WT_SH5_vs_SH2 = intersect_TP_WT_SH5_vs_SH2[,cbind(1:24)]
rownames(intersect_TP_WT_SH5_vs_SH2) = intersect_TP_WT_SH5_vs_SH2$a3
dim(intersect_TP_WT_SH5_vs_SH2)
intersect_TP_WT_SH5_vs_SH2 = intersect_TP_WT_SH5_vs_SH2[,-1]
dim(intersect_TP_WT_SH5_vs_SH2)
#write.csv(intersect_TP_WT_SH5_vs_SH2, "./Result/Result_Annotated/Result_TP_Annotated/intersect_TP_WT_SH5_vs_SH2_Annot.csv")



#########################################################################################################
## Merge Cat2, Cat3 and Cat4 with TP intersect = 03/08/19
#########################################################################################################

setwd("~/OneDrive - University of Warwick/WORK/Results/Trans_Prot_COMBINED/Data/")

geneCategories = read.csv("./SubtiWiki Exports /geneCategories.csv", header = T)
# cat2_coping_with_stress = geneCategories[geneCategories$category2 == "Coping with stress",]
#write.csv(cat2_coping_with_stress, "../SubtiWiki Exports /split_categories/cat_all/cat2_coping_with_stress.csv", row.names = F)

cat2_coping_with_stress = read.csv("./SubtiWiki Exports /split_categories/cat_all/cat2_coping_with_stress.csv", header = T)
cat3_Biofilm_formation = read.csv("../SubtiWiki Exports /split_categories/cat3/cat3_Biofilm_formation.csv", header = T)

intersect_TP_M9 = read.csv("./Result/Result_Annotated/Result_TP_Annotated/intersect_TP_M9.csv", header = T)
intersect_TP_SH2 = read.csv("./Result/Result_Annotated/Result_TP_Annotated/intersect_TP_SH2.csv", header = T)
intersect_TP_SH5 = read.csv("./Result/Result_Annotated/Result_TP_Annotated/intersect_TP_SH5.csv", header = T)
intersect_TP_SH5_vs_SH2 = read.csv("./Result/Result_Annotated/Result_TP_Annotated/intersect_TP_WT_SH5_vs_SH2_Annot.csv", header = T)


intersect_TP_M9_cat2_coping_with_stress = merge(intersect_TP_M9, cat2_coping_with_stress, by.x="X", by.y="gene", all.x = TRUE)
dim(intersect_TP_M9_cat2_coping_with_stress)
intersect_TP_M9_cat2_coping_with_stress = intersect_TP_M9_cat2_coping_with_stress[,-2]
# write.csv(intersect_TP_M9_cat2_coping_with_stress, "./Result/GSEA_SGC/M9/intersect_TP_M9_cat2_coping_with_stress.csv", row.names = F)


#05/08/19

cat2_Sporulation = geneCategories[geneCategories$category2 == "Sporulation",]
# write.csv(cat2_Sporulation, "./SubtiWiki Exports /split_categories/cat2/cat2_Sporulation.csv", row.names = F)

intersect_TP_M9_cat2_Sporulation = merge(intersect_TP_M9, cat2_Sporulation, by.x="X", by.y="gene", all.x = TRUE)
dim(intersect_TP_M9_cat2_Sporulation)
intersect_TP_M9_cat2_Sporulation = intersect_TP_M9_cat2_Sporulation[,-2]
# write.csv(intersect_TP_M9_cat2_Sporulation, "./Result/GSEA_SGC/M9/intersect_TP_M9_cat2_Sporulation.csv", row.names = F)

#12/08/19

intersect_TP_SH2_cat2_Sporulation = merge(intersect_TP_SH2, cat2_Sporulation, by.x="X", by.y="gene", all.x = TRUE)
dim(intersect_TP_SH2_cat2_Sporulation)
intersect_TP_SH2_cat2_Sporulation = na.omit(intersect_TP_SH2_cat2_Sporulation)
dim(intersect_TP_SH2_cat2_Sporulation)
# write.csv(intersect_TP_SH2_cat2_Sporulation, "./Result/GSEA_SGC/SH2/intersect_TP_SH2_cat2_Sporulation.csv", row.names = F)

intersect_TP_SH5_cat2_Sporulation = merge(intersect_TP_SH5, cat2_Sporulation, by.x="X", by.y="gene", all.x = TRUE)
dim(intersect_TP_SH5_cat2_Sporulation)
intersect_TP_SH5_cat2_Sporulation = na.omit(intersect_TP_SH5_cat2_Sporulation)
dim(intersect_TP_SH5_cat2_Sporulation)
# write.csv(intersect_TP_SH5_cat2_Sporulation, "./Result/GSEA_SGC/SH5/intersect_TP_SH5_cat2_Sporulation.csv", row.names = F)

intersect_TP_SH5_vs_SH2_cat2_Sporulation = merge(intersect_TP_SH5_vs_SH2, cat2_Sporulation, by.x="X", by.y="gene", all.x = TRUE)
dim(intersect_TP_SH5_vs_SH2_cat2_Sporulation)
intersect_TP_SH5_vs_SH2_cat2_Sporulation = na.omit(intersect_TP_SH5_vs_SH2_cat2_Sporulation)
dim(intersect_TP_SH5_vs_SH2_cat2_Sporulation)
# write.csv(intersect_TP_SH5_vs_SH2_cat2_Sporulation, "./Result/GSEA_SGC/SH5_vs_SH2/intersect_TP_SH5_vs_SH2_cat2_Sporulation.csv", row.names = F)

#06/08/19

cat3_Translation = geneCategories[geneCategories$category3 == "Translation",]
# write.csv(cat3_Translation, "./SubtiWiki Exports /split_categories/cat3/cat3_Translation.csv", row.names = F)
intersect_TP_SH2_cat3_Translation = merge(intersect_TP_SH2, cat3_Translation, by.x="X", by.y="gene", all.x = TRUE)
dim(intersect_TP_SH2_cat3_Translation)
intersect_TP_SH2_cat3_Translation = intersect_TP_SH2_cat3_Translation[,-2]
# write.csv(intersect_TP_SH2_cat3_Translation, "./Result/GSEA_SGC/SH2/intersect_TP_SH2_cat3_Translation.csv", row.names = F)

#07/08/19

geneRegulations = read.csv("./SubtiWiki Exports /regulations.csv", header = T)

intersect_TP_M9_geneRegulations = merge(intersect_TP_M9, geneRegulations, by.x="X", by.y="locus.tag", all.x = TRUE)
dim(intersect_TP_M9_geneRegulations)
# write.csv(intersect_TP_M9_geneRegulations, "./Result/geneRegulations/intersect_TP_M9_geneRegulations.csv", row.names = F, na = "")

intersect_TP_SH2_geneRegulations = merge(intersect_TP_SH2, geneRegulations, by.x="X", by.y="locus.tag", all.x = TRUE)
dim(intersect_TP_SH2_geneRegulations)
# write.csv(intersect_TP_SH2_geneRegulations, "./Result/geneRegulations/intersect_TP_SH2_geneRegulations.csv", row.names = F, na = "")

intersect_TP_SH5_geneRegulations = merge(intersect_TP_SH5, geneRegulations, by.x="X", by.y="locus.tag", all.x = TRUE)
dim(intersect_TP_SH5_geneRegulations)
# write.csv(intersect_TP_SH5_geneRegulations, "./Result/geneRegulations/intersect_TP_SH5_geneRegulations.csv", row.names = F, na = "")

intersect_TP_SH5_vs_SH2_geneRegulations = merge(intersect_TP_SH5_vs_SH2, geneRegulations, by.x="X", by.y="locus.tag", all.x = TRUE)
dim(intersect_TP_SH5_vs_SH2_geneRegulations)
# write.csv(intersect_TP_SH5_vs_SH2_geneRegulations, "./Result/geneRegulations/intersect_TP_SH5_vs_SH2_geneRegulations.csv", row.names = F, na = "")






































