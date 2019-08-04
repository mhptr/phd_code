#This scripts adds different annotations by merging .csv files


setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/diff.exp.gene/")

geneWizard = read.csv("../SubtiWiki Exports /Genes export wizard.csv", header = T)
geneCategories = read.csv("../SubtiWiki Exports /geneCategories.csv", header=T)
geneRegulators = read.csv("../SubtiWiki Exports /regulations.csv", header=T)

################################
## ELD - LB - ko vs wt 
################################

DEG_ko_LB = read.csv("./DEG_ED/DEG_ko_vs_wt_p0.05_LB_eld.csv", header = T)
DEG_ko_LB_Annot = merge(DEG_ko_LB, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_ko_LB_Annot) = DEG_ko_LB_Annot$X
DEG_ko_LB_Annot = DEG_ko_LB_Annot[,-c(1,18,19)]
#write.csv(DEG_ko_LB_Annot, "./DEG_ED/DEG_ko_LB_Annot_eld.csv")

#DEG_ko_LB_Annot = read.csv("./DEG_ED/DEG_ko_LB_Annot_eld.csv", header = T)
#DEG_ko_LB_Annot_Reg = merge(DEG_ko_LB_Annot, geneRegulators, by.x="X", by.y="locus.tag")
##DEG_ko_LB_Annot_Reg = merge(DEG_ko_LB_Annot, geneRegulators, by.x="X", by.y="regulator.locus.tag")
#rownames(DEG_ko_LB_Annot_Reg) = DEG_ko_LB_Annot_Reg$X
#DEG_ko_LB_Annot_Reg = DEG_ko_LB_Annot_Reg[,-1]
#DEG_ko_LB_Annot_Reg = DEG_ko_LB_Annot_Reg[,-c(1,18,19)]
#write.csv(DEG_ko_LB_Annot_Reg, "./DEG_ED/DEG_ko_LB_Annot_Reg_eld.csv")



################################
## ELD - LB - dc vs wt 
################################

DEG_dc_LB = read.csv("./DEG_ED/resSig.LB.05_dc_surC.csv", header = T)
DEG_dc_LB_Annot = merge(DEG_dc_LB, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_dc_LB_Annot) = DEG_dc_LB_Annot$X
DEG_dc_LB_Annot = DEG_dc_LB_Annot[,-c(1,19,20)]
#write.csv(DEG_dc_LB_Annot, "./DEG_ED/DEG_dc_LB_Annot_eld.csv")



################################
## ELD - LB - ko vs dc 
################################

DEG_ko_vs_dc_LB = read.csv("./DEG_ED/LB/DEG_raw_LB_eld/resSig.LB.05_ko_vs_dc.csv", header = T)
DEG_ko_vs_dc_LB_Annot = merge(DEG_ko_vs_dc_LB, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_ko_vs_dc_LB_Annot) = DEG_ko_vs_dc_LB_Annot$X
DEG_ko_vs_dc_LB_Annot = DEG_ko_vs_dc_LB_Annot[,-c(1,19,20)]
#write.csv(DEG_ko_vs_dc_LB_Annot, "./DEG_ED/LB/DEG_ko_vs_dc_LB_Annot_eld.csv")



################################
## ELD - M9 - ko vs wt 
################################

DEG_ko_M9 = read.csv("./DEG_ED/resSig.M9.05_ko.csv", header = T)
DEG_ko_M9_Annot = merge(DEG_ko_M9, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_ko_M9_Annot) = DEG_ko_M9_Annot$X
DEG_ko_M9_Annot = DEG_ko_M9_Annot[,-c(1,19,20)]
#write.csv(DEG_ko_M9_Annot, "./DEG_ED/DEG_ko_M9_Annot_eld.csv")


################################
## ELD - M9 - dc vs wt 
################################

DEG_dc_M9 = read.csv("./DEG_ED/resSig.M9.05_dc.csv", header = T)
DEG_dc_M9_Annot = merge(DEG_dc_M9, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_dc_M9_Annot) = DEG_dc_M9_Annot$X
DEG_dc_M9_Annot = DEG_dc_M9_Annot[,-c(1,19,20)]
#write.csv(DEG_dc_M9_Annot, "./DEG_ED/DEG_dc_M9_Annot_eld.csv")


################################
## ELD - M9 - ko vs dc 
################################

DEG_ko_vs_dc_M9 = read.csv("./DEG_ED/M9/DEG_raw_M9_eld/resSig.M9.05_ko_vs_dc.csv", header = T)
DEG_ko_vs_dc_M9_Annot = merge(DEG_ko_vs_dc_M9, geneWizard, by.x="X", by.y="locus", all.x = TRUE)
rownames(DEG_ko_vs_dc_M9_Annot) = DEG_ko_vs_dc_M9_Annot$X
DEG_ko_vs_dc_M9_Annot = DEG_ko_vs_dc_M9_Annot[,-c(1,19,20)]
#write.csv(DEG_ko_vs_dc_M9_Annot, "./DEG_ED/M9/DEG_ko_vs_dc_M9_Annot_eld.csv")


