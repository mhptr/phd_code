#Thisscriptsaddsdifferentannotationsbymerging.csvfiles


setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/diff.exp.gene/")


geneWizard=read.csv("../SubtiWiki Exports /Genes export wizard.csv",header=T)
geneCategories=read.csv("../SubtiWiki Exports /geneCategories.csv",header=T)
geneRegulators=read.csv("../SubtiWiki Exports /regulations.csv",header=T)

################################
## WT_SH5_vs_SH2
################################

DEG_WT_SH5_vs_SH2=read.csv("./DEG_SM/Diff_SIG_sm/Diff_SIG_minus_spoVG_upp/DEG_WT_SH5_vs_SH2_p0.05_minus_spoVG_upp.csv",header=T)
DEG_WT_SH5_vs_SH2_Annot=merge(DEG_WT_SH5_vs_SH2,geneWizard,by.x="X",by.y="locus",all.x=TRUE)
rownames(DEG_WT_SH5_vs_SH2_Annot)=DEG_WT_SH5_vs_SH2_Annot$X
DEG_WT_SH5_vs_SH2_Annot=DEG_WT_SH5_vs_SH2_Annot[,-c(1, 19, 20)]
#write.csv(DEG_WT_SH5_vs_SH2_Annot,"./DEG_SM/Annotated_sm/Annot_diff.exp_minus_spoVG_upp/Data/DEG_WT_SH5_vs_SH2_Annot.csv")













