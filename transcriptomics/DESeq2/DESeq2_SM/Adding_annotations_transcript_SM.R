#Thisscriptsaddsdifferentannotationsbymerging.csvfiles


setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/diff.exp.gene/")


geneWizard=read.csv("../SubtiWiki Exports /Genes export wizard.csv",header=T)
geneCategories=read.csv("../SubtiWiki Exports /geneCategories.csv",header=T)
geneRegulators=read.csv("../SubtiWiki Exports /regulations.csv",header=T)

################################
##SM-M9-kovswt
################################

DEG_ko_M9=read.csv("./DEG_SM/Diff_SIG_sm/DEG_ko_vs_wt_p0.05_M9_sm.csv",header=T)
DEG_ko_M9_Annot=merge(DEG_ko_M9,geneWizard,by.x="X",by.y="locus",all.x=TRUE)
rownames(DEG_ko_M9_Annot)=DEG_ko_M9_Annot$X
DEG_ko_M9_Annot=DEG_ko_M9_Annot[,-1]
#write.csv(DEG_ko_M9_Annot,"./DEG_SM/Annotated_sm/DEG_ko_M9_Annot_sm_1.csv")

DEG_ko_M9_Annot=read.csv("./DEG_SM/Annotated_sm/DEG_ko_M9_Annot_sm.csv",header=T)
DEG_ko_M9_Annot_Reg=merge(DEG_ko_M9_Annot,geneRegulators,by.x="X",by.y="locus.tag")
#DEG_ko_M9_Annot_Reg=merge(DEG_ko_M9_Annot,geneRegulators,by.x="X",by.y="regulator.locus.tag")
rownames(DEG_ko_M9_Annot_Reg)=DEG_ko_M9_Annot_Reg$X
DEG_ko_M9_Annot_Reg=DEG_ko_M9_Annot_Reg[,-1]
#write.csv(DEG_ko_M9_Annot_Reg,"./DEG_SM/Annotated_sm/DEG_ko_M9_Annot_Reg_sm.csv")


################################
##SM-SH2-kovswt
################################


DEG_ko_SH2=read.csv("./DEG_SM/DEG_ko_vs_wt_p0.05_SH2_sm.csv",header=T)
DEG_ko_SH2_Annot=merge(DEG_ko_SH2,geneWizard,by.x="X",by.y="locus",all.x=TRUE)
rownames(DEG_ko_SH2_Annot)=DEG_ko_SH2_Annot$X
DEG_ko_SH2_Annot=DEG_ko_SH2_Annot[,-1]
#write.csv(DEG_ko_SH2_Annot,"./DEG_SM/Annotated_sm/DEG_ko_SH2_Annot_sm.csv")

DEG_ko_SH2_Annot=read.csv("./DEG_SM/Annotated_sm/DEG_ko_SH2_Annot_sm.csv",header=T)
DEG_ko_SH2_Annot_Reg=merge(DEG_ko_SH2_Annot,geneRegulators,by.x="X",by.y="locus.tag")
#DEG_ko_SH2_Annot_Reg=merge(DEG_ko_SH2_Annot,geneRegulators,by.x="X",by.y="regulator.locus.tag")
rownames(DEG_ko_SH2_Annot_Reg)=DEG_ko_SH2_Annot_Reg$X
DEG_ko_SH2_Annot_Reg=DEG_ko_SH2_Annot_Reg[,-1]
#write.csv(DEG_ko_SH2_Annot_Reg,"./DEG_SM/DEG_ko_SH2_Annot_Reg_sm.csv")


################################
##SM-SH5-kovswt
################################


DEG_ko_SH5=read.csv("./DEG_SM/Diff_SIG_sm/DEG_ko_vs_wt_p0.05_SH5_sm.csv",header=T)
DEG_ko_SH5_Annot=merge(DEG_ko_SH5,geneWizard,by.x="X",by.y="locus",all.x=TRUE)
rownames(DEG_ko_SH5_Annot)=DEG_ko_SH5_Annot$X
DEG_ko_SH5_Annot=DEG_ko_SH5_Annot[,-1]
#write.csv(DEG_ko_SH5_Annot,"./DEG_SM/Annotated_sm/DEG_ko_SH5_Annot_sm_30062019.csv")

DEG_ko_SH5_Annot=read.csv("./DEG_SM/Annotated_sm/DEG_ko_SH5_Annot_sm.csv",header=T)
DEG_ko_SH5_Annot_Reg=merge(DEG_ko_SH5_Annot,geneRegulators,by.x="X",by.y="locus.tag")
#DEG_ko_SH5_Annot_Reg=merge(DEG_ko_SH5_Annot,geneRegulators,by.x="X",by.y="regulator.locus.tag")
rownames(DEG_ko_SH5_Annot_Reg)=DEG_ko_SH5_Annot_Reg$X
DEG_ko_SH5_Annot_Reg=DEG_ko_SH5_Annot_Reg[,-1]
#write.csv(DEG_ko_SH5_Annot_Reg,"./DEG_SM/Annotated_sm/DEG_ko_SH5_Annot_Reg_sm.csv")



################################################################
##
################################################################



Heat_shock_proteins=read.csv("./combined/output/GSEA_SGC/Heat_shock_proteins.csv",header=T)
Heat_shock_proteins=Heat_shock_proteins[,-1]

diffexp.pr.all.M9_SIG_Annotated=read.csv("./combined/output/diffexp.pr.all.M9_SIG_Annotated.csv",header=T)
Heat_shock_proteins_M9=merge(diffexp.pr.all.M9_SIG_Annotated,Heat_shock_proteins,by.x="X",by.y="gene",all.y=TRUE)
rownames(Heat_shock_proteins_M9)=Heat_shock_proteins_M9$X
colnames(Heat_shock_proteins_M9)=colnames(Heat_shock_proteins_M9)

Heat_shock_proteins_M9=Heat_shock_proteins_M9[,cbind("name","description",
"product","essential","names",
"logFC","adj.P.Val")]

#write.csv(Heat_shock_proteins_M9,"./combined/output/GSEA_SGC/Heat_shock_proteins_M9.csv")



#################
###SH2
#################

Translation=read.csv("./combined/output/GSEA_SGC/Translation.csv",header=T)
Translation=Translation[,-1]

diffexp.pr.all.SH2_SIG_Annotated=read.csv("./combined/output/diffexp.pr.all.SH2_SIG_Annotated.csv",header=T)
Translation_SH2=merge(diffexp.pr.all.SH2_SIG_Annotated,Translation,by.x="X",by.y="gene",all.y=TRUE)
rownames(Translation_SH2)=Translation_SH2$X
colnames(Translation_SH2)=colnames(Translation_SH2)

Translation_SH2=Translation_SH2[,cbind("name","description",
"product","essential","names",
"logFC","adj.P.Val")]

#write.csv(Translation_SH2,"./combined/output/GSEA_SGC/Translation_SH2.csv")




#################
###SH5
#################

Sporulation_proteins=read.csv("./combined/output/GSEA_SGC/Sporulation_proteins.csv",header=T)
Sporulation_proteins=Sporulation_proteins[,-1]

diffexp.pr.all.SH5_SIG_Annotated=read.csv("./combined/output/diffexp.pr.all.SH5_SIG_Annotated.csv",header=T)
Sporulation_proteins_SH5=merge(diffexp.pr.all.SH5_SIG_Annotated,Sporulation_proteins,by.x="X",by.y="gene",all.y=TRUE)
rownames(Sporulation_proteins_SH5)=Sporulation_proteins_SH5$X
colnames(Sporulation_proteins_SH5)=colnames(Sporulation_proteins_SH5)

Sporulation_proteins_SH5=Sporulation_proteins_SH5[,cbind("name","description",
"product","essential","names",
"logFC","adj.P.Val")]

#write.csv(Sporulation_proteins_SH5,"./combined/output/GSEA_SGC/Sporulation_proteins_SH5.csv")
#Sporulation_proteins_SH5=Sporulation_proteins_SH5[!is.na(Sporulation_proteins_SH5)>2,]




################################################################################
################################################################################
################################################################################



#Addingannotationtocommontranscripts


setwd("~/OneDrive-UniversityofWarwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/")

Common_transcripts=read.csv("./diff.exp.gene/DEG_SM/Annotated_sm/Annot_diff.exp_minus_spoVG_upp/Common_transcripts.csv",header=T)
geneWizard=read.csv("./SubtiWikiExports/Genesexportwizard.csv",header=T)


Common_transcripts_Annot=merge(Common_transcripts,geneWizard,by.x="X",by.y="locus",all.x=TRUE)
rownames(Common_transcripts_Annot)=Common_transcripts_Annot$X
Common_transcripts_Annot=Common_transcripts_Annot[,-1]
#write.csv(Common_transcripts_Annot,"./diff.exp.gene/DEG_SM/Annotated_sm/Annot_diff.exp_minus_spoVG_upp/Common_transcripts_Annot.csv")


#mergeLFCandpadj.


Common_transcripts_Annot=read.csv("./diff.exp.gene/DEG_SM/Annotated_sm/Annot_diff.exp_minus_spoVG_upp/Common_transcripts_Annot.csv",header=T)
M9_SIG_p0.05=read.csv("./diff.exp.gene/DEG_SM/Annotated_sm/Annot_diff.exp_minus_spoVG_upp/DEG_ko_M9_Annot_minus_spoVG_upp.csv",header=T)
SH2_SIG_p0.05=read.csv("./diff.exp.gene/DEG_SM/Annotated_sm/Annot_diff.exp_minus_spoVG_upp/DEG_ko_SH2_Annot_minus_spoVG_upp.csv",header=T)
SH5_SIG_p0.05=read.csv("./diff.exp.gene/DEG_SM/Annotated_sm/Annot_diff.exp_minus_spoVG_upp/DEG_ko_SH5_Annot_minus_spoVG_upp.csv",header=T)
Common_M9_SIG=merge(Common_transcripts_Annot,M9_SIG_p0.05,by.x="X",by.y="X",all.x=TRUE)
Common_M9_SIG=Common_M9_SIG[,cbind("X","name.x","description.x","mw.x","proteinLength.x","geneLength.x","function..x",
"product.x","essential.x","names.x","pI.x","ec.x","log2FoldChange","padj")]
Common_SH2_SIG=merge(Common_M9_SIG,SH2_SIG_p0.05,by.x="X",by.y="X",all.x=TRUE)
Common_SH2_SIG = Common_SH2_SIG[, cbind("X", "name.x", "description.x", "mw.x", "proteinLength.x", "geneLength.x", "function..x", 
                                        "product.x", "essential.x", "names.x", "pI.x", "ec.x", "log2FoldChange.x", "padj.x",
                                        "log2FoldChange.y", "padj.y")]
colnames(Common_SH2_SIG) = c("X", "name", "description", "mw", "proteinLength", "geneLength", "function", 
                             "product", "essential", "names", "pI", "ec", "log2FoldChange_M9", "padj_M9",
                             "log2FoldChange_SH2", "padj_SH2")

Common_SH5_SIG = merge(Common_SH2_SIG,SH5_SIG_p0.05, by.x="X",by.y="X",all.x=TRUE)
Common_SH5_SIG = Common_SH5_SIG[, cbind("X", "name.x", "description.x", "mw.x", "proteinLength.x", "geneLength.x", 
                                        "function", "product.x", "essential.x", "names.x", "pI.x", "ec.x",
                                        "log2FoldChange_M9", "padj_M9", "log2FoldChange_SH2", "padj_SH2", 
                                        "log2FoldChange", "padj")]
colnames(Common_SH5_SIG) = c("X", "name", "description", "mw", "proteinLength", "geneLength", "function", 
                             "product", "essential", "names", "pI", "ec", "log2FoldChange_M9", "padj_M9",
                             "log2FoldChange_SH2", "padj_SH2", "log2FoldChange_SH5", "padj_SH5")
rownames(Common_SH5_SIG) = Common_SH5_SIG$X
Common_SH5_SIG = Common_SH5_SIG[,-1]
View(Common_SH5_SIG)
#write.csv(Common_SH5_SIG, "./diff.exp.gene/DEG_SM/Annotated_sm/Annot_diff.exp_minus_spoVG_upp/Common_Conditions_SIG_transcriptomics_minus_spoVG_upp_Annot.csv")





















