#adding annotation to counts table 


setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/diff.exp.gene/")


geneWizard=read.csv("../SubtiWiki Exports /Genes export wizard.csv",header=T)
geneCategories=read.csv("../SubtiWiki Exports /geneCategories.csv",header=T)
geneRegulators=read.csv("../SubtiWiki Exports /regulations.csv",header=T)



# counts_SM = read.csv("../filt_txt_all_2/counts.all_SM_S1583.csv",header=T)
# counts_SM_annot = merge(counts_SM, geneWizard,by.x="X",by.y="locus",all.x=TRUE)
# rownames(counts_SM_annot)=counts_SM_annot$X
# counts_SM_annot=counts_SM_annot[,-c(21:32)]
# write.csv(counts_SM_annot, "../filt_txt_all_2/counts.all_annot_SM.csv", row.names = F)


