
#29/07/19
setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/diff.exp.gene/DEG_SM/")


DEG_KEGG = read.csv("./DEG_KEGG/kegg.enrich.sig_transcriptomics.csv", header = T)
rownames(DEG_KEGG) = DEG_KEGG$X
DEG_KEGG = DEG_KEGG[,-1]

# strsplit(x, split)
# unlist(strsplit("a.b.c", "."))


#under contruction 
my.de.path.gene = list(DEG_KEGG$M9_my.de.path.gene, DEG_KEGG$SH2_my.de.path.gene, DEG_KEGG$SH5_my.de.path.gene)
gene_split = c()
for(i in 1:length(my.de.path.gene)){
  my.de = my.de.path.gene[[i]]
  list_gene = unlist(strsplit(as.character(my.de.path.gene), "\\|"))
  gene_split = cbind()
}
  

#works
# M9_my.de.path.gene = unlist(strsplit(as.character(DEG_KEGG$M9_my.de.path.gene), "\\|"))
# write.csv(M9_my.de.path.gene, "./DEG_KEGG/M9_my.de.path.gene.csv")
# 
















