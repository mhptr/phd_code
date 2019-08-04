
setwd("~/OneDrive - University of Warwick/WORK/Results/Trans_Prot_COMBINED/Data/")

library(VennDiagram)
library(RColorBrewer)
swati.color = c(brewer.pal(8,"Dark2"))

#13/07/2019

venn_intersect = function(DEG,DPP, tag) {
  
  DegFileLocation = paste("./DEG_SIG_minus_spoVG_upp/", DEG, sep = "")
  DppFileLocation = paste("./DPP_SIG_minus_spoVG_upp/", DPP, sep = "")
  DEG_File = read.csv(DegFileLocation, header = T)
  DPP_File = read.csv(DppFileLocation, header = T)
  rownames(DEG_File) = DEG_File$X
  DEG_File = DEG_File[,-1]
  rownames(DPP_File) = DPP_File$X
  DPP_File = DPP_File[,-1]
  
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

  DEG_File_rowname = rownames(DEG_File)
  DPP_File_rowname = rownames(DPP_File)
  conditions = list(DEG_File_rowname, DPP_File_rowname)
  common = calculate.overlap(conditions)
  #intersect DEG and DPP
  intersect_TP = common["a3"]
  #write.csv(intersect_TP, paste("./Result/Test/intersect_TP_", tag, ".csv", sep = ""),row.names=FALSE)
  DEG_minus_intersect_TP = common["a1"]
  #write.csv(DEG_minus_intersect_TP, paste("./Result/Test/DEG_minus_intersect_TP_", tag, ".csv", sep = ""),row.names=FALSE)
  DPP_minus_intersect_TP = common["a2"]
 #write.csv(DPP_minus_intersect_TP, paste("./Result/Test/DPP_minus_intersect_TP_", tag, ".csv", sep = ""),row.names=FALSE)
}

# venn_intersect("DEG_ko_M9_Annot_minus_spoVG_upp.csv", "04062019_diffexp.pr.all_M9_SIG_minus_spoVG_upp.csv", "M9")
# venn_intersect("DEG_ko_SH2_Annot_minus_spoVG_upp.csv", "04062019_diffexp.pr.all_SH2_SIG_minus_spoVG_upp.csv", "SH2")
# venn_intersect("DEG_ko_SH5_Annot_minus_spoVG_upp.csv", "04062019_diffexp.pr.all_SH5_SIG_minus_spoVG_upp.csv", "SH5")
#define the correct numbers in intersection

venn_intersect("DEG_WT_SH5_vs_SH2_Annot.csv", "Diff.pr_SIG_SH5_vs_SH2_Annotated.csv", "WT")  #WT_SH5_vs_SH2



