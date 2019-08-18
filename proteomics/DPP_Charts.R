setwd("~/OneDrive - University of Warwick/WORK/Results/Proteomics/FINAL Result/Analysis - R-Script /Data/")


DPP_LB_regulations = read.csv("./combined/output/geneRegulations/old/DPP_LB_geneRegulations.csv", header = T)
DPP_M9_regulations = read.csv("./combined/output/geneRegulations/old/DPP_M9_geneRegulations.csv", header = T)
DPP_SH2_regulations = read.csv("./combined/output/geneRegulations/old/DPP_SH2_geneRegulations.csv", header = T)
DPP_SH5_regulations = read.csv("./combined/output/geneRegulations/old/DPP_SH5_geneRegulations.csv", header = T)
DPP_SH5_vs_SH2_regulations = read.csv("./combined/output/geneRegulations/old/DPP_SH5_vs_SH2_geneRegulations.csv", header = T)


length(DPP_M9_regulations$regulator)


DPP_M9_regulators = table(DPP_M9_regulations$regulator)
# DPP_M9_regulators = DPP_M9_regulators[!DPP_M9_regulators$Name == "",] #removes the empty row
write.csv(DPP_M9_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", row.names = F)
View(DPP_M9_regulators)
