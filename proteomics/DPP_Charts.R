#making graphs and charts of DPP
#18/08/19

setwd("~/OneDrive - University of Warwick/WORK/Results/Proteomics/FINAL Result/Analysis - R-Script /Data/")
# 
# DPP_LB_regulations = read.csv("./combined/output/geneRegulations/old/DPP_LB_geneRegulations.csv", header = T)
# DPP_M9_regulations = read.csv("./combined/output/geneRegulations/old/DPP_M9_geneRegulations.csv", header = T)
# DPP_SH2_regulations = read.csv("./combined/output/geneRegulations/old/DPP_SH2_geneRegulations.csv", header = T)
# DPP_SH5_regulations = read.csv("./combined/output/geneRegulations/old/DPP_SH5_geneRegulations.csv", header = T)
# DPP_SH5_vs_SH2_regulations = read.csv("./combined/output/geneRegulations/old/DPP_SH5_vs_SH2_geneRegulations.csv", header = T)
# 
# dim(DPP_LB_regulations)
# DPP_LB_regulations = DPP_LB_regulations[!DPP_LB_regulations$regulator == "",] 
# dim(DPP_LB_regulations)
# write.csv(DPP_LB_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_LB_regulations.csv", row.names = F)
# 
# dim(DPP_M9_regulations)
# DPP_M9_regulations = DPP_M9_regulations[!DPP_M9_regulations$regulator == "",] 
# dim(DPP_M9_regulations)
# write.csv(DPP_M9_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_M9_regulations.csv", row.names = F)
# 
# dim(DPP_SH2_regulations)
# DPP_SH2_regulations = DPP_SH2_regulations[!DPP_SH2_regulations$regulator == "",] 
# dim(DPP_SH2_regulations)
# write.csv(DPP_SH2_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_SH2_regulations.csv", row.names = F)
# 
# dim(DPP_SH5_regulations)
# DPP_SH5_regulations = DPP_SH5_regulations[!DPP_SH5_regulations$regulator == "",] 
# dim(DPP_SH5_regulations)
# write.csv(DPP_SH5_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_SH5_regulations.csv", row.names = F)
# 
# dim(DPP_SH5_vs_SH2_regulations)
# DPP_SH5_vs_SH2_regulations = DPP_SH5_vs_SH2_regulations[!DPP_SH5_vs_SH2_regulations$regulator == "",] 
# dim(DPP_SH5_vs_SH2_regulations)
# write.csv(DPP_SH5_vs_SH2_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_SH5_cs_SH2_regulations.csv", row.names = F)
# 
# 
# DPP_LB_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_LB_regulations.csv", header = T)
# DPP_M9_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_M9_regulations.csv", header = T)
# DPP_SH2_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_SH2_regulations.csv", header = T)
# DPP_SH5_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_SH5_regulations.csv", header = T)
# DPP_SH5_vs_SH2_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_SH5_cs_SH2_regulations.csv", header = T)
# 
# # DPP_all = list(DPP_LB_regulations, DPP_M9_regulations, DPP_SH2_regulations, DPP_SH5_regulations, DPP_SH5_vs_SH2_regulations)
# # drop_empty_row = for (i in DPP_all) {
# #   DPP_regulations = DPP_all(i)[!DPP_all(i)$regulator == "",] 
# # }
# 
# 
# 
# DPP_LB_regulators = table(DPP_LB_regulations$regulator)
# DPP_M9_regulators = table(DPP_M9_regulations$regulator)
# DPP_SH2_regulators = table(DPP_SH2_regulations$regulator)
# DPP_SH5_regulators = table(DPP_SH5_regulations$regulator)
# DPP_SH5_vs_SH2_regulators = table(DPP_SH5_vs_SH2_regulations$regulator)
# 
# write.csv(DPP_LB_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_LB_regulators.csv", row.names = F)
# DPP_LB_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_LB_regulators.csv", header = T)
# colnames(DPP_LB_regulators) = c("regulator", "frequency")
# dim(DPP_LB_regulators)
# DPP_LB_regulators = DPP_LB_regulators[!DPP_LB_regulators$regulator == "",]
# dim(DPP_LB_regulators)
# View(DPP_LB_regulators)
# write.csv(DPP_LB_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_LB_regulators.csv", row.names = F)
# 
# # DPP_M9_regulators = DPP_M9_regulators[!DPP_M9_regulators$Name == "",] #removes the empty row
# write.csv(DPP_M9_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", row.names = F)
# DPP_M9_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", header = T)
# colnames(DPP_M9_regulators) = c("regulator", "frequency")
# dim(DPP_M9_regulators)
# DPP_M9_regulators = DPP_M9_regulators[!DPP_M9_regulators$regulator == "",]
# dim(DPP_M9_regulators)
# View(DPP_M9_regulators)
# write.csv(DPP_M9_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", row.names = F)
# 
# write.csv(DPP_SH2_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH2_regulators.csv", row.names = F)
# DPP_SH2_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH2_regulators.csv", header = T)
# colnames(DPP_SH2_regulators) = c("regulator", "frequency")
# dim(DPP_SH2_regulators)
# DPP_SH2_regulators = DPP_SH2_regulators[!DPP_SH2_regulators$regulator == "",]
# dim(DPP_SH2_regulators)
# View(DPP_SH2_regulators)
# write.csv(DPP_SH2_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH2_regulators.csv", row.names = F)
# 
# write.csv(DPP_SH5_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_regulators.csv", row.names = F)
# DPP_SH5_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_regulators.csv", header = T)
# colnames(DPP_SH5_regulators) = c("regulator", "frequency")
# dim(DPP_SH5_regulators)
# DPP_SH5_regulators = DPP_SH5_regulators[!DPP_SH5_regulators$regulator == "",]
# dim(DPP_SH5_regulators)
# View(DPP_SH5_regulators)
# write.csv(DPP_SH5_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_regulators.csv", row.names = F)
# 
# write.csv(DPP_SH5_vs_SH2_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_vs_SH2_regulators.csv", row.names = F)
# DPP_SH5_vs_SH2_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_vs_SH2_regulators.csv", header = T)
# colnames(DPP_SH5_vs_SH2_regulators) = c("regulator", "frequency")
# dim(DPP_SH5_vs_SH2_regulators)
# DPP_SH5_vs_SH2_regulators = DPP_SH5_vs_SH2_regulators[!DPP_SH5_vs_SH2_regulators$regulator == "",]
# dim(DPP_SH5_vs_SH2_regulators)
# View(DPP_SH5_vs_SH2_regulators)
# write.csv(DPP_SH5_vs_SH2_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_vs_SH2_regulators.csv", row.names = F)


DPP_LB_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_LB_regulators.csv", header = T)
DPP_M9_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", header = T)
DPP_SH2_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH2_regulators.csv", header = T)
DPP_SH5_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_regulators.csv", header = T)
DPP_SH5_vs_SH2_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_vs_SH2_regulators.csv", header = T)
























