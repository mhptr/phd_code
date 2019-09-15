#06/09/19

setwd("~/OneDrive - University of Warwick/WORK/Results/Proteomics/FINAL Result/Analysis - R-Script /")

# regulation = read.csv("./SubtiWiki Exports /regulations.csv", header = T)
# regulation_freq = table(regulation$regulator)
# write.csv(regulation_freq, "./SubtiWiki Exports /regulation_freq.csv", row.names = F)

regulation_freq = read.csv("./SubtiWiki Exports /regulation_freq.csv", header = T)
colnames(regulation_freq) = c("BSU_regulators_all", "No_of_total_targets")

find_ratio = function(regulation_freq, my_table_freq) {
  merged_table = merge(my_table_freq, regulation_freq, by.x="regulator", by.y="BSU_regulators_all", all.x = TRUE)
  print(merged_table)
  merged_percentage_table = transform(merged_table, percentage = 100 * merged_table$No_of_targets/merged_table$No_of_total_targets)
  merged_percentage_table = merged_percentage_table[order(merged_percentage_table$percentage,  decreasing = TRUE), ]
  return(merged_percentage_table)
}


#DPP_LB
DPP_LB_geneRegulations = read.csv("./Data/combined/output/geneRegulations/old/regulations/DPP_LB_geneRegulations.csv", header = T)
DPP_Up_LB_geneRegulations = DPP_LB_geneRegulations[(DPP_LB_geneRegulations$logFC)>0,]
DPP_Up_LB_geneRegulations = table(DPP_Up_LB_geneRegulations$regulator)
write.csv(DPP_Up_LB_geneRegulations, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_LB_reg_freq.csv", row.names = F)
DPP_Up_LB_reg_freq = read.csv("./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_LB_reg_freq.csv", header = T)
colnames(DPP_Up_LB_reg_freq) = c("regulator", "No_of_targets")
DPP_Up_LB_reg_freq = DPP_Up_LB_reg_freq[(DPP_Up_LB_reg_freq$No_of_targets)>0,]
DPP_Up_LB_reg_freq = DPP_Up_LB_reg_freq[!DPP_Up_LB_reg_freq$regulator == "",]
DPP_Up_LB_reg_percentage = find_ratio(regulation_freq, DPP_Up_LB_reg_freq)
write.csv(DPP_Up_LB_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/Up/DPP_LB_reg_percentage.csv", row.names = FALSE)
View(DPP_Up_LB_reg_percentage)

#DPP_M9
DPP_M9_geneRegulations = read.csv("./Data/combined/output/geneRegulations/old/regulations/DPP_M9_geneRegulations.csv", header = T)
DPP_Up_M9_geneRegulations = DPP_M9_geneRegulations[(DPP_M9_geneRegulations$logFC)>0,]
DPP_Up_M9_geneRegulations = table(DPP_Up_M9_geneRegulations$regulator)
write.csv(DPP_Up_M9_geneRegulations, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_M9_reg_freq.csv", row.names = F)
DPP_Up_M9_reg_freq = read.csv("./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_M9_reg_freq.csv", header = T)
colnames(DPP_Up_M9_reg_freq) = c("regulator", "No_of_targets")
DPP_Up_M9_reg_freq = DPP_Up_M9_reg_freq[(DPP_Up_M9_reg_freq$No_of_targets)>0,]
DPP_Up_M9_reg_freq = DPP_Up_M9_reg_freq[!DPP_Up_M9_reg_freq$regulator == "",]
DPP_Up_M9_reg_percentage = find_ratio(regulation_freq, DPP_Up_M9_reg_freq)
write.csv(DPP_Up_M9_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/Up/DPP_M9_reg_percentage.csv", row.names = FALSE)
View(DPP_Up_M9_reg_percentage)

#DPP_SH2
DPP_SH2_geneRegulations = read.csv("./Data/combined/output/geneRegulations/old/regulations/DPP_SH2_geneRegulations.csv", header = T)
DPP_Up_SH2_geneRegulations = DPP_SH2_geneRegulations[(DPP_SH2_geneRegulations$logFC)>0,]
DPP_Up_SH2_geneRegulations = table(DPP_Up_SH2_geneRegulations$regulator)
write.csv(DPP_Up_SH2_geneRegulations, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_SH2_reg_freq.csv", row.names = F)
DPP_Up_SH2_reg_freq = read.csv("./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_SH2_reg_freq.csv", header = T)
colnames(DPP_Up_SH2_reg_freq) = c("regulator", "No_of_targets")
DPP_Up_SH2_reg_freq = DPP_Up_SH2_reg_freq[(DPP_Up_SH2_reg_freq$No_of_targets)>0,]
DPP_Up_SH2_reg_freq = DPP_Up_SH2_reg_freq[!DPP_Up_SH2_reg_freq$regulator == "",]
DPP_Up_SH2_reg_percentage = find_ratio(regulation_freq, DPP_Up_SH2_reg_freq)
write.csv(DPP_Up_SH2_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/Up/DPP_SH2_reg_percentage.csv", row.names = FALSE)
View(DPP_Up_SH2_reg_percentage)

#DPP_SH5
DPP_SH5_geneRegulations = read.csv("./Data/combined/output/geneRegulations/old/regulations/DPP_SH5_geneRegulations.csv", header = T)
DPP_Up_SH5_geneRegulations = DPP_SH5_geneRegulations[(DPP_SH5_geneRegulations$logFC)>0,]
DPP_Up_SH5_geneRegulations = table(DPP_Up_SH5_geneRegulations$regulator)
write.csv(DPP_Up_SH5_geneRegulations, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_SH5_reg_freq.csv", row.names = F)
DPP_Up_SH5_reg_freq = read.csv("./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_SH5_reg_freq.csv", header = T)
colnames(DPP_Up_SH5_reg_freq) = c("regulator", "No_of_targets")
DPP_Up_SH5_reg_freq = DPP_Up_SH5_reg_freq[(DPP_Up_SH5_reg_freq$No_of_targets)>0,]
DPP_Up_SH5_reg_freq = DPP_Up_SH5_reg_freq[!DPP_Up_SH5_reg_freq$regulator == "",]
DPP_Up_SH5_reg_percentage = find_ratio(regulation_freq, DPP_Up_SH5_reg_freq)
write.csv(DPP_Up_SH5_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/Up/DPP_SH5_reg_percentage.csv", row.names = FALSE)
View(DPP_Up_SH5_reg_percentage)

#DPP_SH5_vs_SH2
DPP_SH5_vs_SH2_geneRegulations = read.csv("./Data/combined/output/geneRegulations/old/regulations/DPP_SH5_vs_SH2_geneRegulations.csv", header = T)
DPP_Up_SH5_vs_SH2_geneRegulations = DPP_SH5_vs_SH2_geneRegulations[(DPP_SH5_vs_SH2_geneRegulations$logFC)>0,]
DPP_Up_SH5_vs_SH2_geneRegulations = table(DPP_Up_SH5_vs_SH2_geneRegulations$regulator)
write.csv(DPP_Up_SH5_vs_SH2_geneRegulations, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_SH5_vs_SH2_reg_freq.csv", row.names = F)
DPP_Up_SH5_vs_SH2_reg_freq = read.csv("./Data/combined/output/geneRegulations/regulators/percentage_regulator/Up/DPP_Up_SH5_vs_SH2_reg_freq.csv", header = T)
colnames(DPP_Up_SH5_vs_SH2_reg_freq) = c("regulator", "No_of_targets")
DPP_Up_SH5_vs_SH2_reg_freq = DPP_Up_SH5_vs_SH2_reg_freq[(DPP_Up_SH5_vs_SH2_reg_freq$No_of_targets)>0,]
DPP_Up_SH5_vs_SH2_reg_freq = DPP_Up_SH5_vs_SH2_reg_freq[!DPP_Up_SH5_vs_SH2_reg_freq$regulator == "",]
DPP_Up_SH5_vs_SH2_reg_percentage = find_ratio(regulation_freq, DPP_Up_SH5_vs_SH2_reg_freq)
write.csv(DPP_Up_SH5_vs_SH2_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/Up/DPP_SH5_vs_SH2_reg_percentage.csv", row.names = FALSE)
View(DPP_Up_SH5_vs_SH2_reg_percentage)













