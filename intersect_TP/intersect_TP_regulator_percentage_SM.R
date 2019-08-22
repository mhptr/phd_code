#22/08/19

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
DPP_LB_reg_freq = read.csv("./Data/combined/output/geneRegulations/old/regulators/csv/DPP_LB_regulators.csv", header = T)
colnames(DPP_LB_reg_freq) = c("regulator", "No_of_targets")
DPP_LB_reg_freq = DPP_LB_reg_freq[!DPP_LB_reg_freq$regulator == "",]
DPP_LB_reg_percentage = find_ratio(regulation_freq, DPP_LB_reg_freq)
write.csv(DPP_LB_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/DPP_LB_regulators.csv", row.names = FALSE)
View(DPP_LB_reg_percentage)

#DPP_M9
DPP_M9_reg_freq = read.csv("./Data/combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", header = T)
colnames(DPP_M9_reg_freq) = c("regulator", "No_of_targets")
DPP_M9_reg_freq = DPP_M9_reg_freq[!DPP_M9_reg_freq$regulator == "",]
DPP_M9_reg_percentage = find_ratio(regulation_freq, DPP_M9_reg_freq)
write.csv(DPP_M9_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/DPP_M9_regulators.csv", row.names = FALSE)
View(DPP_M9_reg_percentage)

#DPP_SH2
DPP_SH2_reg_freq = read.csv("./Data/combined/output/geneRegulations/old/regulators/csv/DPP_SH2_regulators.csv", header = T)
colnames(DPP_SH2_reg_freq) = c("regulator", "No_of_targets")
DPP_SH2_reg_freq = DPP_SH2_reg_freq[!DPP_SH2_reg_freq$regulator == "",]
DPP_SH2_reg_percentage = find_ratio(regulation_freq, DPP_SH2_reg_freq)
write.csv(DPP_SH2_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/DPP_SH2_regulators.csv", row.names = FALSE)
View(DPP_SH2_reg_percentage)

#DPP_SH5
DPP_SH5_reg_freq = read.csv("./Data/combined/output/geneRegulations/old/regulators/csv/DPP_SH5_regulators.csv", header = T)
colnames(DPP_SH5_reg_freq) = c("regulator", "No_of_targets")
DPP_SH5_reg_freq = DPP_SH5_reg_freq[!DPP_SH5_reg_freq$regulator == "",]
DPP_SH5_reg_percentage = find_ratio(regulation_freq, DPP_SH5_reg_freq)
write.csv(DPP_SH5_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/DPP_SH5_regulators.csv", row.names = FALSE)
View(DPP_SH5_reg_percentage)

#DPP_SH5_vs_SH2
DPP_SH5_vs_SH2_reg_freq = read.csv("./Data/combined/output/geneRegulations/old/regulators/csv/DPP_SH5_vs_SH2_regulators.csv", header = T)
colnames(DPP_SH5_vs_SH2_reg_freq) = c("regulator", "No_of_targets")
DPP_SH5_vs_SH2_reg_freq = DPP_SH5_vs_SH2_reg_freq[!DPP_SH5_vs_SH2_reg_freq$regulator == "",]
DPP_SH5_vs_SH2_reg_percentage = find_ratio(regulation_freq, DPP_SH5_vs_SH2_reg_freq)
write.csv(DPP_SH5_vs_SH2_reg_percentage, "./Data/combined/output/geneRegulations/regulators/percentage_regulator/output/DPP_SH5_vs_SH2_regulators.csv", row.names = FALSE)
View(DPP_SH5_vs_SH2_reg_percentage)


