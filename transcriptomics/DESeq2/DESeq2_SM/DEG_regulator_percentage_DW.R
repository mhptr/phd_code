#22/08/19

setwd("~/OneDrive - University of Warwick/WORK/Results/Trans_Prot_COMBINED/Data/")

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

#DEG_M9
intersect_TP_M9_geneRegulations = read.csv("./Result/geneRegulations/intersect_TP_M9_geneRegulations.csv", header = T)
intersect_TP_M9_reg_freq = table(intersect_TP_M9_geneRegulations$regulator)
write.csv(intersect_TP_M9_reg_freq, "./Result/geneRegulations/regulators/percentage_regulator/intersect_TP_M9_reg_freq.csv", row.names = F)
intersect_TP_M9_reg_freq = read.csv("./Result/geneRegulations/regulators/percentage_regulator/intersect_TP_M9_reg_freq.csv", header = T)
colnames(intersect_TP_M9_reg_freq) = c("regulator", "No_of_targets")
intersect_TP_M9_reg_freq = intersect_TP_M9_reg_freq[!intersect_TP_M9_reg_freq$regulator == "",]
intersect_TP_M9_reg_percentage = find_ratio(regulation_freq, intersect_TP_M9_reg_freq)
write.csv(intersect_TP_M9_reg_percentage, "./Result/geneRegulations/regulators/percentage_regulator/output/intersect_TP_M9_reg_percentage.csv", row.names = FALSE)
View(intersect_TP_M9_reg_percentage)

#DEG_SH2
intersect_TP_SH2_geneRegulations = read.csv("./Result/geneRegulations/intersect_TP_SH2_geneRegulations.csv", header = T)
intersect_TP_SH2_reg_freq = table(intersect_TP_SH2_geneRegulations$regulator)
write.csv(intersect_TP_SH2_reg_freq, "./Result/geneRegulations/regulators/percentage_regulator/intersect_TP_SH2_reg_freq.csv", row.names = F)
intersect_TP_SH2_reg_freq = read.csv("./Result/geneRegulations/regulators/percentage_regulator/intersect_TP_SH2_reg_freq.csv", header = T)
colnames(intersect_TP_SH2_reg_freq) = c("regulator", "No_of_targets")
intersect_TP_SH2_reg_freq = intersect_TP_SH2_reg_freq[!intersect_TP_SH2_reg_freq$regulator == "",]
intersect_TP_SH2_reg_percentage = find_ratio(regulation_freq, intersect_TP_SH2_reg_freq)
write.csv(intersect_TP_SH2_reg_percentage, "./Result/geneRegulations/regulators/percentage_regulator/output/intersect_TP_SH2_reg_percentage.csv", row.names = FALSE)
View(intersect_TP_SH2_reg_percentage)

#DEG_SH5
intersect_TP_SH5_geneRegulations = read.csv("./Result/geneRegulations/intersect_TP_SH5_geneRegulations.csv", header = T)
intersect_TP_SH5_reg_freq = table(intersect_TP_SH5_geneRegulations$regulator)
write.csv(intersect_TP_SH5_reg_freq, "./Result/geneRegulations/regulators/percentage_regulator/intersect_TP_SH5_reg_freq.csv", row.names = F)
intersect_TP_SH5_reg_freq = read.csv("./Result/geneRegulations/regulators/percentage_regulator/intersect_TP_SH5_reg_freq.csv", header = T)
colnames(intersect_TP_SH5_reg_freq) = c("regulator", "No_of_targets")
intersect_TP_SH5_reg_freq = intersect_TP_SH5_reg_freq[!intersect_TP_SH5_reg_freq$regulator == "",]
intersect_TP_SH5_reg_percentage = find_ratio(regulation_freq, intersect_TP_SH5_reg_freq)
write.csv(intersect_TP_SH5_reg_percentage, "./Result/geneRegulations/regulators/percentage_regulator/output/intersect_TP_SH5_reg_percentage.csv", row.names = FALSE)
View(intersect_TP_SH5_reg_percentage)

#DEG_SH5_vs_SH2
intersect_TP_SH5_vs_SH2_geneRegulations = read.csv("./Result/geneRegulations/intersect_TP_SH5_vs_SH2_geneRegulations.csv", header = T)
intersect_TP_SH5_vs_SH2_reg_freq = table(intersect_TP_SH5_vs_SH2_geneRegulations$regulator)
write.csv(intersect_TP_SH5_vs_SH2_reg_freq, "./Result/geneRegulations/regulators/percentage_regulator/intersect_TP_SH5_vs_SH2_reg_freq.csv", row.names = F)
intersect_TP_SH5_vs_SH2_reg_freq = read.csv("./Result/geneRegulations/regulators/percentage_regulator/intersect_TP_SH5_vs_SH2_reg_freq.csv", header = T)
colnames(intersect_TP_SH5_vs_SH2_reg_freq) = c("regulator", "No_of_targets")
intersect_TP_SH5_vs_SH2_reg_freq = intersect_TP_SH5_vs_SH2_reg_freq[!intersect_TP_SH5_vs_SH2_reg_freq$regulator == "",]
intersect_TP_SH5_vs_SH2_reg_percentage = find_ratio(regulation_freq, intersect_TP_SH5_vs_SH2_reg_freq)
write.csv(intersect_TP_SH5_vs_SH2_reg_percentage, "./Result/geneRegulations/regulators/percentage_regulator/output/intersect_TP_SH5_vs_SH2_reg_percentage.csv", row.names = FALSE)
View(intersect_TP_SH5_vs_SH2_reg_percentage)



























