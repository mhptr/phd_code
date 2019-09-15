#06/09/19

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/")

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
DEG_M9_geneRegulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_M9_geneRegulations.csv", header = T)
DEG_Down_M9_geneRegulations = DEG_M9_geneRegulations[(DEG_M9_geneRegulations$log2FoldChange)<0,]
DEG_Down_M9_geneRegulations = table(DEG_Down_M9_geneRegulations$regulator)
write.csv(DEG_Down_M9_geneRegulations, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/Down/DEG_Down_M9_reg_freq.csv", row.names = F)
DEG_Down_M9_reg_freq = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/Down/DEG_Down_M9_reg_freq.csv", header = T)
colnames(DEG_Down_M9_reg_freq) = c("regulator", "No_of_targets")
DEG_Down_M9_reg_freq = DEG_Down_M9_reg_freq[(DEG_Down_M9_reg_freq$No_of_targets)>0,]
DEG_Down_M9_reg_freq = DEG_Down_M9_reg_freq[!DEG_Down_M9_reg_freq$regulator == "",]
DEG_Down_M9_reg_percentage = find_ratio(regulation_freq, DEG_Down_M9_reg_freq)
write.csv(DEG_Down_M9_reg_percentage, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/output/Down/DEG_Down_M9_reg_percentage.csv", row.names = FALSE)
View(DEG_Down_M9_reg_percentage)

#DEG_SH2
DEG_SH2_geneRegulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_SH2_geneRegulations.csv", header = T)
DEG_Down_SH2_geneRegulations = DEG_SH2_geneRegulations[(DEG_SH2_geneRegulations$log2FoldChange)<0,]
DEG_Down_SH2_geneRegulations = table(DEG_Down_SH2_geneRegulations$regulator)
write.csv(DEG_Down_SH2_geneRegulations, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/Down/DEG_Down_SH2_reg_freq.csv", row.names = F)
DEG_Down_SH2_reg_freq = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/Down/DEG_Down_SH2_reg_freq.csv", header = T)
colnames(DEG_Down_SH2_reg_freq) = c("regulator", "No_of_targets")
DEG_Down_SH2_reg_freq = DEG_Down_SH2_reg_freq[(DEG_Down_SH2_reg_freq$No_of_targets)>0,]
DEG_Down_SH2_reg_freq = DEG_Down_SH2_reg_freq[!DEG_Down_SH2_reg_freq$regulator == "",]
DEG_Down_SH2_reg_percentage = find_ratio(regulation_freq, DEG_Down_SH2_reg_freq)
write.csv(DEG_Down_SH2_reg_percentage, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/output/Down/DEG_Down_SH2_reg_percentage.csv", row.names = FALSE)
View(DEG_Down_SH2_reg_percentage)

#DEG_SH5
DEG_SH5_geneRegulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_SH5_geneRegulations.csv", header = T)
DEG_Down_SH5_geneRegulations = DEG_SH5_geneRegulations[(DEG_SH5_geneRegulations$log2FoldChange)<0,]
DEG_Down_SH5_geneRegulations = table(DEG_Down_SH5_geneRegulations$regulator)
write.csv(DEG_Down_SH5_geneRegulations, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/Down/DEG_Down_SH5_reg_freq.csv", row.names = F)
DEG_Down_SH5_reg_freq = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/Down/DEG_Down_SH5_reg_freq.csv", header = T)
colnames(DEG_Down_SH5_reg_freq) = c("regulator", "No_of_targets")
DEG_Down_SH5_reg_freq = DEG_Down_SH5_reg_freq[(DEG_Down_SH5_reg_freq$No_of_targets)>0,]
DEG_Down_SH5_reg_freq = DEG_Down_SH5_reg_freq[!DEG_Down_SH5_reg_freq$regulator == "",]
DEG_Down_SH5_reg_percentage = find_ratio(regulation_freq, DEG_Down_SH5_reg_freq)
write.csv(DEG_Down_SH5_reg_percentage, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/output/Down/DEG_Down_SH5_reg_percentage.csv", row.names = FALSE)
View(DEG_Down_SH5_reg_percentage)

#DEG_SH5_vs_SH2
DEG_SH5_vs_SH5_geneRegulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_SH5_vs_SH2_geneRegulations.csv", header = T)
DEG_Down_SH5_vs_SH5_geneRegulations = DEG_SH5_vs_SH5_geneRegulations[(DEG_SH5_vs_SH5_geneRegulations$log2FoldChange)<0,]
DEG_Down_SH5_vs_SH5_geneRegulations = table(DEG_Down_SH5_vs_SH5_geneRegulations$regulator)
write.csv(DEG_Down_SH5_vs_SH5_geneRegulations, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/Down/DEG_Down_SH5_vs_SH5_reg_freq.csv", row.names = F)
DEG_Down_SH5_vs_SH5_reg_freq = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/Down/DEG_Down_SH5_vs_SH5_reg_freq.csv", header = T)
colnames(DEG_Down_SH5_vs_SH5_reg_freq) = c("regulator", "No_of_targets")
DEG_Down_SH5_vs_SH5_reg_freq = DEG_Down_SH5_vs_SH5_reg_freq[(DEG_Down_SH5_vs_SH5_reg_freq$No_of_targets)>0,]
DEG_Down_SH5_vs_SH5_reg_freq = DEG_Down_SH5_vs_SH5_reg_freq[!DEG_Down_SH5_vs_SH5_reg_freq$regulator == "",]
DEG_Down_SH5_vs_SH5_reg_percentage = find_ratio(regulation_freq, DEG_Down_SH5_vs_SH5_reg_freq)
write.csv(DEG_Down_SH5_vs_SH5_reg_percentage, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/output/Down/DEG_Down_SH5_vs_SH5_reg_percentage.csv", row.names = FALSE)
View(DEG_Down_SH5_vs_SH5_reg_percentage)



























