#22/08/19

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/")


# 1. make a frequency table for each regulator from the gene regulation file
# 2. calculate the number of regulators found in the differential data set 
# 3. calculate the percentage of the regulators changed from the master gene regulation file


create_freq_table = function(csv, colname) {
  # print(csv)
  print(colname)
  freq = table(csv$colname)
  return(freq)
}
a = create_freq_table(DEG_M9_geneRegulations, regulator)



DEG_M9_regulators = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/frequency_table/DEG_M9_regulators_frequency.csv", header = T)

regulation = read.csv("./SubtiWiki Exports /regulations.csv", header = T)
regulation_freq = table(regulation$regulator)
write.csv(regulation_freq, "./SubtiWiki Exports /regulation_freq.csv", row.names = F)

regulation_freq = read.csv("./SubtiWiki Exports /regulation_freq.csv", header = T)
colnames(regulation_freq) = c("regulator", "frequency")
dim(regulation_freq)
regulation_freq = regulation_freq[!regulation_freq$regulator == "",]
dim(regulation_freq)

DEG_M9_regulators_regulation_freq = merge(DEG_M9_regulators, regulation_freq, by.x="regulator", by.y="regulator", all.x = TRUE)
colnames(DEG_M9_regulators_regulation_freq) = c("regulator", "DEG_M9_frequency", "regulation_freq")

DEG_M9_regulators_regulation_freq = transform(DEG_M9_regulators_regulation_freq, DEG_M9_regulators_percentage = 100 * DEG_M9_regulators_regulation_freq$DEG_M9_frequency/DEG_M9_regulators_regulation_freq$regulation_freq)
View(DEG_M9_regulators_regulation_freq)






























