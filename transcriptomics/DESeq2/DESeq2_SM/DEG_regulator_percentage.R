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
write.csv(DEG_M9_regulators_regulation_freq, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/percentage_regulator/DEG_M9_regulators_regulation_freq.csv", row.names = F)
View(DEG_M9_regulators_regulation_freq)


#M9 
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/percentage_regulator/barplot/DEG_barplot_regulators_M9_freq.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DEG_M9_regulators_regulation_freq = DEG_M9_regulators_regulation_freq[(DEG_M9_regulators_regulation_freq$DEG_M9_regulators_percentage)>75,]
DEG_M9_regulators_regulation_freq = ggplot(DEG_M9_regulators_regulation_freq, aes(x="", y=DEG_M9_regulators_percentage, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DEG Regulators \n M9") + 
  xlab("Regulons affected due to del_spoVG in M9") + ylab("No. of genes in each regulon") +
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
DEG_M9_regulators_regulation_freq
coord_fixed()
dev.off() 
#Create a pie chart 
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/percentage_regulator/pie/DEG_pie_regulators_M9_freq.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_M9 <- DEG_M9_regulators_regulation_freq + coord_polar("y", start=0)
pie_M9
coord_fixed()
dev.off() 





####################################################





















