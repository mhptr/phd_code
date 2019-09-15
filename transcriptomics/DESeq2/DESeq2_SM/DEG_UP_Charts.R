#making graphs and charts of DEG
#05/09/19

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/")
# 
# DEG_M9_regulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_M9_geneRegulations.csv", header = T)
# DEG_SH2_regulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_SH2_geneRegulations.csv", header = T)
# DEG_SH5_regulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_SH5_geneRegulations.csv", header = T)
# DEG_SH5_vs_SH2_regulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_SH5_vs_SH2_geneRegulations.csv", header = T)


# dim(DEG_M9_regulations)
# DEG_M9_regulations = DEG_M9_regulations[!DEG_M9_regulations$regulator == "",]
# dim(DEG_M9_regulations)
# write.csv(DEG_M9_regulations, "./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_M9_regulations.csv", row.names = F)
# 
# dim(DEG_SH2_regulations)
# DEG_SH2_regulations = DEG_SH2_regulations[!DEG_SH2_regulations$regulator == "",]
# dim(DEG_SH2_regulations)
# write.csv(DEG_SH2_regulations, "./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_SH2_regulations.csv", row.names = F)
# 
# dim(DEG_SH5_regulations)
# DEG_SH5_regulations = DEG_SH5_regulations[!DEG_SH5_regulations$regulator == "",]
# dim(DEG_SH5_regulations)
# write.csv(DEG_SH5_regulations, "./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_SH5_regulations.csv", row.names = F)
# 
# dim(DEG_SH5_vs_SH2_regulations)
# DEG_SH5_vs_SH2_regulations = DEG_SH5_vs_SH2_regulations[!DEG_SH5_vs_SH2_regulations$regulator == "",]
# dim(DEG_SH5_vs_SH2_regulations)
# write.csv(DEG_SH5_vs_SH2_regulations, "./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_SH5_vs_SH2_regulations.csv", row.names = F)


#M9
DEG_M9_regulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_M9_regulations.csv", header = T)
DEG_M9_regulations_Up = DEG_M9_regulations[(DEG_M9_regulations$log2FoldChange)>0,]
DEG_M9_reg_freq_Up = table(DEG_M9_regulations_Up$regulator)
write.csv(DEG_M9_reg_freq_Up, "./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_regulation_freq_Up/DEG_M9_reg_freq_Up.csv", row.names = F)
DEG_M9_reg_freq_Up = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_regulation_freq_Up/DEG_M9_reg_freq_Up.csv", header = T)
colnames(DEG_M9_reg_freq_Up) = c("regulator", "frequency")
dim(DEG_M9_reg_freq_Up)
DEG_M9_reg_freq_Up = DEG_M9_reg_freq_Up[!DEG_M9_reg_freq_Up$regulator == "",]
dim(DEG_M9_reg_freq_Up)
DEG_M9_reg_freq_Up = DEG_M9_reg_freq_Up[(DEG_M9_reg_freq_Up$frequency)>0,]
write.csv(DEG_M9_reg_freq_Up, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_reg_freq_Up/DEG_M9_reg_freq_Up.csv", row.names = F)
View(DEG_M9_reg_freq_Up)

#SH2
DEG_SH2_regulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_SH2_regulations.csv", header = T)
DEG_SH2_regulations_Up = DEG_SH2_regulations[(DEG_SH2_regulations$log2FoldChange)>0,]
DEG_SH2_reg_freq_Up = table(DEG_SH2_regulations_Up$regulator)
write.csv(DEG_SH2_reg_freq_Up, "./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_regulation_freq_Up/DEG_SH2_reg_freq_Up.csv", row.names = F)
DEG_SH2_reg_freq_Up = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_regulation_freq_Up/DEG_SH2_reg_freq_Up.csv", header = T)
colnames(DEG_SH2_reg_freq_Up) = c("regulator", "frequency")
dim(DEG_SH2_reg_freq_Up)
DEG_SH2_reg_freq_Up = DEG_SH2_reg_freq_Up[!DEG_SH2_reg_freq_Up$regulator == "",]
dim(DEG_SH2_reg_freq_Up)
DEG_SH2_reg_freq_Up = DEG_SH2_reg_freq_Up[(DEG_SH2_reg_freq_Up$frequency)>0,]
write.csv(DEG_SH2_reg_freq_Up, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_reg_freq_Up/DEG_SH2_reg_freq_Up.csv", row.names = F)
View(DEG_SH2_reg_freq_Up)

#SH5
DEG_SH5_regulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_SH5_regulations.csv", header = T)
DEG_SH5_regulations_Up = DEG_SH5_regulations[(DEG_SH5_regulations$log2FoldChange)>0,]
DEG_SH5_reg_freq_Up = table(DEG_SH5_regulations_Up$regulator)
write.csv(DEG_SH5_reg_freq_Up, "./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_regulation_freq_Up/DEG_SH5_reg_freq_Up.csv", row.names = F)
DEG_SH5_reg_freq_Up = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_regulation_freq_Up/DEG_SH5_reg_freq_Up.csv", header = T)
colnames(DEG_SH5_reg_freq_Up) = c("regulator", "frequency")
dim(DEG_SH5_reg_freq_Up)
DEG_SH5_reg_freq_Up = DEG_SH5_reg_freq_Up[!DEG_SH5_reg_freq_Up$regulator == "",]
dim(DEG_SH5_reg_freq_Up)
DEG_SH5_reg_freq_Up = DEG_SH5_reg_freq_Up[(DEG_SH5_reg_freq_Up$frequency)>0,]
write.csv(DEG_SH5_reg_freq_Up, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_reg_freq_Up/DEG_SH5_reg_freq_Up.csv", row.names = F)
View(DEG_SH5_reg_freq_Up)

#SH5_vs_SH2
DEG_SH5_vs_SH2_regulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_SH5_vs_SH2_regulations.csv", header = T)
DEG_SH5_vs_SH2_regulations_Up = DEG_SH5_vs_SH2_regulations[(DEG_SH5_vs_SH2_regulations$log2FoldChange)>0,]
DEG_SH5_vs_SH2_reg_freq_Up = table(DEG_SH5_vs_SH2_regulations_Up$regulator)
write.csv(DEG_SH5_vs_SH2_reg_freq_Up, "./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_regulation_freq_Up/DEG_SH5_vs_SH2_reg_freq_Up.csv", row.names = F)
DEG_SH5_vs_SH2_reg_freq_Up = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulations/csv/DEG_regulation_freq_Up/DEG_SH5_vs_SH2_reg_freq_Up.csv", header = T)
colnames(DEG_SH5_vs_SH2_reg_freq_Up) = c("regulator", "frequency")
dim(DEG_SH5_vs_SH2_reg_freq_Up)
DEG_SH5_vs_SH2_reg_freq_Up = DEG_SH5_vs_SH2_reg_freq_Up[!DEG_SH5_vs_SH2_reg_freq_Up$regulator == "",]
dim(DEG_SH5_vs_SH2_reg_freq_Up)
DEG_SH5_vs_SH2_reg_freq_Up = DEG_SH5_vs_SH2_reg_freq_Up[(DEG_SH5_vs_SH2_reg_freq_Up$frequency)>0,]
write.csv(DEG_SH5_vs_SH2_reg_freq_Up, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_reg_freq_Up/DEG_SH5_vs_SH2_reg_freq_Up.csv", row.names = F)
View(DEG_SH5_vs_SH2_reg_freq_Up)


#frequency tables
DEG_M9_reg_freq_Up = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_reg_freq_Up/DEG_M9_reg_freq_Up.csv", header = T)
DEG_SH2_reg_freq_Up = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_reg_freq_Up/DEG_SH2_reg_freq_Up.csv", header = T)
DEG_SH5_reg_freq_Up = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_reg_freq_Up/DEG_SH5_reg_freq_Up.csv", header = T)
DEG_SH5_vs_SH2_reg_freq_Up = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_reg_freq_Up/DEG_SH5_vs_SH2_reg_freq_Up.csv", header = T)


library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(scales)

swati.color = c(brewer.pal(8,"Dark2"))

#http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization

# Barplot

#M9
DEG_M9_reg_freq_Up_mod = DEG_M9_reg_freq_Up[(DEG_M9_reg_freq_Up$frequency)>2,]
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_Up_barplot/barplot_DEG_M9_reg_Up.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DEG_M9_reg_freq_Up_mod = ggplot(DEG_M9_reg_freq_Up_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DEG M9 reg Up") + 
  xlab("Regulons affected due to del_spoVG in M9") + ylab("No. of genes in each regulon") +
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
DEG_M9_reg_freq_Up_mod
coord_fixed()
dev.off()  
#Create a pie chart 
# pie_M9 = regulators_M9 + coord_polar("y", start=0) + geom_text(aes(y = frequency/length(DEG_M9_regulators_mod$regulator) + c(0, cumsum(frequency)[-length(frequency)]), 
#                                                                    label = percent(frequency/100)), size=5)
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_Up_pie/pie_DEG_M9_reg_Up.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_M9 = DEG_M9_reg_freq_Up_mod + coord_polar("y", start=0)
pie_M9
coord_fixed()
dev.off() 


#SH2
DEG_SH2_reg_freq_Up_mod = DEG_SH2_reg_freq_Up[(DEG_SH2_reg_freq_Up$frequency)>2,]
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_Up_barplot/barplot_DEG_SH2_reg_Up.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DEG_SH2_reg_freq_Up_mod = ggplot(DEG_SH2_reg_freq_Up_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DEG SH2 reg Up") + 
  xlab("Regulons affected due to del_spoVG in SH2") + ylab("No. of genes in each regulon") +
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
DEG_SH2_reg_freq_Up_mod
coord_fixed()
dev.off()  
#Create a pie chart 
# pie_SH2 = regulators_SH2 + coord_polar("y", start=0) + geom_text(aes(y = frequency/length(DEG_SH2_regulators_mod$regulator) + c(0, cumsum(frequency)[-length(frequency)]), 
#                                                                    label = percent(frequency/100)), size=5)
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_Up_pie/pie_DEG_SH2_reg_Up.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_SH2 = DEG_SH2_reg_freq_Up_mod + coord_polar("y", start=0)
pie_SH2
coord_fixed()
dev.off() 


#SH5
DEG_SH5_reg_freq_Up_mod = DEG_SH5_reg_freq_Up[(DEG_SH5_reg_freq_Up$frequency)>2,]
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_Up_barplot/barplot_DEG_SH5_reg_Up.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DEG_SH5_reg_freq_Up_mod = ggplot(DEG_SH5_reg_freq_Up_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DEG SH5 reg Up") + 
  xlab("Regulons affected due to del_spoVG in SH5") + ylab("No. of genes in each regulon") +
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
DEG_SH5_reg_freq_Up_mod
coord_fixed()
dev.off()  
#Create a pie chart 
# pie_SH5 = regulators_SH5 + coord_polar("y", start=0) + geom_text(aes(y = frequency/length(DEG_SH5_regulators_mod$regulator) + c(0, cumsum(frequency)[-length(frequency)]), 
#                                                                    label = percent(frequency/100)), size=5)
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_Up_pie/pie_DEG_SH5_reg_Up.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_SH5 = DEG_SH5_reg_freq_Up_mod + coord_polar("y", start=0)
pie_SH5
coord_fixed()
dev.off() 


#SH5_vs_SH2
DEG_SH5_vs_SH2_reg_freq_Up_mod = DEG_SH5_vs_SH2_reg_freq_Up[(DEG_SH5_vs_SH2_reg_freq_Up$frequency)>18,]
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_Up_barplot/barplot_DEG_SH5_vs_SH2_reg_Up.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DEG_SH5_vs_SH2_reg_freq_Up_mod = ggplot(DEG_SH5_vs_SH2_reg_freq_Up_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DEG SH5_vs_SH2 reg Up") + 
  xlab("Regulons affected due to del_spoVG in SH5_vs_SH2") + ylab("No. of genes in each regulon") +
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
DEG_SH5_vs_SH2_reg_freq_Up_mod
coord_fixed()
dev.off()  
#Create a pie chart 
# pie_SH5_vs_SH2 = regulators_SH5_vs_SH2 + coord_polar("y", start=0) + geom_text(aes(y = frequency/length(DEG_SH5_vs_SH2_regulators_mod$regulator) + c(0, cumsum(frequency)[-length(frequency)]), 
#                                                                    label = percent(frequency/100)), size=5)
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_Up_pie/pie_DEG_SH5_vs_SH2_reg_Up.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_SH5_vs_SH2 = DEG_SH5_vs_SH2_reg_freq_Up_mod + coord_polar("y", start=0)
pie_SH5_vs_SH2
coord_fixed()
dev.off()  


#barplot_all
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_Up_barplot/barplot_DEG_SH5_reg_freq_Up_mod_all.pdf" ,width=5, height=5)      #turn this OFF if just want to see the picture in the Plots
barplot_DEG_SH5_reg_freq_Up_mod_all = grid.arrange(DEG_M9_reg_freq_Up_mod, DEG_SH2_reg_freq_Up_mod,
                           DEG_SH5_reg_freq_Up_mod, DEG_SH5_vs_SH2_reg_freq_Up_mod, nrow = 6)
barplot_DEG_SH5_reg_freq_Up_mod_all
coord_fixed()
dev.off() 

#pie_all
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_Up_pie/pie_DEG_SH5_reg_freq_Up_mod_all.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_DEG_SH5_reg_freq_Up_mod_all = grid.arrange(pie_M9, pie_SH2,
                       pie_SH5, pie_SH5_vs_SH2, nrow = 4)
coord_fixed()
dev.off() 

















############################################################################################################################

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



























############################################################################################################################
############################################################################################################################
############################################################################################################################


hist(DEG_LB_regulators,
     main="Regulons affected due to del_spoVG in LB",
     xlab="Regulators",
     xlim=c(50,100),
     col="No. of genes in each regulon",
     freq=FALSE
)


plot(DEG_LB_regulators)

hist(DEG_LB_regulators$frequency)

# Add a Normal Curve (Thanks to Peter Dalgaard)
x <- DEG_LB_regulators$frequency
h<-hist(x, breaks=10, col="pink", xlab="Regulons affected due to del_spoVG in LB", 
        main="Histogram with Normal Curve", 
        ylab = "No. of genes in each regulon") 
xfit<-seq(min(x),max(x),length=100) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="red", lwd=2)

# Add all the labels and color. Set the y axis indexes horizontal. Set limits for axis and
# Drawing lables on top of bars
hist(DEG_LB_regulators$frequency,main="Regulons affected due to del_spoVG in LB", xlab = "Regulons affected due to del_spoVG in LB", ylab = "No. of genes in each regulon",border="red", col="blue",las=1,xlim=c(1,40),ylim=c(0,60),labels=TRUE)
hist(DEG_SH5_regulators$frequency,main="Regulons affected due to del_spoVG in LB", xlab = "Regulons affected due to del_spoVG in LB", ylab = "No. of genes in each regulon",border="red", col="blue",las=1,xlim=c(1,100),ylim=c(0,100),labels=TRUE)




# #http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization
# library(ggplot2) 
# library(gridExtra)
# p1 <- ggplot(DEG_LB_regulators, aes(x = frequency)) +
#   geom_histogram(binwidth = 1000, color = "grey30", fill = "white") +
#   ggtitle("Bin Width = 1000")
# p2 <- ggplot(DEG_M9_regulators, aes(x = Freq)) +
#   geom_histogram(binwidth = 5000, color = "grey30", fill = "white") +
#   ggtitle("Bin Width = 5000")
# p3 <- ggplot(DEG_SH2_regulators, aes(x = frequency)) +
#   geom_histogram(binwidth = 10000, color = "grey30", fill = "white") +
#   ggtitle("Bin Width = 10000")
# grid.arrange(p1, p2, p3, ncol=3)
# library(ggplot2)
# # Basic histogram
# ggplot(DEG_LB_regulators, aes(x=frequency)) + geom_histogram()
# # Change the width of bins
# ggplot(DEG_LB_regulators, aes(x=frequency)) + 
#   geom_histogram(binwidth=1)
# # Change colors
# p<-ggplot(DEG_LB_regulators, aes(x=frequency)) + 
#   geom_histogram(color="black", fill="white")
# p
# library(plyr)
# mu <- ddply(DEG_LB_regulators, "regulator", summarise, grp.mean=mean(frequency))
# head(mu)
# # Change histogram plot line colors by groUps
# ggplot(DEG_M9_regulators, aes(x=Freq, color=regulator)) +
#   geom_histogram(fill="white")
# # Overlaid histograms
# ggplot(DEG_LB_regulators, aes(x=frequency, color=regulator)) +
#   geom_histogram(fill="white", alpha=0.5, position="identity")
# DEG_M9_regulators = DEG_M9_regulators[(DEG_M9_regulators$Freq)>=3,]
# dim(DEG_M9_regulators)
# rownames(DEG_M9_regulators) = FALSE
# hist(DEG_M9_regulators$Var1)
# ggplot(DEG_M9_regulators, aes(x=Freq, color=Var1)) +
#   geom_histogram(fill="white")






# Apply blank theme
library(scales)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )



pie + scale_fill_grey() +  blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                label = percent(value/100)), size=5)




# Use brewer palette
pie + scale_fill_brewer("Blues") + blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                label = percent(value/100)), size=5)














