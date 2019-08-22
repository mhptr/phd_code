#22/08/19

setwd("~/OneDrive - University of Warwick/WORK/Results/Transcriptomic/intersection_toRNAdo_all/")

library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(scales)

DEG_M9_geneRegulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_M9_geneRegulations.csv", header = T)
DEG_SH2_geneRegulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_SH2_geneRegulations.csv", header = T)
DEG_SH5_geneRegulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_SH5_geneRegulations.csv", header = T)
DEG_SH5_vs_SH2_geneRegulations = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/DEG_SH5_vs_SH2_geneRegulations.csv", header = T)

DEG_M9_regulators = table(DEG_M9_geneRegulations$regulator)
DEG_SH2_regulators = table(DEG_SH2_geneRegulations$regulator)
DEG_SH5_regulators = table(DEG_SH5_geneRegulations$regulator)
DEG_SH5_vs_SH2_regulators = table(DEG_SH5_vs_SH2_geneRegulations$regulator)

#DEG_M9_regulators = DEG_M9_regulators[!DEG_M9_regulators$Name == "",] #removes the empty row
write.csv(DEG_M9_regulators, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_M9_regulators.csv", row.names = F)
write.csv(DEG_SH2_regulators, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_SH2_regulators.csv", row.names = F)
write.csv(DEG_SH5_regulators, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_SH5_regulators.csv", row.names = F)
write.csv(DEG_SH5_vs_SH2_regulators, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_SH5_vs_SH2_regulators.csv", row.names = F)

#M9
dim(DEG_M9_regulators)
DEG_M9_regulators = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_M9_regulators.csv", header = T)
colnames(DEG_M9_regulators) = c("regulator", "frequency")
dim(DEG_M9_regulators)
DEG_M9_regulators = DEG_M9_regulators[!DEG_M9_regulators$regulator == "",]
dim(DEG_M9_regulators)
# DEG_M9_regulators = DEG_M9_regulators[, order(DEG_M9_regulators$frequency)]
write.csv(DEG_M9_regulators, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/frequency_table/DEG_M9_regulators_frequency.csv", row.names = F)
View(DEG_M9_regulators)

#SH2
dim(DEG_SH2_regulators)
DEG_SH2_regulators = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_SH2_regulators.csv", header = T)
colnames(DEG_SH2_regulators) = c("regulator", "frequency")
dim(DEG_SH2_regulators)
DEG_SH2_regulators = DEG_SH2_regulators[!DEG_SH2_regulators$regulator == "",]
dim(DEG_SH2_regulators)
# DEG_SH2_regulators = DEG_SH2_regulators[, order(DEG_SH2_regulators$frequency)]
write.csv(DEG_SH2_regulators, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/frequency_table/DEG_SH2_regulators_frequency.csv", row.names = F)
View(DEG_SH2_regulators)

#SH5
dim(DEG_SH5_regulators)
DEG_SH5_regulators = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_SH5_regulators.csv", header = T)
colnames(DEG_SH5_regulators) = c("regulator", "frequency")
dim(DEG_SH5_regulators)
DEG_SH5_regulators = DEG_SH5_regulators[!DEG_SH5_regulators$regulator == "",]
dim(DEG_SH5_regulators)
# DEG_SH5_regulators = DEG_SH5_regulators[, order(DEG_SH5_regulators$frequency)]
write.csv(DEG_SH5_regulators, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/frequency_table/DEG_SH5_regulators_frequency.csv", row.names = F)
View(DEG_SH5_regulators)

#SH5_vs_SH2
dim(DEG_SH5_vs_SH2_regulators)
DEG_SH5_vs_SH2_regulators = read.csv("./diff.exp.gene/DEG_SM/geneRegulations/regulators/csv/DEG_SH5_vs_SH2_regulators.csv", header = T)
colnames(DEG_SH5_vs_SH2_regulators) = c("regulator", "frequency")
dim(DEG_SH5_vs_SH2_regulators)
DEG_SH5_vs_SH2_regulators = DEG_SH5_vs_SH2_regulators[!DEG_SH5_vs_SH2_regulators$regulator == "",]
dim(DEG_SH5_vs_SH2_regulators)
# DEG_SH5_vs_SH2_regulators = DEG_SH5_vs_SH2_regulators[, order(DEG_SH5_vs_SH2_regulators$frequency)]
write.csv(DEG_SH5_vs_SH2_regulators, "./diff.exp.gene/DEG_SM/geneRegulations/regulators/frequency_table/DEG_SH5_vs_SH2_regulators_frequency.csv", row.names = F)
View(DEG_SH5_vs_SH2_regulators)


# Barplot

#M9
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_barplot_regulators_M9.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DEG_M9_regulators_mod = DEG_M9_regulators[(DEG_M9_regulators$frequency)>3,]
regulators_M9 = ggplot(DEG_M9_regulators_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DEG Regulators \n M9") + 
  xlab("Regulons affected due to del_spoVG in M9") + ylab("No. of genes in each regulon") +
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
regulators_M9
coord_fixed()
dev.off() 
#Create a pie chart 
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_pie_regulators_M9.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_M9 <- regulators_M9 + coord_polar("y", start=0)
pie_M9
coord_fixed()
dev.off() 

#SH2
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_barplot_regulators_SH2.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DEG_SH2_regulators_mod = DEG_SH2_regulators[(DEG_SH2_regulators$frequency)>2,]
regulators_SH2 = ggplot(DEG_SH2_regulators_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DEG Regulators \n SH2") + 
  xlab("Regulons affected due to del_spoVG in SH2") + ylab("No. of genes in each regulon") +
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
regulators_SH2
coord_fixed()
dev.off() 
#Create a pie chart 
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_pie_regulators_SH2.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_SH2 <- regulators_SH2 + coord_polar("y", start=0)
pie_SH2
coord_fixed()
dev.off() 


#SH5
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_barplot_regulators_SH5.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DEG_SH5_regulators_mod = DEG_SH5_regulators[(DEG_SH5_regulators$frequency)>10,]
regulators_SH5 = ggplot(DEG_SH5_regulators_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DEG Regulators \n SH5") + 
  xlab("Regulons affected due to del_spoVG in SH5") + ylab("No. of genes in each regulon") +
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
regulators_SH5
coord_fixed()
dev.off() 
#Create a pie chart 
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_pie_regulators_SH5.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_SH5 <- regulators_SH5 + coord_polar("y", start=0)
pie_SH5
coord_fixed()
dev.off() 

#SH5_vs_SH2
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_barplot_regulators_SH5_vs_SH2.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DEG_SH5_vs_SH2_regulators_mod = DEG_SH5_vs_SH2_regulators[(DEG_SH5_vs_SH2_regulators$frequency)>20,]
regulators_SH5_vs_SH2 = ggplot(DEG_SH5_vs_SH2_regulators_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DEG Regulators \n SH5_vs_SH2") + 
  xlab("Regulons affected due to del_spoVG in SH5_vs_SH2") + ylab("No. of genes in each regulon") +
  theme(
    plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold")
  )
regulators_SH5_vs_SH2
coord_fixed()
dev.off() 
#Create a pie chart 
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_pie_regulators_SH5_vs_SH2.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_SH5_vs_SH2 <- regulators_SH5_vs_SH2 + coord_polar("y", start=0)
pie_SH5_vs_SH2
coord_fixed()
dev.off() 

#barplot_all
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/barplot/DEG_barplot_regulators_all.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
barplot_all = grid.arrange(regulators_M9, regulators_SH2,
                           regulators_SH5, regulators_SH5_vs_SH2, nrow = 4)
coord_fixed()
dev.off() 

#pie_all
pdf("./diff.exp.gene/DEG_SM/Figures/geneRegulations/regulators/pie/DEG_pie_regulators_all.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_all = grid.arrange(pie_M9, pie_SH2,
                       pie_SH5, pie_SH5_vs_SH2, nrow = 4)
coord_fixed()
dev.off() 










