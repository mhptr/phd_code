#making graphs and charts of DPP
#18/08/19

setwd("~/OneDrive - University of Warwick/WORK/Results/Proteomics/FINAL Result/Analysis - R-Script /Data/")

DPP_LB_regulations = read.csv("./combined/output/geneRegulations/old/regulations/DPP_LB_geneRegulations.csv", header = T)
DPP_M9_regulations = read.csv("./combined/output/geneRegulations/old/regulations/DPP_M9_geneRegulations.csv", header = T)
DPP_SH2_regulations = read.csv("./combined/output/geneRegulations/old/regulations/DPP_SH2_geneRegulations.csv", header = T)
DPP_SH5_regulations = read.csv("./combined/output/geneRegulations/old/regulations/DPP_SH5_geneRegulations.csv", header = T)
DPP_SH5_vs_SH2_regulations = read.csv("./combined/output/geneRegulations/old/regulations/DPP_SH5_vs_SH2_geneRegulations.csv", header = T)

# dim(DPP_LB_regulations)
# DPP_LB_regulations = DPP_LB_regulations[!DPP_LB_regulations$regulator == "",]
# dim(DPP_LB_regulations)
# write.csv(DPP_LB_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_LB_regulations.csv", row.names = F)

# dim(DPP_M9_regulations)
# DPP_M9_regulations = DPP_M9_regulations[!DPP_M9_regulations$regulator == "",]
# dim(DPP_M9_regulations)
# write.csv(DPP_M9_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_M9_regulations.csv", row.names = F)

# dim(DPP_SH2_regulations)
# DPP_SH2_regulations = DPP_SH2_regulations[!DPP_SH2_regulations$regulator == "",]
# dim(DPP_SH2_regulations)
# write.csv(DPP_SH2_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_SH2_regulations.csv", row.names = F)

# dim(DPP_SH5_regulations)
# DPP_SH5_regulations = DPP_SH5_regulations[!DPP_SH5_regulations$regulator == "",]
# dim(DPP_SH5_regulations)
# write.csv(DPP_SH5_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_SH5_regulations.csv", row.names = F)
# 
# dim(DPP_SH5_vs_SH2_regulations)
# DPP_SH5_vs_SH2_regulations = DPP_SH5_vs_SH2_regulations[!DPP_SH5_vs_SH2_regulations$regulator == "",]
# dim(DPP_SH5_vs_SH2_regulations)
# write.csv(DPP_SH5_vs_SH2_regulations, "./combined/output/geneRegulations/old/regulations/csv/DPP_SH5_cs_SH2_regulations.csv", row.names = F)


DPP_LB_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_LB_regulations.csv", header = T)
DPP_M9_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_M9_regulations.csv", header = T)
DPP_SH2_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_SH2_regulations.csv", header = T)
DPP_SH5_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_SH5_regulations.csv", header = T)
DPP_SH5_vs_SH2_regulations = read.csv("./combined/output/geneRegulations/old/regulations/csv/DPP_SH5_cs_SH2_regulations.csv", header = T)

# DPP_all = list(DPP_LB_regulations, DPP_M9_regulations, DPP_SH2_regulations, DPP_SH5_regulations, DPP_SH5_vs_SH2_regulations)
# drop_empty_row = for (i in DPP_all) {
#   DPP_regulations = DPP_all(i)[!DPP_all(i)$regulator == "",]
# }

DPP_LB_regulators = table(DPP_LB_regulations$regulator)
DPP_M9_regulators = table(DPP_M9_regulations$regulator)
DPP_SH2_regulators = table(DPP_SH2_regulations$regulator)
DPP_SH5_regulators = table(DPP_SH5_regulations$regulator)
DPP_SH5_vs_SH2_regulators = table(DPP_SH5_vs_SH2_regulations$regulator)

# write.csv(DPP_LB_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_LB_regulators.csv", row.names = F)
# DPP_LB_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_LB_regulators.csv", header = T)
# colnames(DPP_LB_regulators) = c("regulator", "frequency")
# dim(DPP_LB_regulators)
# DPP_LB_regulators = DPP_LB_regulators[!DPP_LB_regulators$regulator == "",]
# dim(DPP_LB_regulators)
# View(DPP_LB_regulators)
# write.csv(DPP_LB_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_LB_regulators.csv", row.names = F)

# DPP_M9_regulators = DPP_M9_regulators[!DPP_M9_regulators$Name == "",] #removes the empty row
# write.csv(DPP_M9_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", row.names = F)
# DPP_M9_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", header = T)
# colnames(DPP_M9_regulators) = c("regulator", "frequency")
# dim(DPP_M9_regulators)
# DPP_M9_regulators = DPP_M9_regulators[!DPP_M9_regulators$regulator == "",]
# dim(DPP_M9_regulators)
# View(DPP_M9_regulators)
# write.csv(DPP_M9_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", row.names = F)

# write.csv(DPP_SH2_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH2_regulators.csv", row.names = F)
# DPP_SH2_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH2_regulators.csv", header = T)
# colnames(DPP_SH2_regulators) = c("regulator", "frequency")
# dim(DPP_SH2_regulators)
# DPP_SH2_regulators = DPP_SH2_regulators[!DPP_SH2_regulators$regulator == "",]
# dim(DPP_SH2_regulators)
# View(DPP_SH2_regulators)
# write.csv(DPP_SH2_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH2_regulators.csv", row.names = F)

# write.csv(DPP_SH5_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_regulators.csv", row.names = F)
# DPP_SH5_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_regulators.csv", header = T)
# colnames(DPP_SH5_regulators) = c("regulator", "frequency")
# dim(DPP_SH5_regulators)
# DPP_SH5_regulators = DPP_SH5_regulators[!DPP_SH5_regulators$regulator == "",]
# dim(DPP_SH5_regulators)
# View(DPP_SH5_regulators)
# write.csv(DPP_SH5_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_regulators.csv", row.names = F)

# write.csv(DPP_SH5_vs_SH2_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_vs_SH2_regulators.csv", row.names = F)
# DPP_SH5_vs_SH2_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_vs_SH2_regulators.csv", header = T)
# colnames(DPP_SH5_vs_SH2_regulators) = c("regulator", "frequency")
# dim(DPP_SH5_vs_SH2_regulators)
# DPP_SH5_vs_SH2_regulators = DPP_SH5_vs_SH2_regulators[!DPP_SH5_vs_SH2_regulators$regulator == "",]
# dim(DPP_SH5_vs_SH2_regulators)
# View(DPP_SH5_vs_SH2_regulators)
# write.csv(DPP_SH5_vs_SH2_regulators, "./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_vs_SH2_regulators.csv", row.names = F)
# 

DPP_LB_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_LB_regulators.csv", header = T)
DPP_M9_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_M9_regulators.csv", header = T)
DPP_SH2_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH2_regulators.csv", header = T)
DPP_SH5_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_regulators.csv", header = T)
DPP_SH5_vs_SH2_regulators = read.csv("./combined/output/geneRegulations/old/regulators/csv/DPP_SH5_vs_SH2_regulators.csv", header = T)


library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(scales)

swati.color = c(brewer.pal(8,"Dark2"))

#http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization

# Barplot

#LB
DPP_LB_regulators_mod = DPP_LB_regulators[(DPP_LB_regulators$frequency)>2,]
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/barplot/barplot_regulators_LB.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
regulators_LB = ggplot(DPP_LB_regulators_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DPP Regulators \n LB") + 
  xlab("Regulons affected due to del_spoVG in LB") + ylab("No. of genes in each regulon") +
  theme(
       plot.title = element_text(hjust = 0.5, color="darkBlue", face="bold", size=14),
       axis.title.x = element_text(color="black", size=12, face="bold"),
       axis.title.y = element_text(color="black", size=12, face="bold")
     )
regulators_LB
coord_fixed()
dev.off()  
#Create a pie chart 
# pie_LB = regulators_LB + coord_polar("y", start=0) + geom_text(aes(y = frequency/length(DPP_LB_regulators_mod$regulator) + c(0, cumsum(frequency)[-length(frequency)]), 
#                                                                    label = percent(frequency/100)), size=5)
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/pie/pie_regulators_LB.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_LB = regulators_LB + coord_polar("y", start=0)
pie_LB
coord_fixed()
dev.off() 
# # use brewer color palettes
# pie + scale_fill_brewer(palette="Blues") #need 54 

#M9
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/barplot/barplot_regulators_M9.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DPP_M9_regulators_mod = DPP_M9_regulators[(DPP_M9_regulators$frequency)>3,]
regulators_M9 = ggplot(DPP_M9_regulators_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DPP Regulators \n M9") + 
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
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/pie/pie_regulators_M9.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_M9 <- regulators_M9 + coord_polar("y", start=0)
pie_M9
coord_fixed()
dev.off() 

#SH2
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/barplot/barplot_regulators_SH2.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DPP_SH2_regulators_mod = DPP_SH2_regulators[(DPP_SH2_regulators$frequency)>20,]
regulators_SH2 = ggplot(DPP_SH2_regulators_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DPP Regulators \n SH2") + 
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
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/pie/pie_regulators_SH2.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_SH2 <- regulators_SH2 + coord_polar("y", start=0)
pie_SH2
coord_fixed()
dev.off() 


#SH5
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/barplot/barplot_regulators_SH5.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DPP_SH5_regulators_mod = DPP_SH5_regulators[(DPP_SH5_regulators$frequency)>10,]
regulators_SH5 = ggplot(DPP_SH5_regulators_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DPP Regulators \n SH5") + 
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
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/pie/pie_regulators_SH5.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_SH5 <- regulators_SH5 + coord_polar("y", start=0)
pie_SH5
coord_fixed()
dev.off() 

#SH5_vs_SH2
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/barplot/barplot_regulators_SH5_vs_SH2.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
DPP_SH5_vs_SH2_regulators_mod = DPP_SH5_vs_SH2_regulators[(DPP_SH5_vs_SH2_regulators$frequency)>20,]
regulators_SH5_vs_SH2 = ggplot(DPP_SH5_vs_SH2_regulators_mod, aes(x="", y=frequency, fill=regulator))+
  geom_bar(width = 2, stat = "identity") + ggtitle("DPP Regulators \n SH5_vs_SH2") + 
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
pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/pie/pie_regulators_SH5_vs_SH2.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_SH5_vs_SH2 <- regulators_SH5_vs_SH2 + coord_polar("y", start=0)
pie_SH5_vs_SH2
coord_fixed()
dev.off() 


pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/barplot/barplot_regulators_all.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
barplot_all = grid.arrange(regulators_LB , regulators_M9, regulators_SH2,
                           regulators_SH5, regulators_SH5_vs_SH2, nrow = 4)
coord_fixed()
dev.off() 

pdf("../Figures/combined/20190110_combined/adj.P.Val/Final/geneRulations/regulators/pie/pie_regulators_all.pdf" ,width=10, height=10)      #turn this OFF if just want to see the picture in the Plots
pie_all = grid.arrange(pie_LB , pie_M9, pie_SH2,
             pie_SH5, pie_SH5_vs_SH2, nrow = 3)
coord_fixed()
dev.off() 
















hist(DPP_LB_regulators,
     main="Regulons affected due to del_spoVG in LB",
     xlab="Regulators",
     xlim=c(50,100),
     col="No. of genes in each regulon",
     freq=FALSE
)


plot(DPP_LB_regulators)

hist(DPP_LB_regulators$frequency)

# Add a Normal Curve (Thanks to Peter Dalgaard)
x <- DPP_LB_regulators$frequency
h<-hist(x, breaks=10, col="pink", xlab="Regulons affected due to del_spoVG in LB", 
        main="Histogram with Normal Curve", 
        ylab = "No. of genes in each regulon") 
xfit<-seq(min(x),max(x),length=100) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="red", lwd=2)

# Add all the labels and color. Set the y axis indexes horizontal. Set limits for axis and
# Drawing lables on top of bars
hist(DPP_LB_regulators$frequency,main="Regulons affected due to del_spoVG in LB", xlab = "Regulons affected due to del_spoVG in LB", ylab = "No. of genes in each regulon",border="red", col="blue",las=1,xlim=c(1,40),ylim=c(0,60),labels=TRUE)
hist(DPP_SH5_regulators$frequency,main="Regulons affected due to del_spoVG in LB", xlab = "Regulons affected due to del_spoVG in LB", ylab = "No. of genes in each regulon",border="red", col="blue",las=1,xlim=c(1,100),ylim=c(0,100),labels=TRUE)




# #http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization
# library(ggplot2) 
# library(gridExtra)
# p1 <- ggplot(DPP_LB_regulators, aes(x = frequency)) +
#   geom_histogram(binwidth = 1000, color = "grey30", fill = "white") +
#   ggtitle("Bin Width = 1000")
# p2 <- ggplot(DPP_M9_regulators, aes(x = Freq)) +
#   geom_histogram(binwidth = 5000, color = "grey30", fill = "white") +
#   ggtitle("Bin Width = 5000")
# p3 <- ggplot(DPP_SH2_regulators, aes(x = frequency)) +
#   geom_histogram(binwidth = 10000, color = "grey30", fill = "white") +
#   ggtitle("Bin Width = 10000")
# grid.arrange(p1, p2, p3, ncol=3)
# library(ggplot2)
# # Basic histogram
# ggplot(DPP_LB_regulators, aes(x=frequency)) + geom_histogram()
# # Change the width of bins
# ggplot(DPP_LB_regulators, aes(x=frequency)) + 
#   geom_histogram(binwidth=1)
# # Change colors
# p<-ggplot(DPP_LB_regulators, aes(x=frequency)) + 
#   geom_histogram(color="black", fill="white")
# p
# library(plyr)
# mu <- ddply(DPP_LB_regulators, "regulator", summarise, grp.mean=mean(frequency))
# head(mu)
# # Change histogram plot line colors by groups
# ggplot(DPP_M9_regulators, aes(x=Freq, color=regulator)) +
#   geom_histogram(fill="white")
# # Overlaid histograms
# ggplot(DPP_LB_regulators, aes(x=frequency, color=regulator)) +
#   geom_histogram(fill="white", alpha=0.5, position="identity")
# DPP_M9_regulators = DPP_M9_regulators[(DPP_M9_regulators$Freq)>=3,]
# dim(DPP_M9_regulators)
# rownames(DPP_M9_regulators) = FALSE
# hist(DPP_M9_regulators$Var1)
# ggplot(DPP_M9_regulators, aes(x=Freq, color=Var1)) +
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














