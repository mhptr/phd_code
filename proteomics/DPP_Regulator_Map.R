
#Load library and set working dir

library(igraph)

setwd("~/OneDrive - University of Warwick/WORK/RESULTS/PROTEOMICS/FINAL Result")


diffexp.pr.all.LB = read.csv("./Analysis - R-Script /Data/combined/output/diff.exp./diff.exp._minus_spoVG_upp/20190514_diffexp.pr.all_LB_minus_spoVG_upp.csv", header = T)
rownames(diffexp.pr.all.LB) = diffexp.pr.all.LB$X
diffexp.pr.all.LB = diffexp.pr.all.LB[,-1]

################################################
#### REGULATORS Interaction Maps
################################################

regulations = read.csv("./Analysis - R-Script /SubtiWiki Exports /regulations.csv")
regulations2 = regulations[!as.character(regulations$regulator.locus.tag)==as.character(regulations$locus.tag),]  #removes auotoregulat


###############################
##REGULATORS - Universe
###############################

regulations2.U = regulations2[regulations2$regulator.locus.tag%in%rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05] | 
                                regulations2$locus.tag%in%rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05], ]


#regulations2.U = regulations2 
g.U = graph.data.frame(cbind(as.character(regulations2.U$regulator.locus.tag),
                             as.character(regulations2.U$locus.tag)), directed=F)  #no direction decided at this step
g.U = igraph::simplify(g.U)  #sort of finds the unique interactions
nodes.U  = rownames(as.matrix(V(g.U))) #V is the vertex which is the node
edge.list.U = as.data.frame(get.edgelist(g.U))
degree.U = as.matrix(igraph::degree(g.U),ncol=1) #prints the number of edges per vertex
degree.U = as.data.frame(degree.U)
colnames(degree.U) = "Degree"

#plot(g.U,vertex.size=3,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0,
#     mark.border="grey",
#     #vertex.color=dim(degree.U),
#     pch=19,vertex.label=NA)    #check pch 


##########################
##REGULATORS - LB - Non-directional
##########################


regulations2.LB = regulations2[regulations2$regulator.locus.tag%in%rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05] | 
                                 regulations2$locus.tag%in%rownames(diffexp.pr.all.LB)[diffexp.pr.all.LB$adj.P.Val<0.05], ]

g.LB = graph.data.frame(cbind(as.character(regulations2.LB$regulator.locus.tag),
                              as.character(regulations2.LB$locus.tag)), directed=T)  #no direction decided at this step

#add gene names with the bsu numbers
Gene = read.csv("./Analysis - R-Script /SubtiWiki Exports /Gene.csv", header =T)
rownames(Gene) = rownames(Gene$locus)

g.LB = igraph::simplify(g.LB)  #sort of finds the unique interactions

nodes.LB  = rownames(as.matrix(V(g.LB))) #V is the vertex which is the node
edge.list.LB = as.data.frame(get.edgelist(g.LB))
degree.LB = as.matrix(igraph::degree(g.LB),ncol=1) #print the number of edges per vertex
degree.LB = as.data.frame(degree.LB)
colnames(degree.LB) = "Degree"


#Self exploratory 17/Sep/19
plot(g.LB,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0.4, 
     mark.border="red",pch=19,vertex.label=degree.LB$Degree, #rownames(degree.LB[(degree.LB$Degree)>10 ,]) = gives BSU numbers in the name 
     xlim=c(-0.8,0.8), ylim = c(-0.5,1.5),
     #vertex.color=colorRampPalette(c("white", "red")) (10),
     vertex.color=degree.U$Degree,        #gives multiple colours 
     #vertex.color=dim(degree.U),  #gives two colours : orange and blue
     vertex.label.color="black", vertex.label.dist=1.5,
     #vertex.shape=diffexp.pr.all.LB_SIG$P.Value,
     vertex.shape="circle",  #how to set two shapes for two kinds of expression data
     vertex.size=degree.LB$Degree, scale=degree.LB$Degree)    #check pch #[degree.LB$Degree<300 = 800,]
#how to make the size of each each in a certain margin length eg all value >20 =20
title("Regulatory Networks in LB - Proteomics",cex.main=1,col.main="black")


























#Self exploratory 18/Jan/19

#degree.LB = degree.LB[!apply(degree.LB, function(x) is.na(degree.LB), ...)]
#degree.LB[(as.numeric(degree.LB$Degree)<500), ] = "blue"
#degree.LB[(as.numeric(degree.LB$Degree)>500), ] = "red"
#kegg.enrich.sig.pv.log10[kegg.enrich.sig.pv.log10 < -log10(0.05)] = 0
#kegg.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10)


#Self exploratory 17/Jan/19
plot(g.LB,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0.4, 
     mark.border="red",pch=19,vertex.label=degree.LB$Degree, #rownames(degree.LB[(degree.LB$Degree)>10 ,]) = gives BSU numbers in the name 
     xlim=c(-0.8,0.8), ylim = c(-0.5,1.5),
     #vertex.color=colorRampPalette(c("white", "red")) (10),
     vertex.color=degree.U$Degree,        #gives multiple colours 
     #vertex.color=dim(degree.U),  #gives two colours : orange and blue
     vertex.label.color="black", vertex.label.dist=1.5,
     #vertex.shape=diffexp.pr.all.LB_SIG$P.Value,
     vertex.shape="circle",  #how to set two shapes for two kinds of expression data
     vertex.size=degree.LB$Degree, scale=degree.LB$Degree)    #check pch #[degree.LB$Degree<300 = 800,]
#how to make the size of each each in a certain margin length eg all value >20 =20
title("Regulatory Networks in LB - Proteomics",cex.main=1,col.main="black")

#biggest dot is BSU25200(30) = sigA row 12
#second biggest  is BSU25490 (9) = hrcA row 11
#third BSU00370 (8) = abrB row 47



#Done with Tauqeer
# plot(g.LB,vertex.size=3,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0,
#       mark.border="grey",pch=19,vertex.label=NA)    #check pch
