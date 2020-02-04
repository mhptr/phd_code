####################################################################
### 1. Include Packages 
####################################################################

library(RColorBrewer)
library(igraph)

####################################################################
### 2. Set Working Directory 
####################################################################

setwd("~/OneDrive - University of Warwick/WORK/RESULTS/Transcriptomic/intersection_toRNAdo_all/filt_txt_all_2/")

################################################
######REGULATORS Interaction Maps
################################################

library(igraph)
regulations = read.csv("../SubtiWiki Exports /regulations.csv")
#this removes autoregulations
regulations2 = regulations[!as.character(regulations$regulator.locus.tag)==as.character(regulations$locus.tag),]  #removes auotoregulat

###############################
##REGULATORS - Universe
###############################

DEG.M9 = read.csv("../diff.exp.gene/DEG_SM/Diff_all_sm/Diff_all_sm_minus_spoVG_upp/DEG_ko_vs_wt_all_M9_minus_spoVG_upp.csv", header = T)
rownames(DEG.M9) = DEG.M9$X
DEG.M9 = DEG.M9[,-1]

# either regulator or the target is differentially expressed
regulations2.U = regulations2[regulations2$regulator.locus.tag%in%rownames(DEG.M9)[DEG.M9$padj<0.05] | 
                                 regulations2$locus.tag%in%rownames(DEG.M9)[DEG.M9$padj<0.05], ]
#regulations2.U = regulations2 
g.U = graph.data.frame(cbind(as.character(regulations2.U$regulator.locus.tag),
                              as.character(regulations2.U$locus.tag)), directed=F)  #no direction decided at this step
g.U = igraph::simplify(g.U)  #sort of finds the unique interactions
nodes.U  = rownames(as.matrix(V(g.U))) #V is the vertex which is the node
edge.list.U = as.data.frame(get.edgelist(g.U))
degree.U = as.matrix(igraph::degree(g.U),ncol=1) #print the number of edges per vertex
degree.U = as.data.frame(degree.U)
colnames(degree.U) = "Degree"
#plot(g.U,vertex.size=3,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0,
#     mark.border="grey",
#     #vertex.color=dim(degree.U),
#     pch=19,vertex.label=NA)    #check pch 


##############################################################################
##REGULATORS - M9 - Non-directional
##############################################################################

regulations2.M9 = regulations2[regulations2$regulator.locus.tag%in%rownames(DEG.M9)[DEG.M9$padj<0.05] | 
  regulations2$locus.tag%in%rownames(DEG.M9)[DEG.M9$padj<0.05], ]
g.M9 = graph.data.frame(cbind(as.character(regulations2.M9$regulator.locus.tag),
                            as.character(regulations2.M9$locus.tag)), directed=T)  #no direction decided at this step
#add gene names with the bsu numbers
Gene = read.csv("../SubtiWiki Exports /Gene.csv", header =T)
rownames(Gene) = rownames(Gene$locus)
g.M9 = igraph::simplify(g.M9)  #sort of finds the unique interactions
nodes.M9  = rownames(as.matrix(V(g.M9))) #V is the vertex which is the node
edge.list.M9 = as.data.frame(get.edgelist(g.M9))
degree.M9 = as.matrix(igraph::degree(g.M9),ncol=1) #print the number of edges per vertex
degree.M9 = as.data.frame(degree.M9)
colnames(degree.M9) = "Degree"

#Self exploratory 18/Jan/19
#degree.LB = degree.LB[!apply(degree.LB, function(x) is.na(degree.LB), ...)]
#degree.LB[(as.numeric(degree.LB$Degree)<500), ] = "blue"
#degree.LB[(as.numeric(degree.LB$Degree)>500), ] = "red"
#kegg.enrich.sig.pv.log10[kegg.enrich.sig.pv.log10 < -log10(0.05)] = 0
#kegg.enrich.sig.pv.log10,color = colorRampPalette(c("white", "blue")) (10)

#Self exploratory 17/Jan/19
# plot(g.LB,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0.4, 
#      mark.border="red",pch=19,vertex.label=degree.LB$Degree, #rownames(degree.LB[(degree.LB$Degree)>10 ,]) = gives BSU numbers in the name 
#      xlim=c(-0.8,0.8), ylim = c(-0.5,1.5),
#      #vertex.color=colorRampPalette(c("white", "red")) (10),
#      vertex.color=degree.U$Degree,        #gives multiple colours 
#      #vertex.color=dim(degree.U),  #gives two colours : orange and blue
#      vertex.label.color="black", vertex.label.dist=1.5,
#      #vertex.shape=DEG.LB_SIG$P.Value,
#      vertex.shape="circle",  #how to set two shapes for two kinds of expression data
#      vertex.size=degree.LB$Degree, scale=degree.LB$Degree)    #check pch #[degree.LB$Degree<300 = 800,]
# #how to make the size of each each in a certain margin length eg all value >20 =20
# title("Regulatory Networks in LB - Proteomics",cex.main=1,col.main="black")
#biggest dot is BSU25200(30) = sigA row 12
#second biggest  is BSU25490 (9) = hrcA row 11
#third BSU00370 (8) = abrB row 47


###########################################
#Done with Tauqeer
# plot(g.LB,vertex.size=3,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0,
#       mark.border="grey",pch=19,vertex.label=NA)    #check pch
###########################################

#Self exploratory 30/Sep/19
l <- layout_with_fr(g.M9)
plot(g.M9, layout = l, vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0.2, 
     mark.border="red",pch=19,vertex.label=degree.M9$Degree, #rownames(degree.M9[(degree.M9$Degree)>10 ,]) = gives BSU numbers in the name 
     xlim=c(-0.8,0.8), ylim = c(-0.5,1.5),
     #vertex.color=colorRampPalette(c("white", "red")) (10),
     vertex.color=degree.U$Degree,        #gives multiple colours 
     #vertex.color=dim(degree.U),  #gives two colours : orange and blue
     vertex.label.color="black", vertex.label.dist=1.5,
     #vertex.shape=DEG.M9_SIG$P.Value,
     vertex.shape="circle",  #how to set two shapes for two kinds of expression data
     vertex.size=degree.M9$Degree, 
     scale=degree.M9$Degree)    #check pch #[degree.M9$Degree<300 = 800,]
#how to make the size of each each in a certain margin length eg all value >20 =20
#title("Regulatory Networks in M9 - Transcriptomic",cex.main=1,col.main="black")


#simple graph
plot(g.M9,vertex.size=3,vertex.frame.color=NA,vertex.label.dist=0.5, edge.arrow.size=0,
          mark.border="grey",pch=19,vertex.label=NA)






























#28/0919 
#John A's script

# data1=read.csv("all_vsall.csv")
# CNV_data=read.csv(args[2])
# library("igraph") 
# g  <- graph.adjacency(as.matrix(data1), weighted=T,mode="undirected",diag=F)
# clp2=cluster_fast_greedy(g)
# pol= fastgreedy.community(g)
# l <- layout_with_fr(g)
# plot(g,layout=l,vertex.size=10,curved=0,edge.width=1.5,vertex.label=NA)


regulations = read.csv("./Analysis - R-Script /SubtiWiki Exports /regulations.csv")
regulations2 = regulations[!as.character(regulations$regulator.locus.tag)==as.character(regulations$locus.tag),]  #removes auotoregulat

regulations2.U = regulations2[regulations2$regulator.locus.tag%in%rownames(DEG.LB)[DEG.LB$padj<0.05] | 
                                regulations2$locus.tag%in%rownames(DEG.LB)[DEG.LB$padj<0.05], ]

DEG.LB = read.csv("./Analysis - R-Script /Data/combined/output/diff.exp./diff.exp._minus_spoVG_upp/20190514_DEG_LB_minus_spoVG_upp.csv", header = T)
rownames(DEG.LB) = DEG.LB$X
DEG.LB = DEG.LB[,-1]

# regulations2 = regulations2[,-c(1,3,4,6)]
# data1 = get.adjacency(graph.edgelist(as.matrix(regulations2), directed=FALSE))
# CNV_data=DEG.LB
# library("igraph") 
# g  <- graph.adjacency(as.matrix(data1), weighted=T,mode="undirected",diag=F)
# clp2=cluster_fast_greedy(g)
# pol= fastgreedy.community(g)
# # Make sure the nodes stay in place in both plot
# l <- layout_with_fr(g)
# #plot(g,layout=l,vertex.size=10,curved=0,edge.width=1.5,vertex.label=NA)
# 
# # # Generate colors based on media type
# # colrs <- c("gray50", "tomato", "gold")
# # V(g)$color <- colrs[V(g)$regulator]# 
# # # Set edge width based on weight:# 
# # E(net)$width <- E(net)$weight/6
# 
# 
# plot(g, main = "DPP : LB",  edge.width = 2, edge.arrow.size=.5, vertex.size=10, vertex.label.color="black", vertex.label.dist=1.5,
#      vertex.color=c( "pink", "skyblue", "red", "yellow"))























