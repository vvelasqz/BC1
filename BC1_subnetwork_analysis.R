
arab_interactome<- read.csv("Network analysis/arab_interactome.csv")
BC1_interactors<- read.csv("Network analysis/BC1_interactorsV2.csv")
sum(duplicated(BC1_interactors$Arabidopsis.homologue))
#one duplicated
bc1_interactions<- data.frame(Source= "BC1", Target=BC1_interactors$Arabidopsis.homologue)
nrow(arab_interactome[arab_interactome$Source %in% BC1_interactors$Arabidopsis.homologue | arab_interactome$Target %in% BC1_interactors$Arabidopsis.homologue,])
arab_interactome_bc1<- rbind(arab_interactome, bc1_interactions)
bc1_net<- arab_interactome_bc1[arab_interactome_bc1$Source %in% BC1_interactors$Arabidopsis.homologue | arab_interactome_bc1$Target %in% BC1_interactors$Arabidopsis.homologue,]
write.csv(bc1_net, file = "Network analysis/BC1_net_firstNeighV2.csv")

length(unique(c(bc1_net$Source, bc1_net$Target)))

#calculate node stats
calculate_net_stats<- function(weighted_network, annotation){
  arab_net0<- weighted_network
  V(arab_net0)$degree <- degree(graph = arab_net0)
  V(arab_net0)$closeness <- closeness(graph = arab_net0)
  V(arab_net0)$betweenness <- betweenness(graph = arab_net0)
  
  graphRankingStats_arab_net0 <- data.frame(
    degree = V(arab_net0)$degree,
    betweenness = V(arab_net0)$betweenness,
    closeness = V(arab_net0)$closeness,
    row.names = V(arab_net0)$name)
  
  graphRankingStats_arab_net0$node<- rownames(graphRankingStats_arab_net0)
  graphRankingStats_arab_net0<- merge(graphRankingStats_arab_net0, annotation, by.x = "node", by.y="gene", all.x = T)
  return(graphRankingStats_arab_net0)
}

library(igraph)
bc1_net0 <- graph_from_data_frame(bc1_net, directed = FALSE, vertices = NULL)
bc1_net0<-igraph::simplify(bc1_net0)

annotation<- read.csv("Network analysis/arab_annotation.csv", stringsAsFactors = F, quote = "\"")[,c(1,3)]
colnames(annotation)<-c("gene","description")
net_stats_arab<- calculate_net_stats(arab_net0, annotation)
net_stats_arab_bc1<- calculate_net_stats(arab_net1, annotation)
net_stats_bc1<- calculate_net_stats(bc1_net0, annotation)
Nodes_bc1_net0<- data.frame(node= unique(c(bc1_net$Source, bc1_net$Target)))
Nodes_bc1_net0<- merge(Nodes_bc1_net0, net_stats_arab_bc1, by="node", all.x = T)
Nodes_bc1_net0<- merge(Nodes_bc1_net0, BC1_interactors, by.x="node", by.y = "Arabidopsis.homologue", all.x = T)
write.csv(Nodes_bc1_net0, file = "Network analysis/BC1_nodes_firstNeighV2.csv")

write.csv(net_stats_arab_bc1, "Network analysis/Node_stats_arabInt_bc1V2.csv")
write.csv(net_stats_bc1, "Network analysis/Node_stats_bc1_net0V2.csv")

#####################################################################
###################################################################################################################

library(RCy3)
cytoscapePing()
#network with only the first neighbors

Edges_bc1_net0<- read.csv(file = "Network analysis/BC1_net_firstNeighV2.csv", row.names = 1)
colnames(Edges_bc1_net0)<- c("source", "target")
Nodes_bc1_net0<- read.csv(file = "Network analysis/BC1_nodes_firstNeighV2.csv", row.names = 1)
hist(Nodes_bc1_net0[!is.na(Nodes_bc1_net0$Type), "degree"], breaks = length(unique(Nodes_bc1_net0[!is.na(Nodes_bc1_net0$Type), "degree"])))
names(Nodes_bc1_net0)[1]<- "id"
Nodes_bc1_net0[Nodes_bc1_net0$id=="BC1", "Type"]<- "BC1"
Nodes_bc1_net0[is.na(Nodes_bc1_net0$Type), "Type"]<- "AtInt"

Nodes_bc1_net0$log_degree<- 10*log10(Nodes_bc1_net0$degree)+10

style.name = "BC1V2"
defaults <- list(NODE_SHAPE="ellipse",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_FILL_COLOR="#AAAAAA",
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','Type','d', c("Co_IP", "Y2H", "Co_IP;Y2H","BC1", "AtInt"),
                               c("#990000","#3333FF", "#800080","#00994C", "#AAAAAA"))
nodeSize <- mapVisualProperty('node size','log_degree','p')
createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,nodeSize))

createNetworkFromDataFrames(Nodes_bc1_net0,Edges_bc1_net0, title="BC1 first neighbors V2", collection="DataFrame Example")
setVisualStyle(style.name)

#in cytoscape layout group attribute layout, cluster





