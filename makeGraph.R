plotTable = tableToPlot %>% 
  filter(Genome %in% filteredGenomes) 

data = plotTable[,-1]

graph_constructor=function(data, threshold=0.002){
  ## create vertices (nodes)
  num_samples <- ncol(data)-1  #number of samples
  num_OTU <- nrow(data)        #number of identified (bacterial) taxa
  num_nodes <- num_samples+num_OTU    #number of nodes
  
  nodes <- c(names(data)[2:length(names(data))], as.character(data[,1]))   #vector of nodes = their label
  id <- c(1:num_nodes)          #id of nodes (just numbers)
  type <- c(rep(1,num_samples), rep(2,num_OTU))   # type of node, you can change size of nodes in Gephi with this (larger samples, smaller OTUs)
  
  table_nodes <- data.frame(nodes,id,nodes,type)
  names(table_nodes) <- c("Nodes","Id","Label","Type")
  
  write.table(table_nodes,"graph/nodes.txt", quote = FALSE, row.names = FALSE)
  write.csv(table_nodes,"graph/nodes.csv", quote = FALSE, row.names = FALSE)
  ## create edges
  data_num <- data[,1:num_samples+1]     #only numbers from OTU table
  ind <- 1
  edge_source <- c()
  edge_target <- c()
  edge_id <- c()
  edge_weight <- c()
  for(i in 1:num_OTU){                     #for every taxon
    for(j in 1:num_samples){               #for every sample
      if(data_num[i,j]>threshold){         #if abundance of i-th taxon in j-th sample is higher that threshold
        edge_source[ind] <- i+num_samples   #source node (taxon)
        edge_target[ind] <- j               #target node (sample)
        edge_id[ind] <- ind                 #edge id (just number)
        edge_weight[ind] <- ceiling(10*(data_num[i,j]/max(data_num[i,])))    #edge weighting
        ind <- ind+1
      }
    }
  }
  edge_type <- rep("undirected", (ind-1))   #all edges are undirected
  edge_label <- rep("xxx", (ind-1))         #there is no need for actual labeling of edges but gephi probably need this
  
  table_edges <- data.frame(edge_source, edge_target, edge_type, edge_id, edge_label, edge_weight)
  names(table_edges) <- c("Source","Target","Type","Id","Label","Weight") 
  
  write.table(table_edges,"graph/edges.txt", quote = FALSE, row.names = FALSE)
  write.csv(table_edges,"graph/edges.csv", quote = FALSE, row.names = FALSE)
  
}



# recreate the table
cytNodes = table_nodes[,-3]

cytEdges = table_edges
cytEdges$Source = table_nodes$Nodes[table_edges$Source]  
cytEdges$Target = table_nodes$Nodes[table_edges$Target]  
cytEdges = cytEdges[,-c(3,4,5)]
colnames(cytEdges) = c("from", "to", "weight")

cytNodes = cytNodes[,-2]


ig <- graph_from_data_frame(cytEdges, directed=F, vertices= cytNodes)
createNetworkFromIgraph(ig,"myIgraph") 

