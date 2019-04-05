library(RCy3)
library(igraph)
# Take data from dataPreparation.R script

#dataToPlot = as.data.frame(finalBestHit)
constructGraph = function(dataToPlot, 
                          threshold = 0.001,
                          saveNodesAndEdges = F){
  
  # select the rows where at least one of the organisms reached the threshold - drop the others
  data = dataToPlot[apply(dataToPlot[,-c(1,2)], 1, function(x) !all(x < threshold)),-1] 
  rownamesData = data$Organism
  data = data[, -1]
  rownames(data) = rownamesData

  
  ## create vertices (nodes)
  numSamples = ncol(data)  #number of samples
  numOTU = nrow(data)        #number of identified (bacterial) taxa
  numNodes = numSamples+numOTU    #number of nodes
  
  nodes = c(colnames(data), rownames(data))   #vector of nodes = their label
  
  rownames(metadata) = metadata$ID
  # type of node- for node coloring in cytoscape - for kids do Banglades vs UVA /  Organisms - just do Organism
  typeNode = c(as.character(metadata[colnames(data), "site"]), rep("Organism",numOTU))  
  type = c(rep(0, numSamples), rep(1,numOTU))  
  
  # get the medians / means / max across rows - median of the viral abundance - do log2 transformation and shift 
  # up to eliminate negative numbers -> the lowest abundance = 1
  medSize = log2(apply(data, 1, median) * 100)
  medSize = ceiling(medSize + abs(min(medSize)) + 10)
  
  meanSize = log2(apply(data, 1, mean) * 100)
  meanSize = ceiling(meanSize + abs(min(meanSize)) + 10)
  
  maxSize = log2(apply(data, 1, max) * 100)
  maxSize = ceiling(maxSize + abs(min(maxSize)) + 10)
  
  # set the kid nodes as the biggest and all the same size
  typeMedian = c(rep(max(medSize + 5), numSamples), medSize)  
  typeMean = c(rep(max(meanSize + 5), numSamples), meanSize)
  typeMax = c(rep(max(maxSize + 5), numSamples), maxSize)
  
  tableNodes = data.frame(nodes, type, typeNode, typeMedian, typeMean, typeMax)
  names(tableNodes) = c("Nodes", "type","typeNode", "sizeLog2Median", "sizeLog2Mean", "sizelog2Max")
  
  if(saveNodesAndEdges){
    write.table(tableNodes,"nodes.txt", quote = FALSE, row.names = FALSE)
    write.csv(tableNodes,"nodes.csv", quote = FALSE, row.names = FALSE)
  }

  ## create edges
  edgeSource = c()
  edgeTarget = c()
  edgeWeight = c()
  for(i in 1:numOTU){                     #for every taxon
    for(j in 1:numSamples){               #for every sample
      if(data[i,j]>threshold){         #if abundance of i-th taxon in j-th sample is higher that threshold
        edgeSource = c(edgeSource, i+numSamples)   #source node (taxon - virus)
        edgeTarget = c(edgeTarget, j)               #target node (sample)
        edgeWeight = c(edgeWeight, ceiling(5*(data[i,j]/max(data[i,]))))    #edge weighting
      }
    }
  }
  
  tableEdges = data.frame(edgeSource, edgeTarget, edgeWeight)
  names(tableEdges) = c("Source","Target","Weight") 
  
  if(saveNodesAndEdges){
    write.table(tableEdges,"graph/edges.txt", quote = FALSE, row.names = FALSE)
    write.csv(tableEdges,"graph/edges.csv", quote = FALSE, row.names = FALSE)
  }
  
  
  graphList = list(tableNodes = tableNodes, tableEdges=tableEdges)
  return(graphList)
}



# run the function
graphList = constructGraph(as.data.frame(finalBestHit))


# transport to cytoscape
tableNodes = graphList$tableNodes
tableEdges = graphList$tableEdges

cytNodes = tableNodes
cytEdges = tableEdges

# name the nodes
cytEdges$Source = tableNodes$Nodes[tableEdges$Source]  
cytEdges$Target = tableNodes$Nodes[tableEdges$Target]  
#rename edgets
colnames(cytEdges) = c("from", "to", "weight")

#asign te metadata to edges - extract the site - to color the edges
cytEdges$to = as.character(cytEdges$to)
colorVec = left_join(cytEdges, metadata, by = c("to" = "ID"))
cytEdges$site = colorVec$site


# create graph and export to cytoscape
ig = graph_from_data_frame(cytEdges, directed=F, vertices= cytNodes)
createNetworkFromIgraph(ig,"finalBestHit0_001") 




