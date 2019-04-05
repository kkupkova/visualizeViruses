rm(list=ls())
library(data.table)
library(ggplot2)
library(reshape2)

filePath = "results_with_taxa_names/"

# get the tax IDs in all of the TSV files - starts with virus, ends with .tsv
tsvFiles = list.files(path = filePath)

# get the ID numbers from the TSV files
ID = regmatches(tsvFiles, gregexpr("[[:digit:]]+", tsvFiles))
ID = unlist(ID)

condition = rep("Dhaka", length(ID))

# get the first digit in the ID number
ID1 = substr(ID, 1, 1)
condition[ID1 == "7"] = "UVA"

condition = as.factor(condition)

# upload the sample table, I make a dummy example here for these two IDs
sample = data.frame(ID = ID, area = condition)

rm(condition, ID, ID1)

# tsvFiles = the files to be uploaded
# viralNameColumn = the column which contains the information about the viruses,
#                   the column name must match the column name in the file exactly
# valueColumn = the name of the column containing information about the virus abundance
# threshold = if at least one of the samples crosses this value, the virus will be plotted

tsvFiles = list.files(path = filePath, full.names = T)
plotVirus = function(tsvFiles, viralNameColumn, valueColumn, sample, threshold){
  
  #make a list containng data.tables with Genome and selected column
  fileList = list()
  
  for (i in tsvFiles){
    
    #uploa tables one by one
    helpFile = fread(i, skip = 1, sep = "\t", header = T)
    
    # select desired columns
    myTable = helpFile[,mget(c(viralNameColumn,valueColumn))]
    
    setkeyv(myTable, viralNameColumn)
    
    # add to the list
    IDi = unlist(regmatches(i, gregexpr("[[:digit:]]+", i)))
    
    # make specific column for each sample, so they can be succesfully merged later
    colnames(myTable) = c(viralNameColumn, IDi)
    
    fileList[[IDi]] = myTable
    
    rm(myTable)
  }
  
  # merge the datasets and replace NAs with 0
  joinTable = plyr::join_all(fileList, by = viralNameColumn, type = "full")
  joinTable[is.na(joinTable)] = 0
  
  #select columns, where at least one sample hit the predefined threshold
  filteredTable = joinTable[apply(joinTable[,-1], 1, function(x) !all(x < threshold)),]

  # reshape the final table for ggplot
  plotTable = melt(filteredTable, id.vars = viralNameColumn)
  plotTable[, viralNameColumn] = as.factor(plotTable[, viralNameColumn])
  
  colnames(plotTable) = c(viralNameColumn, "ID", "value")
  
  plotTable = merge(plotTable, sample, by="ID")
  # add log of the value
  plotTable$logVal = log(plotTable$value)
  
  return(list(plotTable = plotTable, filteredTable = filteredTable))
  
}



viralNameColumn = "Organism"
valueColumn = "Final Best Hit"
# run the function:
tableList = plotVirus(tsvFiles = tsvFiles, 
                      viralNameColumn = viralNameColumn, 
                      valueColumn = valueColumn, 
                      sample = sample, 
                      threshold = 0.005)

plotTable = tableList$plotTable
filteredTable = tableList$filteredTable


# reorder viral labels according to the mean abundance of the in Dhaka kids
Dhaka = filteredTable[, c(viralNameColumn, sample$ID[sample$area == "Dhaka"])]
means = apply(Dhaka[,-1], 1, mean)
factorOrder = filteredTable[order(means, decreasing = T),1]

plotTable[,viralNameColumn] = factor(plotTable[,viralNameColumn], levels = factorOrder)

p = ggplot(plotTable, aes(x= get(viralNameColumn), y = ID, size = value, fill = area))
p + geom_point(shape = 21) + theme_bw() +
  scale_size_continuous(range=c(2,15)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(viralNameColumn) +
  ggtitle(valueColumn)


# plot logarithm of the percanetage
p = ggplot(plotTable, aes(x= get(viralNameColumn), y = ID, size = logVal, fill = area))
p + geom_point(shape = 21) + theme_bw() +
  scale_size_continuous(range=c(2,15)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(viralNameColumn) +
  ggtitle(paste0("Log of ", valueColumn))


# remove the most abundant column to renormalize the data
plotTable_nonAlpha = plotTable[plotTable[, viralNameColumn] != factorOrder[1],]

p = ggplot(plotTable_nonAlpha, aes(x= get(viralNameColumn), y = ID, size = value, fill = area))
p + geom_point(shape = 21) + theme_bw() +
  scale_size_continuous(range=c(2,15)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(viralNameColumn) +
  ggtitle(valueColumn)

