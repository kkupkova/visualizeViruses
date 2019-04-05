rm(list=ls())
library(data.table)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(gplots)
library(RColorBrewer)
############# DATA PREPARATION ######
# uploda metadata
metadata = read.csv("metadata.csv")
colnames(metadata)[1] = "ID"
metadata$ID = as.character(metadata$ID)

# path to the folder with tsv files
filePath = "results_with_taxa_names/"

# get the tax IDs in all of the TSV files - starts with virus, ends with .tsv
tsvFiles = list.files(path = filePath)

# get the ID numbers from the TSV files
ID = regmatches(tsvFiles, gregexpr("[[:digit:]]+", tsvFiles))
ID = unlist(ID)

# filter out only the used metadata
metadata = metadata[metadata$ID %in% ID, ] 


# upload files and make the tables
tsvFiles = list.files(path = filePath, full.names = T)
for (i in 2:length(tsvFiles)){
  
  keepCols = c("Genome", "Final Guess", "Final Best Hit", "Final Best Hit Read Numbers", "Final High Confidence Hits", "Organism")
  if(i == 2){
    
    #upload table - select interesting columns - give them ID - merge by Genome and Organism
    firstFile = fread(tsvFiles[i-1], skip = 1, sep = "\t", header = T) %>% 
      select(keepCols)
    firstID =  unlist(regmatches(tsvFiles[i-1], gregexpr("[[:digit:]]+", tsvFiles[i-1])))
    
    secondFile = fread(tsvFiles[i], skip = 1, sep = "\t", header = T)%>% 
      select(keepCols)
    secondID =  unlist(regmatches(tsvFiles[i], gregexpr("[[:digit:]]+", tsvFiles[i])))
    
    setkey(firstFile, "Genome", "Organism")
    setkey(secondFile, "Genome", "Organism")
    
    origFirstColnames = colnames(firstFile)[!colnames(firstFile) %in% c("Genome", "Organism")]
    colnames(firstFile)[!colnames(firstFile) %in% c("Genome", "Organism")] = paste(origFirstColnames, firstID, sep = ".")
    
    origSecondColnames = colnames(secondFile)[!colnames(secondFile) %in% c("Genome", "Organism")]
    colnames(secondFile)[!colnames(secondFile) %in% c("Genome", "Organism")] = paste(origSecondColnames, secondID, sep = ".")
    
    myTable = merge(firstFile, secondFile, all = TRUE)
    
  } else {
    newFile = fread(tsvFiles[i], skip = 1, sep = "\t", header = T)%>% 
      select(keepCols)
    newID = unlist(regmatches(tsvFiles[i], gregexpr("[[:digit:]]+", tsvFiles[i])))
    
    setkey(newFile, "Genome", "Organism")
    
    origNewColnames = colnames(newFile)[!colnames(newFile) %in% c("Genome", "Organism")]
    colnames(newFile)[!colnames(newFile) %in% c("Genome", "Organism")] = paste(origNewColnames, newID, sep = ".")
    
    myTable = merge(myTable, newFile, all = TRUE)
    
  }
}
  


myTable[is.na(myTable)] = 0
columnOptions = origNewColnames


# separate the table into separate tables and give the columns the kids' IDs
#--final guess table
finalGuess = myTable %>% 
  select(Genome, Organism, starts_with("Final Guess"))

colnames(finalGuess)[-c(1,2)] = c(unlist(regmatches(colnames(finalGuess), gregexpr("[[:digit:]]+", colnames(finalGuess)))))

#--final best hit
finalBestHit = myTable %>% 
  select(Genome, Organism, starts_with("Final Best Hit."))

colnames(finalBestHit)[-c(1,2)] = c(unlist(regmatches(colnames(finalBestHit), gregexpr("[[:digit:]]+", colnames(finalBestHit)))))


#-final best hit read number table
finalBestHitReads = myTable %>% 
  select(Genome, Organism, starts_with("Final Best Hit Read Numbers"))

colnames(finalBestHitReads)[-c(1,2)] = c(unlist(regmatches(colnames(finalBestHitReads), gregexpr("[[:digit:]]+", colnames(finalBestHitReads)))))

rm(list=setdiff(ls(), c("metadata", "finalGuess", "finalBestHit", "finalBestHitReads")))



################## PLOTTING ############################################################

# ******** function plotVirus ********
# input: *tableToPlot - which of the results do you want plotted (e.g. FinalBestHit or FinalGuess)
#        *criteriaTable - which table is used to filter the plotted genomes (e.g. with finalBestHitReads and the
#                      default threshold minReadThreshold = 100, Organisms where at least one had 100 or more 
#                      reads mapped will be kept, the rest is tossed)
#         *metadata - table with metadata
#         *minReadThreshold - used for the filtering of organisms with at least some abundance
#         *lazOrdered - should the IDs in the final plot be organized by laz score? T/F
#         *savePlots - should the resulting plots be saved as pdf? T/F (will create a folder "figures")

plotVirus = function(tableToPlot, 
                     criteriaTable=finalBestHitReads,
                     metadata=metadata,
                     minReadThreshold=100, 
                     maxPerPlot = 20,
                     lazOrdered=F, 
                     savePlots=T)
{
  if (savePlots & !file.exists(file.path(getwd(), "figures"))){
    dir.create(file.path(getwd(), "figures"))
  }
  
  # 1) == filtering and sorting ==
  # select only the viruses which pass given criteria - default is to take the finalBestHitReadsTable and filter
  # out ony the rows where at least one virus has at least 100 mapped reads
  filteredGenomes = criteriaTable[apply(criteriaTable[,-c(1,2)], 1, function(x) !all(x < minReadThreshold)),1]$Genome
  
  # take the selected table to plot - extract the viruses selected in previous step and sort according to the 
  # mean abundance across samples
  plotTable = tableToPlot %>% 
    filter(Genome %in% filteredGenomes) %>% 
    group_by(Genome, Organism) %>% 
    ungroup() %>% 
    mutate(meanVal = rowMeans(.[,3:ncol(finalBestHit)])) %>% 
    gather(key = "ID", value = "abundance", 3:ncol(tableToPlot)) %>% 
    left_join(metadata, by = "ID") %>% 
    mutate_if(is.character, as.factor) 
  
  
  # 2) == plot abundance across samples in facets ==
  
  # get IDs sorted by the mean abundance - to plot the facets from the most important to the least important
  sortedID = plotTable %>% 
    select(Genome, Organism, meanVal) %>% 
    distinct() %>% 
    arrange(desc(meanVal)) %>% 
    mutate(Organism = as.character(Organism)) %>% 
    mutate(Organism = factor(Organism, Organism))
  
  # relevel the plo table according to the sorted IDs
  plotTable = plotTable %>% 
    mutate(Organism = factor(Organism, levels(sortedID$Organism)))
  i = 0
  while(!is.null(sortedID)){
    i = i+1
    if (nrow(sortedID) > maxPerPlot){
      
      # was the function passed and argument lazOrdered?
      if(lazOrdered==T){
        # LAZ ordered
        p = plotTable %>% 
          filter(plotTable$Genome %in% sortedID$Genome[1:maxPerPlot]) %>% 
          ggplot(aes(x= reorder(ID, laz), y = abundance, color = site))
      }else{
        # just groups
        p = plotTable %>% 
          filter(plotTable$Genome %in% sortedID$Genome[1:maxPerPlot]) %>%
          ggplot(aes(x= ID, y = abundance, color = site))
      }
      
      p = p + geom_point() + theme_bw() + facet_wrap(~Organism, scales = "free") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      print(p)
      #if save plot mode is on- save the plot
      if(savePlots){
        ggsave(paste0("figures/AbundanceFacetPlot_", i, ".pdf"), width = 40, height = 20, units = "cm")
      }
      
      sortedID = sortedID[-seq(1,20),]
    }else{
      
      # was the function passed and argument lazOrdered?
      if(lazOrdered==T){
        # LAZ ordered
        p = plotTable %>% 
          filter(plotTable$Genome %in% sortedID$Genome) %>% 
          ggplot(aes(x= reorder(ID, laz), y = abundance, color = site))
      }else{
        # just groups
        p = plotTable %>% 
          filter(plotTable$Genome %in% sortedID$Genome) %>%
          ggplot(aes(x= ID, y = abundance, color = site))
      }
      
      p = p + geom_point() + theme_bw() + facet_wrap(~Organism, scales = "free") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      
      print(p)
      #if save plot mode is on- save the plot
      if(savePlots){
        ggsave(paste0("figures/AbundanceFacetPlot_", i, ".pdf"), width = 40, height = 20, units = "cm")
      }
      
      sortedID = NULL
    }
  }
  
  
  
  # 3) == make circle plots ==
  
  # plot as is
  p= ggplot(plotTable,aes(x= Organism, y = ID, size = abundance, fill = site))
  
  p = p + geom_point(shape = 21) + theme_bw() +
    scale_size_continuous(range=c(0.2,15)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Circle plot - unscaled")
  print(p)
  #if save plot mode is on- save the plot
  if(savePlots){
    ggsave("figures/DotsOnGridPlot_unscaled.pdf", width = 40, height = 20, units = "cm")
  }
  
  # Filter out the most abundant
  p = plotTable %>% 
    filter(Organism != levels(plotTable$Organism)[1]) %>% 
    ggplot(aes(x= Organism, y = ID, size = abundance, fill = site))
  
  p = p + geom_point(shape = 21) + theme_bw() +
    scale_size_continuous(range=c(0.2,15)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Circle plot - unscaled / Human alphherpesvirus 1 fltered out")
  print(p)
  
  #if save plot mode is on- save the plot
  if(savePlots){
    ggsave("figures/DotsOnGridPlot_unscaled_topOrganismFilteredOut.pdf", width = 40, height = 20, units = "cm")
  }
  
  #  do z-score normalization
  p = plotTable %>% 
    group_by(Organism) %>% 
    mutate(z_score = scale(abundance)) %>% 
    ggplot(aes(x= Organism, y = ID, size = z_score, fill = site, color= site))
  
  p = p + geom_point(shape = 21,alpha = 0.8) + theme_bw() +
    scale_size_continuous(range=c(0.1,7)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
    ggtitle("Circle plot - z-scaled")
  print(p)
  #if save plot mode is on- save the plot
  if(savePlots){
    ggsave("figures/DotsOnGridPlot_zScaled.pdf", width = 40, height = 20, units = "cm")
  }
  
  
}

graphics.off()
plotVirus(finalBestHit, finalBestHitReads, metadata, savePlots = T)
plotVirus(finalGuess, finalBestHitReads, metadata, savePlots = T)

#plotVirus(finalBestHit, finalBestHit, metadata = metadata, minReadThreshold = 0.005)
#plotVirus(finalGuess, finalGuess, metadata = metadata, minReadThreshold = 0.005)




# =================HEATMAP===================

filteredGenomes = finalBestHitReads[apply(finalBestHitReads[,-c(1,2)], 1, function(x) !all(x < 100)),1]$Genome

heatmapTable =finalBestHit %>% 
  filter(Genome %in% filteredGenomes)
rownames(heatmapTable) = heatmapTable$Organism

heatmapTable =finalGuess %>% 
  filter(Genome %in% filteredGenomes)
rownames(heatmapTable) = heatmapTable$Organism

# sort IDs bu LAZ score
sortedIDs = metadata$ID[order(metadata$laz)]
heatmapTable = as.matrix(heatmapTable[,3:ncol(heatmapTable)])
heatmapTable = heatmapTable[,sortedIDs]
#heatmap colors: https://davetang.org/muse/2010/12/06/making-a-heatmap-with-r/
graphics.off()  
par(mar=c(10,4,4,10))
colfunc <- colorRampPalette(c("black", "white", "red"))

# plot pink heatmap - default clustering
heatmap.2(heatmapTable, scale = "row", col=brewer.pal(11,"RdBu"),
          cexRow=0.7,cexCol=1,trace="none",srtCol=90, margins=c(7,12), 
          Colv = F, dendrogram = "row")

# use different clustering : 
distF = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
clustM = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" ,  "centroid")
for (distFunc in distF){
  for (clustMet in clustM){
    pdf(paste0("figures/heatmaps/",distFunc,"_", clustMet,".pdf"), height=10, width=5)
    heatmap.2(heatmapTable, scale = "row", col=brewer.pal(11,"RdBu"),
              cexRow=0.7,cexCol=1,trace="none",srtCol=90, margins=c(7,12), 
              dendrogram = "row", #Colv = F,
              distfun=function(x) dist(x, method=distFunc), 
              hclustfun=function(x) hclust(x, method= clustMet))
    dev.off()
  }
  
}






