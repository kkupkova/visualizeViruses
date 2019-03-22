tableToPlot = finalBestHit
criteriaTable = finalBestHitReads



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

plotVirus(finalBestHit, finalBestHitReads, metadata = metadata)

