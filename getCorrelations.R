library(cowplot)
library(tidyverse)
# 1) run DataPreparation script
finalBestHit = as.data.frame(finalBestHit)

# tableToCorrelate = as.data.frame(finalBestHit[, -c(1,2)])
# 
# tableToCorrelate = finalBestHit
# vectorToCorrelate = metadata$laz
# method = "pearson"


correlate = function(tableToCorrelate, vectorToCorrelate, method){
  # select just the counts to correlate and make sure that 
  # columns are in the correct order
  editedTable = as.data.frame(tableToCorrelate[, -c(1,2)])
  editedTable = editedTable[, metadata$ID]
  
  correlation =apply(editedTable, 1, function(x) cor(x, vectorToCorrelate, method = method))
  
  
  plotTable = cbind(tableToCorrelate, correlation)
  
  
   p = ggplot(plotTable, aes(reorder(Organism, -correlation), correlation))
   p= p + geom_point() + theme(axis.title.x =element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            text = element_text(size=10)) + 
     geom_hline(yintercept = 0)+ xlab(" ") + ylab("correlation of abundance with LAz")
   print(p)
   
   # plot the correlations with the organism name
   q = plotTable %>% 
     gather("ID", "abundance", 3:(ncol(plotTable)-1)) %>%
     left_join(metadata, by = "ID") %>% 
     ggplot(aes(x= reorder(Organism, -correlation), y = ID, size = abundance, fill = site))
   
   # plot the dot plot
   q = q + geom_point(shape = 21, alpha = 0.8) + theme_bw() +
     scale_size_continuous(range=c(0.2,15)) + 
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     xlab(" ") + ylab("abundance")
   q
   print(q)
   
   
   # put those together - correlation on top, dot plot at the bottom
   ggdraw() +
     draw_plot(p , 0, 0.7, 0.91, 0.3) +
     draw_plot(q + theme(legend.justification = "bottom") , 0.01, 0, 1, 0.7) +
     draw_plot_label(c(" ", " "), c(0, 0.5), c(1, 0.92), size = 15)
   
  names(correlation) = tableToCorrelate$Organism
  return(correlation)
  
}

# filter the organisms that mapped by at least 0.001 in at least one child
tableToCorrelate = finalBestHit[apply(finalBestHit[,-c(1,2)], 1, function(x) !all(x < 0.0005)),]
# run the function
# !! HERE SELECT THE TBALE TO PLOT 
#corVals = correlate(finalBestHit, metadata$laz, "pearson")  # whole finalBestHitTable
corVals = correlate(tableToCorrelate, metadata$laz, "pearson") # only when at least passes 0.001 threshold

# get the viruses, which correlate with the LAZ most in pozitive or negative direction
sorted = sort(corVals, decreasing = T)

# add metadata to the table and make it ggplot usable
rownames(finalBestHit)  = finalBestHit$Organism
corPlot = finalBestHit[names(sorted),] %>% 
  add_column(corVal = sorted) %>% 
  gather("ID", "abundance", 3:(ncol(finalBestHit))) %>% 
  left_join(metadata, by = "ID") 

# relevel the factors based on the correlation value
corPlot$Organism = factor(corPlot$Organism, levels = finalBestHit[names(sorted),"Organism"])
corPlot$ID = factor(corPlot$ID, levels = metadata$ID[order(metadata$laz)])

corPlot$percent = corPlot$abundance*100
# make facet plot
p = ggplot(corPlot,aes(laz, percent))
p + geom_point(aes(color = site)) + facet_wrap(~Organism, scales = "free")+
  ylab("Mapped proportions [%]") + xlab("LAZ")

#plot the correlation values
a = ggplot(corPlot, aes(Organism, corVal))
a + geom_point() + theme(axis.text.x = element_text(hjust = 1, angle = 90)) + geom_hline(yintercept = 0, color = "gray60")





####
# plot box plots
p = finalBestHit %>% 
  gather("ID", "abundance", 3:ncol(finalBestHit)) %>% 
  left_join(metadata, by = "ID") %>% 
  ggplot(aes(x = site, y = abundance))
p + geom_boxplot(aes(fill= site), alpha = 0.6) + geom_point() + facet_wrap(~Organism, scales = "free")



