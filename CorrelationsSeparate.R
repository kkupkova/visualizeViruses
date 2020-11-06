library(cowplot)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
# 1) run DataPreparation script
finalBestHit = as.data.frame(finalBestHit)
finalBestHitReads = as.data.frame(finalBestHitReads)
finalGuess = as.data.frame(finalGuess)




correlate = function(tableToCorrelate, vectorToCorrelate, method = "pearson"){
  # select just the counts to correlate and make sure that 
  # columns are in the correct order
  editedTable = as.data.frame(tableToCorrelate[, -c(1,2)])
  
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
# -------------------- Bangladesh only -----------------------
Bangladesh = cbind(tableToCorrelate[,c(1,2)], tableToCorrelate[,metadata$ID[metadata$site=="Icddr,b"]]) 
corVals = correlate(Bangladesh, metadata$laz[metadata$site=="Icddr,b"])

# get the viruses, which correlate with the LAZ most in pozitive or negative direction
sorted = sort(corVals, decreasing = T)

# add metadata to the table and make it ggplot usable
rownames(finalBestHit)  = finalBestHit$Organism
corPlot = finalBestHit[names(sorted),] %>% 
  add_column(corVal = sorted) %>% 
  gather("ID", "abundance", 3:(ncol(finalBestHit))) %>% 
  left_join(metadata, by = "ID") %>% 
  filter(site == "Icddr,b")
corPlot$percent = corPlot$abundance*100

# relevel the factors based on the correlation value
corPlot$Organism = factor(corPlot$Organism, levels = finalBestHit[names(sorted),"Organism"])
corPlot$ID = factor(corPlot$ID, levels = metadata$ID[order(metadata$laz)])

# filter out the 1007 outlier
#corPlot = corPlot %>% 
#  filter(ID != "1007")

# make facet plot
p = ggplot(corPlot,aes(laz, percent))
p + geom_point(color = "darkred", size = 2) + facet_wrap(~Organism, scales = "free", nrow = 8)+
  ylab("Mapped proportions [%]") + xlab("LAZ") + ggtitle("Bangladesh") +  
  theme(strip.text.x = element_text(size = 13)) + geom_smooth(method = "lm", se = F, color = "gray75")
ggsave(filename = "figures/corPlot_individualAbundancesVsLAZ_BangladesOnly_FinalBEstHit0_0005.pdf", width = 20, height = 20)

#plot the correlation values
a = ggplot(corPlot, aes(Organism, corVal))
a + geom_point(color = "darkred") + theme(axis.text.x = element_text(hjust = 1, angle = 90)) + geom_hline(yintercept = 0, color = "gray60") +
  ggtitle("Bangladesh") + ylab("Pearson correlation") + xlab(" ")
ggsave(filename = "figures/corVals_BangladeshOnlyFinalBestHit0_0005.pdf")

#------------** bar plot or pie chart **--------

n <- length(levels(corPlot$Organism))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p = ggplot(corPlot, aes(ID, y = percent, fill= Organism))
p + geom_bar(stat = "identity") + guides(fill=guide_legend(ncol=1)) + 
  scale_fill_manual(values=sample(col_vector, n))

# -------------------UVA only------------------
UVA = cbind(tableToCorrelate[,c(1,2)], tableToCorrelate[,metadata$ID[metadata$site=="Uva"]]) 
corVals = correlate(UVA, metadata$laz[metadata$site=="Uva"])

# get the viruses, which correlate with the LAZ most in pozitive or negative direction
sorted = sort(corVals, decreasing = T)

# add metadata to the table and make it ggplot usable
rownames(finalBestHit)  = finalBestHit$Organism
corPlot = finalBestHit[names(sorted),] %>% 
  add_column(corVal = sorted) %>% 
  gather("ID", "abundance", 3:(ncol(finalBestHit))) %>% 
  left_join(metadata, by = "ID") %>% 
  filter(site == "Uva")
corPlot$percent = corPlot$abundance*100

# relevel the factors based on the correlation value
corPlot$Organism = factor(corPlot$Organism, levels = finalBestHit[names(sorted),"Organism"])
corPlot$ID = factor(corPlot$ID, levels = metadata$ID[order(metadata$laz)])
# make facet plot
p = ggplot(corPlot,aes(laz, percent))
p + geom_point(size = 2, color = "darkblue") + facet_wrap(~Organism, scales = "free")+
  ylab("Mapped proportions [%]") + xlab("LAZ") + ggtitle("UVA") +  
  theme(strip.text.x = element_text(size = 13)) + geom_smooth(method = "lm", color = "gray75", se = F)
ggsave(filename = "figures/corPlot_individualAbundancesVsLAZ_UVAOnly_FinalBEstHit0_0005.pdf", width = 20, height = 20)

#plot the correlation values
a = ggplot(corPlot, aes(Organism, corVal))
a + geom_point(color = "darkblue") + theme(axis.text.x = element_text(hjust = 1, angle = 90)) + 
  geom_hline(yintercept = 0, color = "gray60") + ggtitle("UVA")+ ylab("Pearson correlation") + xlab(" ")
ggsave(filename = "figures/corVals_UVAOnlyFinalBestHit0_0005.pdf")



#### BOX PLOT UVA vs Bangladesh
# plot box plots
p = tableToCorrelate %>% 
  gather("ID", "abundance", 3:ncol(finalBestHit)) %>% 
  left_join(metadata, by = "ID") %>% 
  ggplot(aes(x = site, y = abundance))
p + geom_boxplot(aes(color= site), alpha = 0.6) + geom_point() + facet_wrap(~Organism, scales = "free")+ 
  stat_compare_means()+ 
  scale_color_manual(values=c("darkred", "darkblue"))
ggsave(filename = "figures/boxplotWithPVals_BangaldeshVsUVA_finalBestHit0_0005.pdf", width = 20, height = 15)



# ------------------------- correlation heat map ------------

library(corrplot)
library(gplots)
library(RColorBrewer)
library(ggbiplot)
corMat = cor(tableToCorrelate[, 3:ncol(tableToCorrelate)])

Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)
heatmap.2(corMat, scale = "none", col=Colors, density.info="none", trace="none")
heat = heatmap.2(corMat, scale = "none", col=Colors, density.info="none", trace="none")














