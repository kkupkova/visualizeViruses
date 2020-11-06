library(cowplot)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ggpmisc)
# 1) run DataPreparation script
finalBestHit = as.data.frame(finalBestHit)
finalBestHitReads = as.data.frame(finalBestHitReads)

# eliminate sample 1007 - it behaves weird
finalBestHit = finalBestHit[,colnames(finalBestHit) != "1007"]
finalBestHitReads = finalBestHitReads[,colnames(finalBestHitReads) != "1007"]
metadata = metadata[metadata$ID != "1007",]

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

# upload the taxa table with the info - this is based on the organisms that passed tableToCorrelate 
taxaTable = read.delim("/Users/lilcrusher/epigen_mal/viruses/hierarchy/taxID_organismName_filtered_host.csv", 
                       sep = ",", as.is = T, na.strings=c("","NA"))

# replace NA with x and change column names
taxaTable[is.na(taxaTable)]  = "x"
colnames(taxaTable)  = c("host", "family1", "taxID", "Organism", "species", "genus","family2", "keep")

# convert the first letter to uppercase- just precaution
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
taxaTable = as.data.frame(apply(taxaTable, 2, function(x) firstup(x))) 

# filter out the undesired organisms and recalculate the proportions
rownames(finalBestHitReads) = finalBestHitReads$Organism
keepOrganisms = taxaTable %>% 
  filter(keep != "Omit")
finalBestHitReads = finalBestHitReads[as.character(keepOrganisms$Organism),]

# get the number of the mapped reads after filtering
mappedReads = colSums(finalBestHitReads[-c(1,2)])

# recalculate the proportions
helpTab = sweep(finalBestHitReads[,-c(1,2)], 2, mappedReads, "/")
finalBestHit = cbind(finalBestHitReads[,c(1,2)], helpTab)

# run the function
# -------------------- Bangladesh only -----------------------
tableToCorrelate = finalBestHit
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


# make facet plot
p = ggplot(corPlot,aes(laz, percent))
p + geom_point(color = "darkred", size = 2) + facet_wrap(~Organism, scales = "free")+
  ylab("Mapped proportions [%]") + xlab("LAZ") + ggtitle("Bangladesh") +  
  theme(strip.text.x = element_text(size = 13)) + geom_smooth(method = "lm", se = F, color = "gray75")
#ggsave(filename = "figures/fiteredOrganisms_corPlot_individualAbundancesVsLAZ_BangladesOnly_FinalBEstHit0_0005.pdf", width = 20, height = 20)

#plot the correlation values
a = ggplot(corPlot, aes(Organism, corVal))
a + geom_point(color = "darkred") + theme(axis.text.x = element_text(hjust = 1, angle = 90)) + geom_hline(yintercept = 0, color = "gray60") +
  ggtitle("Bangladesh") + ylab("Pearson correlation") + xlab(" ")
##ggsave(filename = "figures/fiteredOrganisms_corVals_BangladeshOnlyFinalBestHit0_0005.pdf")



# ********* collapse the groups *************

# factor to char, make organism a rowname to make sure we have the same order and append the genus column to the finalBestHit table
# if there is no genus information, replace it with family2
keepOrganisms = taxaTable %>% 
  filter(keep != "Omit") %>%
  mutate_if(is.factor, as.character)
keepOrganisms$genus[keepOrganisms$genus == "X"] = keepOrganisms$family2[keepOrganisms$genus == "X"]
rownames(keepOrganisms) = keepOrganisms$Organism

keepOrganisms = keepOrganisms[finalBestHit$Organism,]

# get all the categories
genusLevels = levels(factor(keepOrganisms$genus))

# go through the genus level - if there are more organisms - sum the percentage - give the row genus name
for (i in genusLevels){
  organisms = keepOrganisms$Organism[keepOrganisms$genus == i]
  helpTab = Bangladesh[organisms,]
  
  if (nrow(helpTab) > 1){
   sumVec = colSums(helpTab[,-c(1,2)]) 
   helpTab = as.data.frame(t(sumVec))
  }else{
    helpTab = helpTab[,-c(1,2)]
  }
  rownames(helpTab) = i
  
  if (i == genusLevels[1]) {
    collapsedTable = helpTab
  }else{
    collapsedTable = rbind(collapsedTable, helpTab)
  }
  rm(helpTab)
}


#collapsedTable = collapsedTable[,colnames(collapsedTable) != "1007"]
# make sure we are correlating the corresponding columns
rownames(metadata) = metadata$ID
usedMetadata = metadata[colnames(collapsedTable),]



corVals = apply(collapsedTable, 1, function(x) cor(x, usedMetadata$laz, method = "p"))

# get the viruses, which correlate with the LAZ most in pozitive or negative direction
sorted = sort(corVals, decreasing = T)

# add metadata to the table and make it ggplot usable
corPlot = collapsedTable[names(sorted),] %>% 
  rownames_to_column(var= "genus") %>% 
  add_column(corVal = sorted) %>% 
  gather("ID", "abundance", 2:(ncol(collapsedTable)+1)) %>% 
  left_join(metadata, by = "ID") %>% 
  mutate(percent = abundance * 100)

# relevel the factors based on the correlation value
corPlot$genus = factor(corPlot$genus, levels = names(sorted))
corPlot$ID = factor(corPlot$ID, levels = usedMetadata$ID[order(usedMetadata$laz)])


# make facet plot
formula = y ~ x 
p = ggplot(corPlot,aes(laz, percent))
p + geom_point(color = "darkred", size = 2) + facet_wrap(~genus, scales = "free")+
  ylab("Mapped proportions [%]") + xlab("LAZ") + ggtitle("Bangladesh") +  
  theme(strip.text.x = element_text(size = 13)) #+ 
  #geom_smooth(method = "lm", se = F, color = "gray75", formula = formula) +
  #stat_poly_eq(aes(label = paste(..rr.label..)), 
               # label.x.npc = "right", label.y.npc = 0.95,
               # formula = formula, parse = TRUE, size = 3) 
#ggsave(filename = "figures/genus_no1007_fiteredOrganisms_corPlot_individualAbundancesVsLAZ_BangladesOnly_FinalBEstHit0_0005.pdf", width = 16, height =9)

#plot the correlation values
a = ggplot(corPlot, aes(genus, corVal))
a + geom_point(color = "darkred") + theme(axis.text.x = element_text(hjust = 1, angle = 90)) + 
  geom_hline(yintercept = 0, color = "gray60") +
  ggtitle("Bangladesh") + ylab("Pearson correlation") + xlab(" ")
#ggsave(filename = "figures/genus_no1007_fiteredOrganisms_corVals_BangladeshOnlyFinalBestHit0_0005.pdf")


# plot onlu the interesting ones: 

names(corVals)

interesting = c("Cardiovirus", "Hepacivirus", "Retroviridae", "Orthohepevirus", "Rotavirus")

corPlotSelected = corPlot[corPlot$genus %in% interesting,]
corPlotSelected$genus = factor(corPlotSelected$genus, levels = interesting)

# make facet plot
formula = y ~ x 
p = ggplot(corPlotSelected,aes(laz, percent))
p + geom_point(color = "darkred", size = 2) + facet_wrap(~genus, scales = "free", ncol = 5)+
  ylab("Mapped proportions [%]") + xlab("LAZ") + ggtitle("Icddr,b") +  
  theme(strip.text.x = element_text(size = 13)) + 
  geom_smooth(method = "lm", se = F, color = "gray75", formula = formula) +
  stat_poly_eq(aes(label = paste(..rr.label..)), 
               label.x.npc = "right", label.y.npc = 0.95,
               formula = formula, parse = TRUE, size = 3) 
#ggsave(filename = "figures/topSelected_genus_no1007_fiteredOrganisms_corPlot_individualAbundancesVsLAZ_BangladesOnly_FinalBEstHit0_0005.pdf", width = 15, height =3)


#------------** bar plot or pie chart **--------

n <- length(levels(corPlot$genus))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


means = apply(collapsedTable, 1, mean)
means = sort(means)
corPlot$genus = factor(corPlot$genus, levels = names(means))


myColors = sample(col_vector, n)

myColors = c("#E5C494", "#F2F2F2", "#D9D9D9", "#FCCDE5", "#F4CAE4", "#CBD5E8", "#7FC97F",
             "#FFFFCC", "#F781BF", "#8DA0CB", "#E6F5C9", "#E31A1C", "#386CB0", "#D95F02",
             "#F0027F", "#FC8D62", "#B3CDE3", "#FFF2AE")
# plot the bars ordered by LAZ - best to worst
p = ggplot(corPlot, aes(reorder(ID, laz), y = percent, fill= genus))
p + geom_bar(stat = "identity", alpha = 0.7) + guides(fill=guide_legend(ncol=1)) + 
  scale_fill_manual(values=myColors) +
  ylab("Percent relative abundance") + 
  xlab("ID")
#ggsave("figures/barplot_genusLevel_sortedByLAZ_worstToBetter.pdf", height = 5, width = 8)

# plot the bars on LAZ
p = ggplot(corPlot, aes(laz, y = percent, fill= genus))
p + geom_bar(stat = "identity", alpha = 0.7) + 
  guides(fill=guide_legend(ncol=1)) + 
  scale_fill_manual(values=myColors)+
  ylab("Percent relative abundance") + 
  xlab("LAZ")
#ggsave("figures/barplot_genusLevel_vsLAZ.pdf", height = 5, width = 8)


p = corPlot %>% 
  filter(genus != "Simplexvirus") %>% 
  ggplot(aes(reorder(ID, laz), y = percent, fill= genus))

p + geom_bar(stat = "identity", alpha = 0.5) + 
  guides(fill=guide_legend(ncol=1)) + 
  scale_fill_manual(values=myColors)+
  ylab("Percent relative abundance") + 
  xlab("ID")
#ggsave("figures/barplot_genusLevel_simplexRemoved_sortedByLAZ_worstToBetter.pdf", height = 5, width = 8)

# make the bars vs LAZ just to see
p = corPlot %>% 
  filter(genus != "Simplexvirus") %>% 
  ggplot(aes(laz, y = percent, fill= genus))

p + geom_bar(stat = "identity", alpha = 0.5) + guides(fill=guide_legend(ncol=1)) + 
  scale_fill_manual(values=myColors)+
  ylab("Percent relative abundance") + 
  xlab("LAZ")
#ggsave("figures/barplot_genusLevel_simplexRemoved_plotOnLAZ.pdf", height = 5, width = 8)


# make dot plot

orderNames = names(sort(rowMeans(collapsedTable), decreasing = T))
orderID = metadata[colnames(collapsedTable),] %>% 
  arrange((desc(laz))) %>% 
  unite("IDlabel", ID, sex, laz, sep = "_", remove = F)

orderedCollapsed = collapsedTable[names(sort(rowMeans(collapsedTable), decreasing = T)),]

orderedCollapsed = orderedCollapsed %>% 
  rownames_to_column(var = "genus") %>% 
  gather(key = "ID", value  = "abundance", 2:(ncol(orderedCollapsed)+1)) %>% 
  left_join(orderID, by = "ID") %>% 
  mutate(genus = factor(genus, levels = orderNames)) %>% 
  mutate(ID = factor(ID, levels = orderID$ID)) %>% 
  mutate(IDlabel = factor(IDlabel, levels = orderID$IDlabel))

p = ggplot(filter(orderedCollapsed, genus != "Simplexvirus"), aes(x = genus, y = IDlabel, size = abundance, fill = site)) + 
  geom_point(shape = 21) + theme_bw() +
  scale_size_continuous(range=c(2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
q = ggplot(filter(orderedCollapsed, genus == "Simplexvirus"), aes(x = genus, y = IDlabel, size = abundance, fill = site))+ 
  geom_point(shape = 21) + theme_bw() +
  scale_size_continuous(range=c(2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# add color scaling to the dots
p = ggplot(filter(orderedCollapsed, genus != "Simplexvirus"), 
           aes(x = genus, y = IDlabel, size = abundance, fill = log(abundance))) + 
  geom_point(shape = 21) + theme_bw() +
  scale_size_continuous(range=c(2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_gradient2(midpoint = -8)
p
q = ggplot(filter(orderedCollapsed, genus == "Simplexvirus"), aes(x = genus, y = IDlabel, size = abundance, fill = site))+ 
  geom_point(shape = 21) + theme_bw() +
  scale_size_continuous(range=c(2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))











