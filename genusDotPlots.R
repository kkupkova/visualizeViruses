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
  helpTab = finalBestHit[organisms,]
  
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





# make dot plot

orderNames = names(sort(rowMeans(collapsedTable), decreasing = T))
orderID = metadata[colnames(collapsedTable),] %>% 
  arrange(site,(desc(laz))) %>% 
  unite("IDlabel", ID, sex, laz, sep = "_", remove = F)

orderedCollapsed = collapsedTable[names(sort(rowMeans(collapsedTable), decreasing = T)),]

orderedCollapsed = orderedCollapsed %>% 
  rownames_to_column(var = "genus") %>% 
  gather(key = "ID", value  = "abundance", 2:(ncol(orderedCollapsed)+1)) %>% 
  left_join(orderID, by = "ID") %>% 
  mutate(genus = factor(genus, levels = orderNames)) %>% 
  mutate(ID = factor(ID, levels = orderID$ID)) %>% 
  mutate(IDlabel = factor(IDlabel, levels = orderID$IDlabel))



# Filter out the most abundant and plot separately
p = ggplot(filter(orderedCollapsed, genus != "Simplexvirus"), aes(x = genus, y = ID, size = abundance, fill = site)) + 
  geom_point(shape = 21, alpha = 0.8) + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15)) +
  xlab(" ") + ylab(" ") + scale_fill_discrete(guide = 'none')
p

q = ggplot(filter(orderedCollapsed, genus == "Simplexvirus"), aes(x = genus, y = ID, size = abundance, fill = site))+ 
  geom_point(shape = 21, alpha = 0.7) + theme_bw() +
  scale_size_continuous(range=c(5,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15),
        axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab(" ") + ylab(" ")




# plot alpha herpes next to the first plot
ggdraw() +
  draw_plot(p + theme(legend.justification = "bottom"), 0, 0, 0.85, 1) +
  draw_plot(q + theme(legend.justification = "bottom") , 0.85, 0.16, 0.13, 0.84) +
  draw_plot_label(c(" ", " "), c(0, 0.5), c(1, 0.92), size = 15)
ggsave("figures/genusDotPlot/UVA_vsBangladesh_LAZsorted.pdf", width = 16, height = 9)


#########
# Filter out the most abundant and plot separately
p = ggplot(filter(orderedCollapsed, genus != "Simplexvirus"), aes(x = genus, y = IDlabel, size = abundance, fill = site)) + 
  geom_point(shape = 21, alpha = 0.8) + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15)) +
  xlab(" ") + ylab(" ") + scale_fill_discrete(guide = 'none')
p

q = ggplot(filter(orderedCollapsed, genus == "Simplexvirus"), aes(x = genus, y = IDlabel, size = abundance, fill = site))+ 
  geom_point(shape = 21, alpha = 0.7) + theme_bw() +
  scale_size_continuous(range=c(5,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15),
        axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab(" ") + ylab(" ")




# plot alpha herpes next to the first plot
ggdraw() +
  draw_plot(p + theme(legend.justification = "bottom"), 0, 0, 0.85, 1) +
  draw_plot(q + theme(legend.justification = "bottom") , 0.85, 0.16, 0.13, 0.84) +
  draw_plot_label(c(" ", " "), c(0, 0.5), c(1, 0.92), size = 15)
ggsave("figures/genusDotPlot/UVA_vsBangladesh_LAZsorted_IDlabel.pdf", width = 16, height = 9)





######### BANGLADESH
# Filter out the most abundant and plot separately
p = ggplot(filter(orderedCollapsed, genus != "Simplexvirus", site == "Icddr,b"), aes(x = genus, y = IDlabel, size = abundance)) + 
  geom_point(shape = 21, alpha = 0.8, fill= "darkred") + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15)) +
  xlab(" ") + ylab(" ") + scale_fill_discrete(guide = 'none')
p

q = ggplot(filter(orderedCollapsed, genus == "Simplexvirus", site == "Icddr,b"), aes(x = genus, y = IDlabel, size = abundance))+ 
  geom_point(shape = 21, alpha = 0.8, fill= "darkred") + theme_bw() +
  scale_size_continuous(range=c(5,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15),
        axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab(" ") + ylab(" ")




# plot alpha herpes next to the first plot
ggdraw() +
  draw_plot(p + theme(legend.justification = "bottom"), 0, 0, 0.85, 1) +
  draw_plot(q + theme(legend.justification = "bottom") , 0.85, 0.16, 0.13, 0.84) +
  draw_plot_label(c(" ", " "), c(0, 0.5), c(1, 0.92), size = 15)
ggsave("figures/genusDotPlot/BangladeshOnly_LAZsorted_IDlabel.pdf", width = 16, height = 9)



# Filter out the most abundant and plot separately
p = ggplot(filter(orderedCollapsed, genus != "Simplexvirus", site == "Icddr,b"), aes(x = genus, y = ID, size = abundance)) + 
  geom_point(shape = 21, alpha = 0.8, fill= "darkred") + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15)) +
  xlab(" ") + ylab(" ") + scale_fill_discrete(guide = 'none')
p

q = ggplot(filter(orderedCollapsed, genus == "Simplexvirus", site == "Icddr,b"), aes(x = genus, y = ID, size = abundance))+ 
  geom_point(shape = 21, alpha = 0.8, fill= "darkred") + theme_bw() +
  scale_size_continuous(range=c(5,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15),
        axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab(" ") + ylab(" ")




# plot alpha herpes next to the first plot
ggdraw() +
  draw_plot(p + theme(legend.justification = "bottom"), 0, 0, 0.85, 1) +
  draw_plot(q + theme(legend.justification = "bottom") , 0.85, 0.16, 0.13, 0.84) +
  draw_plot_label(c(" ", " "), c(0, 0.5), c(1, 0.92), size = 15)
ggsave("figures/genusDotPlot/BangladeshOnly_LAZsorted.pdf", width = 16, height = 9)



############# UVA

# Filter out the most abundant and plot separately
p = ggplot(filter(orderedCollapsed, genus != "Simplexvirus", site == "Uva"), aes(x = genus, y = IDlabel, size = abundance)) + 
  geom_point(shape = 21, alpha = 0.8, fill= "darkblue") + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15)) +
  xlab(" ") + ylab(" ") + scale_fill_discrete(guide = 'none')
p

q = ggplot(filter(orderedCollapsed, genus == "Simplexvirus", site == "Uva"), aes(x = genus, y = IDlabel, size = abundance))+ 
  geom_point(shape = 21, alpha = 0.8, fill= "darkblue") + theme_bw() +
  scale_size_continuous(range=c(5,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15),
        axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab(" ") + ylab(" ")




# plot alpha herpes next to the first plot
ggdraw() +
  draw_plot(p + theme(legend.justification = "bottom"), 0, 0, 0.85, 1) +
  draw_plot(q + theme(legend.justification = "bottom") , 0.85, 0.16, 0.13, 0.84) +
  draw_plot_label(c(" ", " "), c(0, 0.5), c(1, 0.92), size = 15)
ggsave("figures/genusDotPlot/UVAOnly_LAZsorted_IDlabel.pdf", width = 16, height = 9)



# Filter out the most abundant and plot separately
p = ggplot(filter(orderedCollapsed, genus != "Simplexvirus", site == "Uva"), aes(x = genus, y = ID, size = abundance)) + 
  geom_point(shape = 21, alpha = 0.8, fill= "darkblue") + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15)) +
  xlab(" ") + ylab(" ") + scale_fill_discrete(guide = 'none')
p

q = ggplot(filter(orderedCollapsed, genus == "Simplexvirus", site == "Uva"), aes(x = genus, y = ID, size = abundance))+ 
  geom_point(shape = 21, alpha = 0.8, fill= "darkblue") + theme_bw() +
  scale_size_continuous(range=c(5,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15),
        axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab(" ") + ylab(" ")




# plot alpha herpes next to the first plot
ggdraw() +
  draw_plot(p + theme(legend.justification = "bottom"), 0, 0, 0.85, 1) +
  draw_plot(q + theme(legend.justification = "bottom") , 0.85, 0.16, 0.13, 0.84) +
  draw_plot_label(c(" ", " "), c(0, 0.5), c(1, 0.92), size = 15)
ggsave("figures/genusDotPlot/UVAOnly_LAZsorted.pdf", width = 16, height = 9)


