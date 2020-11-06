# run dataPereparation.R
library(ggbiplot)
library(tidyverse)
library(DESeq2)
library(plotly)

pcaData = finalBestHit
label = pcaData$Organism
pcaData = pcaData[,-c(1,2)]
pcaData = t(pcaData)
colnames(pcaData) = label

PCA = prcomp(pcaData)

ggbiplot(PCA, var.axes = F, labels = rownames(pcaData))


# 1007 is weird - we saw that even previously

pcaDataFilt = pcaData[rownames(pcaData) != "1007",]

PCA = prcomp(pcaDataFilt)

plotPCA = as.data.frame(PCA$x)
plotPCA = plotPCA %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(metadata, by = "ID")

plotPCA$group = factor(plotPCA$group, levels = c("stunted", "at_risk", "control"))

p = ggplot(plotPCA, aes(x = PC1, y = PC2, color = group, size = laz, shape = site))
p + geom_point(alpha = 0.7) + 
  theme_bw() + 
  scale_color_brewer(palette = "Set1")+
  coord_fixed()
#ggsave("figures/PCA_without1007.pdf", width = 6, height = 4)

# plot only bangladesh samples: 

metadataBanglades = metadata %>% 
  filter(site == "Icddr,b") %>% 
  filter(ID != "1007")
pcaBangladesh = pcaData[metadataBanglades$ID,]

PCA = prcomp(pcaBangladesh)

plotPCA = as.data.frame(PCA$x)
plotPCA = plotPCA %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(metadataBanglades, by = "ID")

plotPCA$group = factor(plotPCA$group, levels = c("stunted", "at_risk"))

PCAsummary = summary(PCA)
variance = PCAsummary$importance[2,]

p = ggplot(plotPCA, aes(x = PC1, y = PC2, fill = laz))
p + geom_point(alpha = 1, shape = 21, color = "gray50", size = 7) + 
  theme_bw() + 
  coord_fixed() + 
  xlab(paste0("PC1 (", round(variance["PC1"]*100, digits = 1), " %)"))+
  ylab(paste0("PC2 (", round(variance["PC2"]*100, digits = 1), " %)")) +
  scale_fill_gradient2(midpoint = -3, high = "#034B61", low = "#57233A", mid = "white")

#ggsave("figures/PCA_without1007_BangladeshOnly.pdf", width = 6, height = 4)

p <- plot_ly(plotPCA, x = ~PC1, y = ~PC2, z = ~PC3, color = ~laz, 
             marker = list(size = 10)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))
p
