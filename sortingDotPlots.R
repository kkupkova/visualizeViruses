rm(list=ls())
library(data.table)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(magrittr)
library(cowplot)
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


tableToPlot= finalBestHit
criteriaTable=finalBestHit
metadata = metadata
minReadThreshold = 0.001
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


# get IDs sorted by the mean abundance - to plot the facets from the most important to the least important
sortedID = plotTable %>% 
  select(Genome, Organism, meanVal) %>% 
  distinct() %>% 
  arrange(desc(meanVal)) %>% 
  mutate(Organism = as.character(Organism)) %>% 
  mutate(Organism = factor(Organism, Organism))

# relevel the plot table according to the sorted IDs
plotTable$Organism = factor(plotTable$Organism, levels = sortedID$Organism)

# plot as is
p= ggplot(plotTable,aes(x= Organism, y = ID, size = abundance, fill = site))

p = p + geom_point(shape = 21, alpha = 0.8) + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Circle plot - unscaled") + xlab(" ") + ylab(" ")
p

# Filter out the most abundant
p = plotTable %>% 
  filter(Organism != "Human alphaherpesvirus 1") %>% 
  ggplot(aes(x= Organism, y = ID, size = abundance, fill = site))

p = p + geom_point(shape = 21, alpha = 0.8) + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15)) +
  xlab(" ") + ylab(" ") + scale_fill_discrete(guide = 'none')
p


# Filter out the most abundant
q = plotTable %>% 
  filter(Organism == "Human alphaherpesvirus 1") %>% 
  ggplot(aes(x= Organism, y = ID, size = abundance, fill = site))

q = q + geom_point(shape = 21, alpha = 0.7) + theme_bw() +
  scale_size_continuous(range=c(5,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=15),
        axis.title.y =element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
  xlab(" ") + ylab(" ")
q


# plot alpha herpes next to the first plot
ggdraw() +
  draw_plot(p + theme(legend.justification = "bottom"), 0, 0, 0.85, 1) +
  draw_plot(q + theme(legend.justification = "bottom") , 0.85, 0.12, 0.13, 0.88) +
  draw_plot_label(c(" ", " "), c(0, 0.5), c(1, 0.92), size = 15)

# --------order by log2fold change mean---------

# filter the genomes that passed criteria
# melt the abundance columns
# add metadata to have the site information
# group by organism and site
# and get a mean value 
meanAbundance = tableToPlot %>% 
  filter(Genome %in% filteredGenomes) %>% 
  gather(key = "ID", value = "abundance", 3:ncol(tableToPlot)) %>% 
  left_join(metadata, by = "ID") %>% 
  group_by(Organism, site) %>% 
  summarise(avg = mean(abundance)) 
  #or geometric mean
  #summarise(avg = exp(mean(log(abundance)))) # these are the mean abundances in the groups

# then get fold changes:
lfcSorted = meanAbundance %>% 
  spread(site, avg) %>% 
  set_colnames(c("Organism", "Bangladesh", "UVA")) %>% 
  mutate(log2foldChange = log2(Bangladesh / UVA)) %>% 
  arrange(desc(log2foldChange))


# relevel the plot table according to log2foldChanges
plotTable$Organism = factor(plotTable$Organism, levels = lfcSorted$Organism)

# plot as is
p= ggplot(plotTable,aes(x= Organism, y = ID, size = abundance, fill = site))

p = p + geom_point(shape = 21, alpha = 0.8) + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Circle plot - unscaled / sorted by log2(fold change of mean)")
p

# Filter out the most abundant
p = plotTable %>% 
  filter(Organism != "Human alphaherpesvirus 1") %>% 
  ggplot(aes(x= Organism, y = ID, size = abundance, fill = site))

p = p + geom_point(shape = 21, alpha = 0.8) + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Circle plot - unscaled / Human alphaherpesvirus 1 filtered out / sorted by log2(fold change of mean)")
p

alpha  = plotTable %>% 
  filter(Organism == "Human alphaherpesvirus 1") 


# =================HEATMAP===================

heatmapTable =finalBestHit %>% 
  filter(Genome %in% filteredGenomes)
rownames(heatmapTable) = heatmapTable$Organism

# sort IDs by LAZ score
##sortedIDs = metadata$ID[order(metadata$laz)]
heatmapTable = as.matrix(heatmapTable[,3:ncol(heatmapTable)])
##heatmapTable = heatmapTable[,sortedIDs]
#heatmap colors: https://davetang.org/muse/2010/12/06/making-a-heatmap-with-r/
graphics.off()  
par(mar=c(10,4,4,10))
colfunc <- colorRampPalette(c("black", "white", "red"))


hm = heatmap.2(heatmapTable, scale = "row", col=brewer.pal(11,"RdBu"),
               cexRow=0.7,cexCol=1,trace="none",srtCol=90, margins=c(7,12), 
               dendrogram = "row", #Colv = F,
               distfun=function(x) dist(x, method="euclidean"), 
               hclustfun=function(x) hclust(x, method= "median"))

#extract the order of the labels
rowOrder = hm$rowInd
colOrder = hm$colInd


#-------- reorder according to the obtained indexes--------
plotTable$Genome = factor(plotTable$Genome, levels = rownames(heatmapTable)[rowOrder])
plotTable$ID = factor(plotTable$ID, levels = colnames(heatmapTable)[colOrder])
# plot as is
p= ggplot(plotTable,aes(x= Organism, y = ID, size = abundance, fill = site))

p = p + geom_point(shape = 21, alpha = 0.8) + theme_bw() +
  scale_size_continuous(range=c(0.2,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Circle plot - unscaled")
p

# Filter out the most abundant
p = plotTable %>% 
  filter(Organism != "Human alphaherpesvirus 1") %>% 
  ggplot(aes(x= Organism, y = ID, size = abundance, fill = site))

p = p + geom_point(shape = 21, alpha = 0.8) + theme_bw() +
  scale_size_continuous(range=c(0.4,15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Circle plot - unscaled / Human alphherpesvirus 1 fltered out")
p
