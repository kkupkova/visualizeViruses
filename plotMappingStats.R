# get number of mapped reads
rm(list=ls())
library(ggplot2)
library(GGally)
library(cowplot)
library(ggpubr)

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
for (i in 1:length(tsvFiles)){
  
  
  if(i == 1){
    
    #upload table - select interesting columns - give them ID - merge by Genome and Organism
    fileTable = fread(tsvFiles[i], nrows = 1)
    fileID =  unlist(regmatches(tsvFiles[i], gregexpr("[[:digit:]]+", tsvFiles[i])))
    
  } else {
    fileTable = rbind(fileTable, fread(tsvFiles[i], nrows = 1))
    fileID = c(fileID, unlist(regmatches(tsvFiles[i], gregexpr("[[:digit:]]+", tsvFiles[i]))))
  }
}


# put everythong into data frame
statsTable = data.frame(mappedReads = fileTable$V2, noOfGenomes = fileTable$V4, ID = factor(fileID))
statsTable = merge(statsTable, metadata, by = "ID") %>% 
  mutate_if(is.character, as.factor)

# plot number of mapped reads 
a = ggplot(statsTable, aes(reorder(ID, mappedReads), mappedReads, fill = site))+  geom_bar(stat = "identity", alpha = 0.8) + 
  xlab("ID") + theme(axis.text.x = element_text(angle = 90))

b = ggplot(statsTable, aes(reorder(ID, mappedReads), mappedReads, fill = laz)) + geom_bar(stat = "identity", alpha = 0.8)  + 
  xlab("ID") + theme(axis.text.x = element_text(angle = 90))

# plot number of genomes - viral genomes it mapped to
c = ggplot(statsTable, aes(reorder(ID, noOfGenomes), noOfGenomes, fill = site)) + geom_bar(stat = "identity", alpha = 0.8)  + 
  xlab("ID") + theme(axis.text.x = element_text(angle = 90))

d = ggplot(statsTable, aes(reorder(ID, noOfGenomes), noOfGenomes, fill = laz)) + geom_bar(stat = "identity", alpha = 0.8)  + 
  xlab("ID") + theme(axis.text.x = element_text(angle = 90))

# is the data normal?
qqline(statsTable$noOfGenomes)

# plot number of mapped reads and number of genomes as boxplots and get p-vals
wilkG = ggplot(statsTable, aes(site, noOfGenomes, fill = site)) + geom_boxplot()+ geom_point() + stat_compare_means(vjust = -0.5)
wilkG

wilkR = ggplot(statsTable, aes(site, mappedReads, fill = site)) + geom_boxplot() + geom_point() + stat_compare_means(vjust = -0.5)
wilkR

tG = ggplot(statsTable, aes(site, noOfGenomes, fill = site)) + geom_boxplot()+ geom_point() + stat_compare_means(method = "t.test", vjust = -0.5)
tG

tR = ggplot(statsTable, aes(site, mappedReads, fill = site)) + geom_boxplot()+ geom_point() + stat_compare_means(method = "t.test", vjust = -0.5)
tR

# plot all into one grid
plot_grid(a,c,b,d, labels = c("A", "B", "C", "D"))
plot_grid(tG, tR)


p = ggpairs(statsTable, aes(color = site, alpha = 0.5), columns = c(2,3,8,6) )
p
