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



