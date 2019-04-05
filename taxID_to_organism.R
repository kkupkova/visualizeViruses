rm(list=ls())
library(data.table)
library(stringr)   

# get the tax IDs in all of the TSV files
tsvFiles = list.files(path = "results_files_sent_to_Kristyna_022519/", pattern = "^virus.*.tsv$", full.names = T)

ti = c()
for (i in tsvFiles){
  fileX = read.table(i, skip = 1, sep = "\t", header = T)
  
  helpVec = unlist(strsplit(fileX$Genome, split = "|", fixed = T))
  
  helpVec = helpVec[seq(2, length(helpVec), by = 2)]
  
  ti = c(ti, helpVec)
  
}

taxID = unique(ti)

# write the taxIDs into a text file - upload it to following link and hit "Save in file"
write.table(taxID, file = "taxIDs.txt", sep = "\n", col.names = F, row.names = F, quote = F)


# paste the text file to following url:
# https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi


# upload the generated text file = dictionary and add the extra column to the originals
taxIDdictionary = fread("tax_report.txt")


taxIDdictionary = data.frame(ti = taxIDdictionary$`primary taxid`, organism = taxIDdictionary$taxname)
rownames(taxIDdictionary) = paste0("ti|",taxIDdictionary$ti)


for (i in tsvFiles){
  
  con <- file(i,"r")
  FXheader = readLines(con,n=1)
  fileX = read.table(i, skip = 1, sep = "\t", header = T)
  close(con)
  
  rownames(fileX) = fileX$Genome
  
  
  fileX$Organism = taxIDdictionary[rownames(fileX), 2]
  
  colnames(fileX) <- str_replace_all(colnames(fileX), "[:punct:]", " ")
  
  cat(FXheader, "\n",file=paste0("name_", basename(i)))
  write.table(fileX, paste0("name_", basename(i)),sep="\t",append=TRUE, col.names=T, row.names = F, quote = F)
}






