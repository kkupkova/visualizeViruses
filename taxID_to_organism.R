# the script uploads all the tsv file obtained from Pathoscope - just change the path to the folder with files 
# and at theis point the regular expression is set in a way thet it will extract only files which start with
# virus and end with .tsv - just change the pattern to whatever needed

# it extracts all the taxonomic ID numbers and gets only the unique ones - take those and paste to the link provided

# once we have the "dictionary" with the taxonomic ID and the organism name - rewtite the tsv files and save - 
# the new tsv files have now prefix name_

rm(list=ls())
library(data.table)
library(stringr)   

# get the tax IDs in all of the TSV files
tsvFiles = list.files(path = "results_files_sent_to_Kristyna_022519", pattern = "^virus.*.tsv$", full.names = T)

ti = c()
for (i in tsvFiles){
  fileX = read.table(i, skip = 1, sep = "\t", header = T, as.is = T)
  
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






