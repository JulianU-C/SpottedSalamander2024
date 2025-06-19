setwd("C:/Users/Owner/OneDrive - Smithsonian Institution/Spotted Salamander 2024/#R-files/spotted 16s analysis/Compare to Osborne24")

library(dada2)
library(DECIPHER)
# IF had multiple runs, read in second dada2 output here of other run to merge multiple runs
st1 <- readRDS("seqtab_run1.rds")
st2 <- readRDS("seqtab_run2.rds")
st3 <- readRDS("spotted_seqtab.rds")
## You get an error message "Duplicated sample names detected in rownames", but this is ok
## dada2 is just letting you know this
seqtab <- mergeSequenceTables(st1, st2, st3, repeats = "sum")
dim(seqtab)

## Distribution of amplicon sizes in bp, if merging multiple runs, call seqtab 'st.all' instead
# remove seqeunces much shorter than desired: result of nonspecific priming
table(nchar(getSequences(seqtab)))
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 355:405]

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=TRUE, verbose=TRUE)

#dim(seqtab.nochim)

## Still retain 93.3% of sequences after chimera removal
#sum(seqtab.nochim)/sum(seqtab)
#rowSums(seqtab.nochim)

## This won't work if you start at line 167, which is fine
#getN <- function(x) sum(getUniques(x))

## Make sure you change 'out1' to whatever output you end up selecting
#track <- cbind(out3, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

#rownames(track) <- sample.names

#head(track)

#write.csv(track, "dada2_output_spotted16S.csv")


## https://benjjneb.github.io/dada2/training.html

taxa <- assignTaxonomy(seqtab.nochim,
                       "C:/Users/Owner/OneDrive - Smithsonian Institution/Spotted Salamander 2024/#R-files/spotted 16s analysis/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
                       , multithread=TRUE)

## this one can give an error message about memory
## Found a work around: https://github.com/benjjneb/dada2/issues/239
## WORKAROUND
#https://github.com/chuvpne/dada2-pipeline/commit/7964cd67da52faadd31f8d93da6385741984360b
chunk.size <- 10000  # size of taxonomy increments
taxonomy.species <- do.call(rbind,
                            lapply(split(c(1:nrow(taxa)), sort(c(1:nrow(taxa))%%ceiling(nrow(taxa)/chunk.size))),
                                   function(x){
                                     return(addSpecies(taxa[x, ], "C:/Users/Owner/OneDrive - Smithsonian Institution/Spotted Salamander 2024/#R-files/spotted 16s analysis/silva_v138.2_assignSpecies.fa.gz"))
                                   }))
# ran from 12:45 - memory limit (somehow worked fine before??)
#taxa <- addSpecies(taxa, "/Users/Julian/Library/CloudStorage/OneDrive-SmithsonianInstitution/Spotted Salamander 2024/#R-files/spotted 16s analysis/silva_v138.2_assignSpecies.fa.gz")

#  inspect the taxonomic assignments:
taxa.print <- taxonomy.species # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


library(phyloseq)
packageVersion("phyloseq") # 1.50.0
## combine feature table and taxonomy table in same order
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps

## rename ASVs to numbers
new.names <- paste0("ASV", seq(ntaxa(ps))) # Define new names ASV1, ASV2, ...
seqs <- taxa_names(ps) # Store sequences
names(seqs) <- new.names # Make map from ASV1 to full sequence
taxa_names(ps) <- new.names # Rename to human-friendly format
head(tax_table(ps))

## convert feature table to matrix
site_species <-as(otu_table(ps), "matrix")

## need to change this to match mapping file later
rownames(site_species)

## transpose to make a species by site matrix
species_site <- t(site_species)

# taxon table 
tax <- as(tax_table(ps), "matrix")

write.csv(species_site, "osborne16S_feature_table.csv")

write.csv(tax, "osborne16S_taxonomy.csv")
write.csv(seqs, "osborne16S_DNAsequences.csv")

#this replaces seqRFLP function that is no longer available
writeFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))  }
  
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

seqsDf <- as.data.frame(seqs)
seqsDf$ASV <- row.names(seqsDf)

seqsDf <- seqsDf[ ,c(2,1)]

colnames(seqsDf) <- c('name','seq')

writeFasta(seqsDf, "osborne16S_DNAsequences.fasta")
