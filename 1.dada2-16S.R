### By: Carly Muletz Wolz

### DADA2 pipeline and creating feature table, taxonomy table and sequence file for later analyses


### INITIAL processing of files in terminal to get this into dada2 format 

## From basespace, the files are downloaded with each sample having a folder 
# and within that folder are the forward and reverse reads

## To make it compatible with dada2 and R we need to copy all the fastq files to one folder 
# Then you can delete the old folder that used to hold the reverse and forward reads. 

## You need to navigate to the project folder (most likely called FASTQ_Generation...) in terminal 
# cd /Users/Carly/Documents/Basespace/Glyphus_16S-427389917/FASTQ_Generation_2024-08-24_07_30_26Z-766968204 
# Move all the files within each folders up one folder to the FASTQ files, which will then hold all
# of the fastq files in one main folder as opposed to per sample folders
# then remove the folders with nothing in them now, and unzip

# mv -v *_L001*/* .     ### move all files within those folders that have L001 in them up one folder (forward and reverse reads) 

# rm -r *_L001*/        ### remove all the empty folders, optional

### YOU MUST UNZIP all of the files

# gunzip *_L001*  

###  INSTALL DADA2 if you don't already have it.
### Follow tutorial on how to install (follow 1. and 2.)
# https://benjjneb.github.io/dada2/dada-installation.html
library(dada2)
packageVersion("dada2") #1.34.0 (1.30.0 on windows)

## if updates required, i always update all and say no to from source

# DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html 
#windows
setwd("C:/BaseSpace/SpottedSals-438710834/FASTQ_Generation_2024-11-23_00_49_30Z-788401622")
#mac desktop
setwd("/Users/Julian/Library/CloudStorage/OneDrive-SmithsonianInstitution/Spotted Salamander 2024/#R-files/FASTQ_Generation_2024-11-23_00_49_30Z-788401622")
# add to path
path <- "/Users/Julian/Library/CloudStorage/OneDrive-SmithsonianInstitution/Spotted Salamander 2024/#R-files/FASTQ_Generation_2024-11-23_00_49_30Z-788401622"
path <- "C:/BaseSpace/SpottedSals-438710834/FASTQ_Generation_2024-11-23_00_49_30Z-788401622"
list.files(path)

#FILTER AND TRIM 

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#INSPECTION OF QUALITY PROFILES (this takes a while!)
plotQualityProfile(fnFs[8:16]) # forward
plotQualityProfile(fnRs[8:16]) # reverse

# Quality looks ok, starts dropping on reverse around 180
# Since need to trim forward primer off 19 bp, will do 275. We want 450 bp merged

#FILTERING AND TRIMMING 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## parameters for filtering data
# Not uncommon to lose 6 to 10k sequences between reads.in and reads.out
# 'maxEE' to 2,5 with low quality samples. This loosens the parameters with the reverse reads
# mutlithread doesn't work in windows, will default to false
#out1 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(275,175),
 #                     maxN=0, maxEE=c(2,5), trimLeft = 19, trimRight = 23, 
  #                    truncQ=2, rm.phix=TRUE,
   #                   compress=TRUE, multithread=TRUE)

#out2 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(275,180), 
 #                     maxN=0, maxEE=c(3,3), trimLeft = 19, trimRight = 23,
  #                    truncQ=2, rm.phix=TRUE,
   #                   compress=TRUE, multithread=TRUE)

out3 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(275,180), 
                      maxN=0, maxEE=c(2,2), trimLeft = 19, trimRight = 23,
                      truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE)

#out4 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(275,180), 
 #                     maxN=0, maxEE=c(2,5), trimLeft = 19, trimRight = 23,
  #                    truncQ=2, rm.phix=TRUE,
   #                   compress=TRUE, multithread=TRUE)

head(out1)
str(out1)

mean(out1[,1])
mean(out1[,2])

mean(out2[,1])
mean(out2[,2])

# {Glyphus} an average pf 85% with this is good. Going with out2, not much difference anyway
# {Spotted} average PF is 81%, try to improve with out2
mean(out1[,2])/mean(out1[,1])

# average PF is 84%, sounds good to me
mean(out2[,2])/mean(out2[,1])

# 78%
mean(out3[, 2]) / mean(out3[, 1])

mean(out4[, 2]) / mean(out4[, 1])
# NEED to run whichever parameter last so that filtFs and filtRs are from this
# GOING with out3
errF <- learnErrors(filtFs, multithread=TRUE)
# 100256256 total bases in 391626 reads from 34 samples will be used for learning the error rates
errR <- learnErrors(filtRs, multithread=TRUE)
# 93214170 total bases in 578970 reads from 49 samples will be used for learning the error rates.
#errors look reasonable I think, negative assoc. error freq vs. q-score
plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]


mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

setwd("/Users/Julian/Library/CloudStorage/OneDrive-SmithsonianInstitution/Spotted Salamander 2024/#R-files/spotted 16s analysis")
seqtabspot <- "spotted_seqtab.rds"
## Used out2 parameters
saveRDS(seqtab, "spotted_seqtab.rds")
#can be used to merge with a second run if two runs were necessary, if not, at least saved at this point

library(DECIPHER)
# IF had multiple runs, read in second dada2 output here of other run to merge multiple runs
#st1 <- readRDS("XengK2022exp2_seqtab.rds")
#st2 <- readRDS("XengK2022exp2run2_seqtab.rds")
## You get an error message "Duplicated sample names detected in rownames", but this is ok
## dada2 is just letting you know this
#st.all <- mergeSequenceTables(st1, st2, repeats = "sum")
#dim(st.all)

## Distribution of amplicon sizes in bp, if merging multiple runs, call seqtab 'st.all' instead
# remove seqeunces much shorter than desired: result of nonspecific priming
table(nchar(getSequences(seqtab)))
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 355:405]

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

## Still retain 93.3% of sequences after chimera removal
sum(seqtab.nochim)/sum(seqtab)
rowSums(seqtab.nochim)

## This won't work if you start at line 167, which is fine
getN <- function(x) sum(getUniques(x))

## Make sure you change 'out1' to whatever output you end up selecting
track <- cbind(out3, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)

write.csv(track, "dada2_output_spotted16S.csv")


## https://benjjneb.github.io/dada2/training.html

taxa <- assignTaxonomy(seqtab.nochim,
                       "/Users/Julian/Library/CloudStorage/OneDrive-SmithsonianInstitution/Spotted Salamander 2024/#R-files/spotted 16s analysis/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
                         , multithread=TRUE)

## this one can give an error message about memory
## Found a work around: https://github.com/benjjneb/dada2/issues/239
## WORKAROUND
#https://github.com/chuvpne/dada2-pipeline/commit/7964cd67da52faadd31f8d93da6385741984360b
chunk.size <- 10000  # size of taxonomy increments
taxonomy.species <- do.call(rbind,
                            lapply(split(c(1:nrow(taxa)), sort(c(1:nrow(taxa))%%ceiling(nrow(taxa)/chunk.size))),
                                   function(x){
                                     return(addSpecies(taxa[x, ], "/Users/Julian/Library/CloudStorage/OneDrive-SmithsonianInstitution/Spotted Salamander 2024/#R-files/spotted 16s analysis/silva_v138.2_assignSpecies.fa.gz"))
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

write.csv(species_site, "spotted16S_feature_table.csv")

write.csv(tax, "spotted16S_taxonomy.csv")
write.csv(seqs, "spotted16S_DNAsequences.csv")

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

writeFasta(seqsDf, "spotted16S_DNAsequences.fasta")










