library(dada2)

setwd("/scratch/marie62/16S")
path <- ("/scratch/marie62/16S")

list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

### si message d'erreur après avoir lancé plotQualityProfile(fnFs[2])
### --> Erreur : BiocParallel errors
### 1 remote errors, element index: 1
### 0 unevaluated and other errors
### first remote error:
### Error in density.default(qscore): 'x' contient une valeur manquante

### Alors regarder si il y des 0-cycle reads in the sample avec la fonction readFastq() du package ShortRead
##library(ShortRead)
##readFastq(fnFs[1])
### Si > readFastq(fnFs)
###class: ShortReadQ
###length: 933277 reads; width: 0..250 cycles

### Alors --> A simple workaround would be to remove all zero-length reads 
###(or perhaps even those below some minimum length cutoff) 
###before doing the quality profile plotting. 
### See the minLen parameters in filterAndTrim for one way to do that.
### because : 
### The underlying methods we use to read in fastq files break when zero-length reads are present.


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs = TRUE) # On Windows set multithread=FALSE


## Si erreur -->> Erreur dans filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(230, 230),  : 
###These are the errors (up to 5) encountered in individual cores...
###Error in (function (fn, fout, maxN = c(0, 0), truncQ = c(2, 2), truncLen = c(0,  : 
###Mismatched forward and reverse sequence files: 61878, 61874.
###Error in (function (fn, fout, maxN = c(0, 0), truncQ = c(2, 2), truncLen = c(0,  : 
###Mismatched forward and reverse sequence files: 100000, 99973
###De plus : Message d'avis :
###Dans mclapply(seq_len(n), do_one, mc.preschedule = mc.preschedule,  :
###les coeurs programmés 1, 2 ont rencontré des erreurs dans le code utilisateur, 
###toutes les valeurs pour ces tâches seront affectées

### Alors dans FilterAndTrim (..., matchIDs=TRUE)


head(out)

#                                                 reads.in reads.out
#  MI.M03737_0016.001.N724---S502.A1_R1.fastq.gz    261878    248988
#  MI.M03737_0016.001.N724---S503.A1G_R1.fastq.gz   539962    472307


# Learn the Error Rates
# Tool to visualize the frequency of error rate as a function of quality score.
# Necessary for the algortith - see the paper
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST = 20)

#plotErrors(errF, nominalQ=TRUE)


# Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding 
# “abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces 
# computation time by eliminating redundant comparisons.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#We are now ready to apply the core sample inference algorithm to the dereplicated data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# dadaFs[[1]]
# dadaRs[[1]]


# Merging paired ends
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# We can now construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file = "Summary.csv")

unite.ref <- "/scratch/marie62/16S/silva_nr99_v138.1_train_set.fa.gz"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
taxa <- addSpecies(taxa, "~/scratch/marie62/16S/silva_species_assignment_v138.1.fa.gz")

asv_seqs_16S <- colnames(seqtab.nochim)
asv_headers_16S <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers_16S[i] <- paste(">ASV_16S", i, sep="_")
}

# fasta:
asv_fasta_16S <- c(rbind(asv_headers_16S, asv_seqs_16S))
write(asv_fasta_16S, "16S_ASVs.fa")

# count table:
asv_tab_16S <- t(seqtab.nochim)
row.names(asv_tab_16S) <- sub(">", "", asv_headers_16S)
write.table(asv_tab_16S, "16S_ASVs_counts.txt", sep="\t", quote=F)

# tax table:
asv_tax_16S <- taxa
row.names(asv_tax_16S) <- sub(">", "", asv_headers_16S)
write.table(asv_tax_16S, "16S_ASVs_taxonomy.txt", sep="\t", quote=F)

ASV.sample <- seqtab.nochim
colnames(ASV.sample) <- sub(">", "", asv_headers_16S)

ASV.sample.sort <- ASV.sample
for (i in ncol(ASV.sample.sort):1) {
  ASV.sample.sort <- ASV.sample.sort[order(-ASV.sample.sort[,i]),]
}

#write.table(ASV.sample.sort, "Seq1B1_ASV_samples.txt", sep="\t", quote=F)  #Change the name !
write.table(ASV.sample.sort, file="Taxonomy_echantillons_16S.txt", quote=F, col.names=NA) #Change the name !



