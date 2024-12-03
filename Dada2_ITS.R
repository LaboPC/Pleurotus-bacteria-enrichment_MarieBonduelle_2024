library(dada2)

setwd("/scratch/marie62/ITS")
path <- ("/scratch/marie62/ITS")

list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

plotQualityProfile(fnFs[1])
#plotQualityProfile(fnRs[1])

### The plotted lines show positional summary statistics: 
### green is the mean
### orange is the median
### the dashed orange lines are the 25th and 75th quantiles.

### If the sequences vary in length, a red line will be plotted 
### showing the percentage of reads that extend to at least that position.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))


### Note that, unlike in the 16S Tutorial Workflow, we will not be truncating the reads 
### to a fixed length, as the ITS region has significant biological length variation 
### that is lost by such an appraoch.

#For this dataset, we will use standard filtering paraments: maxN=0 (DADA2 requires sequences contain no Ns)
#truncQ = 2, rm.phix = TRUE and maxEE=2. 
#The maxEE parameter sets the maximum number of “expected errors” allowed in a read
#which is a better filter than simply averaging quality scores. 
#Note: We enforce a minLen here, to get rid of spurious very low-length sequences. 
#This was not needed in the 16S Tutorial Workflow because truncLen already served that purpose.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)

###---> pour truncLen=c(230,230)
###                                                  reads.in reads.out
###MI.M03737_0016.001.N727---S505.A1C_R1.fastq.gz    223891       856
###MI.M03737_0016.001.N727---S506.A1GC_R1.fastq.gz   389911      1031
###MI.M03737_0016.001.N727---S507.B1C_R1.fastq.gz    398205      1082

###---> pour truncLen=c(200,200)
###                                                  reads.in reads.out
###MI.M03737_0016.001.N727---S505.A1C_R1.fastq.gz    223891     83691
###MI.M03737_0016.001.N727---S506.A1GC_R1.fastq.gz   389911    117659
###MI.M03737_0016.001.N727---S507.B1C_R1.fastq.gz    398205    140503


# Learn the Error Rates
# Tool to visualize the frequency of error rate as a function of quality score.
# Necessary for the algortith - see the paper
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST = 20)
#plotErrors(errF, nominalQ=TRUE)

##**What to do with it :** 
## - These plots show the error rates for each possible transition (A->C, A->G, etc.)  

##- Legend : 
##  --> Points = observed error rates for each consensus quality score 
##--> Black lines = estimate error after convergence of the matrice-learning algorithm 
##--> Red lines = shows the error rates expected under the nominal def of the Q-score  

##- Basically, as long as points and black line are a fit + error rates drop, it's good ! 


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
dadaFs[[1]]
dadaRs[[1]]


# Merging paired ends

####- Mergeing of the foward and reverse reads together to obtain a full denoised sequences 
####--> Merging is performed by aligning the denoised F and R reads and construct a merged "contig" sequences 
###- Merdge only if at least 12 bases overlap and are identical to each other 
###in the overlap region That can be changed via function argument if wanted

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
###                                     input filtered denoisedF denoisedR merged nonchim
###MI.M03737_0016.001.N727---S505.A1C  223891    83691     83135     83165  26101   25845
###MI.M03737_0016.001.N727---S506.A1GC 389911   117659    117232    116863  46178   45897
###MI.M03737_0016.001.N727---S507.B1C  398205   140503    140229    139410  48288   47819

unite.ref <- "/scratch/marie62/ITS/sh_general_release_dynamic_29.11.2022.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

asv_seqs_ITS <- colnames(seqtab.nochim)
asv_headers_ITS <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers_ITS[i] <- paste(">ASV_ITS", i, sep="_")
}

# fasta:
asv_fasta_ITS <- c(rbind(asv_headers_ITS, asv_seqs_ITS))
write(asv_fasta_ITS, "ITS_ASVs.fa")

# count table:
asv_tab_ITS <- t(seqtab.nochim)
row.names(asv_tab_ITS) <- sub(">", "", asv_headers_ITS)
write.table(asv_tab_ITS, "ITS_ASVs_counts.txt", sep="\t", quote=F)

# tax table:
asv_tax_ITS <- taxa
row.names(asv_tax_ITS) <- sub(">", "", asv_headers_ITS)
write.table(asv_tax_ITS, "ITS_ASVs_taxonomy.txt", sep="\t", quote=F)

ASV.sample <- seqtab.nochim
colnames(ASV.sample) <- sub(">", "", asv_headers_ITS)

ASV.sample.sort <- ASV.sample
for (i in ncol(ASV.sample.sort):1) {
  ASV.sample.sort <- ASV.sample.sort[order(-ASV.sample.sort[,i]),]
}

#write.table(ASV.sample.sort, "Seq1B1_ASV_samples.txt", sep="\t", quote=F)  #Change the name !
write.table(ASV.sample.sort, file="Taxonomy_echantillons_ITS.txt", quote=F, col.names=NA) #Change the name !



