library(phyloseq) ;packageVersion("phyloseq")
library(Biostrings) ;packageVersion("Biostrings")
library(ggplot2) ; packageVersion("ggplot2") #est dans phyloseq ?
library(tidyverse) ;packageVersion("tidyverse")
library(dplyr) ;packageVersion("dplyr")
library(vegan) ;packageVersion("vegan")
library(BiocGenerics)
library(dendextend) ;packageVersion("dendextend")

setwd("/Users/marie/Desktop")
getwd()
list.files()

# Load the data --------------------------------------------------------

count_tab <- read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/ITS_ASVs_counts.txt", header=T, row.names=1, check.names=F, sep="\t")
tax_tab <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/ITS_ASVs_taxonomy.txt", header=T, row.names=1, check.names=F, sep="\t"))

# Filter at 0.005 % -------------------------------

# Sum the ASVs ----

Sum_ASV <- cbind(count_tab, Sum = rowSums(count_tab))
Sum_total <- sum(Sum_ASV$Sum)

# Filter the count table ----

Sum_ASV_Filter <- filter(Sum_ASV, Sum_ASV$Sum >= (0.00005*Sum_total))

# Create a new file without the sum column
Count_tab_filter <- Sum_ASV_Filter %>% select(-(Sum))

# Changer les noms des colonnes
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N727---S505.A1C"] <- "A1"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N727---S506.A1GC"] <- "A1G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N727---S507.B1C"] <- "B1"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N727---S510.B1GC"] <- "B1G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N727---S511.C1C"] <- "C1"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N727---S513.C1GC"] <- "C1G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N728---S502.A2C"] <- "A2"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N728---S503.A2GC"] <- "A2G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N728---S505.B2C"] <- "B2"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N728---S506.B2GC"] <- "B2G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N728---S507.C2C"] <- "C2"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N728---S510.C2GC"] <- "C2G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N728---S511.A3C"] <- "A3"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N728---S513.A3GC"] <- "A3G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N729---S502.B3C"] <- "B3"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N729---S503.B3GC"] <- "B3G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N729---S505.C3C"] <- "C3"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N729---S506.C3GC"] <- "C3G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N723---S505.1its"] <- "A1C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N723---S506.2its"] <- "A1CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N723---S507.3its"] <- "B1C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N723---S510.4its"] <- "B1CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N723---S511.5its"] <- "C1C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N723---S513.6its"] <- "C1CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N724---S502.7its"] <- "A2C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N724---S503.8its"] <- "A2CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N724---S505.9its"] <- "B2C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N724---S506.10its"] <- "B2CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N724---S507.11its"] <- "C2C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N724---S510.12its"] <- "C2CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N724---S511.13its"] <- "A3C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N724---S513.14its"] <- "A3CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N726---S502.15its"] <- "B3C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N726---S503.16its"] <- "B3CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N726---S505.17its"] <- "C3C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N726---S506.18its"] <- "C3CG"

# Filter the tax table ----

### The tax table doesn't have the ASV count, therefor, we'll use the count 
#table filtered to keep the taxonomy needed.

tax_tab_filter <- tax_tab[rownames(Count_tab_filter),]


write.table(Count_tab_filter, file="/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/Count_tab_filter.tsv", sep="\t", quote=F, col.names=NA)
write.table(tax_tab_filter, file="/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/tax_tab_filter.tsv", sep="\t", quote=F, col.names=NA)

# 4. Create a data frame with all your information and count --------------


# 4.1. Create a table of information on your samples ---

### The info table will contain your variables according to the name of your sample.
# --> For each, I indicate, the number of tract and which tree. The 2 
#variables need to be of the same size as the subject (= your sample)
# --> Here, my variables are the number of tracts and the tree.

samples.out <- rownames(t(Count_tab_filter))
sample <- sapply(strsplit(samples.out, "D"), `[`, 1)
co2 <- c("ATM", "20%", "ATM","20%","ATM","20%","ATM","20%","ATM","20%","ATM","20%", "ATM", "20%","ATM","20%","ATM","20%","ATM","20%","ATM", "20%","ATM","20%","ATM","20%","ATM" ,"20%","ATM","20%","ATM","20%","ATM","20%","ATM","20%")
milieu <- c("MSM","MSM", "PCS", "PCS","R2A","R2A", "MSM","MSM","PCS","PCS","R2A","R2A","MSM","MSM","PCS","PCS","R2A","R2A","MSM","MSM","PCS","PCS","R2A","R2A","MSM", "MSM", "PCS", "PCS", "R2A", "R2A","MSM","MSM","PCS","PCS","R2A","R2A")
etat <- c("Enrichissement", "Enrichissement", "Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement", "Enrichissement", "Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Colonisation","Colonisation","Colonisation", "Colonisation","Colonisation","Colonisation","Colonisation","Colonisation","Colonisation" ,"Colonisation","Colonisation","Colonisation","Colonisation","Colonisation","Colonisation","Colonisation","Colonisation","Colonisation")
sample_info <- data.frame(Sample=sample, CO2=co2, Milieu=milieu, Etat=etat)
head(sample_info)
rownames(sample_info) <- samples.out

write.table(sample_info, file="/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/sample_info.tsv", sep="\t", quote=F, col.names=NA)
