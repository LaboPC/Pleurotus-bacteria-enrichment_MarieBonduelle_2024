install.packages("phyloseq")
install.packages("Biostrings")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("vegan")
install.packages("dendextend")

library(phyloseq) ;packageVersion("phyloseq")
library(Biostrings) ;packageVersion("Biostrings")
library(ggplot2) ; packageVersion("ggplot2") 
library(tidyverse) ;packageVersion("tidyverse")
library(dplyr) ;packageVersion("dplyr")
library(vegan) ;packageVersion("vegan")
library(BiocGenerics)
library(dendextend) ;packageVersion("dendextend")
library(ape)

setwd("/Users/marie/Desktop")
getwd()
list.files()

# Load the data --------------------------------------------------------

count_tab <- read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/16S_ASVs_counts.txt", header=T, row.names=1, check.names=F, sep="\t")
tax_tab <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/16S_ASVs_taxonomy.txt", header=T, row.names=1, check.names=F, sep="\t"))

# Filter at 0.005 % -------------------------------

# Sum the ASVs ----

Sum_ASV <- cbind(count_tab, Sum = rowSums(count_tab))
Sum_total <- sum(Sum_ASV$Sum)

# Filter the count table ----

Sum_ASV_Filter <- filter(Sum_ASV, Sum_ASV$Sum >= (0.00005*Sum_total))

# Create a new file without the sum column
Count_tab_filter <- Sum_ASV_Filter %>% select(-(Sum))

# Changer les noms des colonnes
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N724---S502.A1"] <- "A1"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N724---S503.A1G"] <- "A1G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N724---S505.B1"] <- "B1"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N724---S506.B1G"] <- "B1G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N724---S507.C1"] <- "C1"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N724---S510.C1G"] <- "C1G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N724---S511.A2"] <- "A2"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N724---S513.A2G"] <- "A2G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N726---S502.B2"] <- "B2"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N726---S503.B2G"] <- "B2G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N726---S505.C2"] <- "C2"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N726---S506.C2G"] <- "C2G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N726---S507.A3"] <- "A3"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N726---S510.A3G"] <- "A3G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N726---S511.B3"] <- "B3"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N726---S513.B3G"] <- "B3G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N727---S502.C3"] <- "C3"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M03737_0016.001.N727---S503.C3G"] <- "C3G"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N721---S502.A16s"] <- "A1C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N721---S503.B16s"] <- "A1CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N721---S505.C16s"] <- "B1C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N721---S506.D16s"] <- "B1CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N721---S507.E16s"] <- "C1C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N721---S510.F16s"] <- "C1CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N721---S511.G16s"] <- "A2C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N721---S513.H16s"] <- "A2CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N722---S502.I16s"] <- "B2C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N722---S503.J16s"] <- "B2CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N722---S505.K16s"] <- "C2C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N722---S506.L16s"] <- "C2CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N722---S507.M16s"] <- "A3C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N722---S510.N16s"] <- "A3CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N722---S511.O16s"] <- "B3C"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N722---S513.P16s"] <- "B3CG"
names(Count_tab_filter)[names(Count_tab_filter) == "MI.M05812_0333.001.N723---S502.Q16s"] <- "C3C"


# Filter the tax table ----

### The tax table doesn't have the ASV count, therefor, we'll use the count 
#table filtered to keep the taxonomy needed.

tax_tab_filter <- tax_tab[rownames(Count_tab_filter),]


write.table(Count_tab_filter, file="/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/count_tab_filter.tsv", sep="\t", quote=F, col.names=NA)
write.table(tax_tab_filter, file="/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/tax_tab_filter.tsv", sep="\t", quote=F, col.names=NA)

# 4. Create a data frame with all your information and count --------------


# 4.1. Create a table of information on your samples ---

### The info table will contain your variables according to the name of your sample.
# --> For each, I indicate, the number of tract and which tree. The 2 
#variables need to be of the same size as the subject (= your sample)
# --> Here, my variables are the number of tracts and the tree.

samples.out <- rownames(t(Count_tab_filter))
sample <- sapply(strsplit(samples.out, "D"), `[`, 1)
co2 <- c("ATM", "20%", "ATM","20%","ATM","20%","ATM","20%","ATM","20%","ATM","20%", "ATM", "20%","ATM","20%","ATM","20%","ATM","20%","ATM", "20%","ATM","20%","ATM","20%","ATM" ,"20%","ATM","20%","ATM","20%","ATM","20%","ATM")
milieu <- c("MSM","MSM", "PCS", "PCS","R2A","R2A", "MSM","MSM","PCS","PCS","R2A","R2A","MSM","MSM","PCS","PCS","R2A","R2A","MSM","MSM","PCS","PCS","R2A","R2A","MSM", "MSM", "PCS", "PCS", "R2A", "R2A","MSM","MSM","PCS","PCS","R2A")
etat <- c("Enrichissement", "Enrichissement", "Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement", "Enrichissement", "Enrichissement","Enrichissement","Enrichissement","Enrichissement","Enrichissement","Colonisation","Colonisation","Colonisation", "Colonisation","Colonisation","Colonisation","Colonisation","Colonisation","Colonisation" ,"Colonisation","Colonisation","Colonisation","Colonisation","Colonisation","Colonisation","Colonisation","Colonisation")
sample_info <- data.frame(Sample=sample, CO2=co2, Milieu=milieu, Etat=etat)
head(sample_info)
rownames(sample_info) <- samples.out

write.table(sample_info, file="/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/sample_info.tsv", sep="\t", quote=F, col.names=NA)
