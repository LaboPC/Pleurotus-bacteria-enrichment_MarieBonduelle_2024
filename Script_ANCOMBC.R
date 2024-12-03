library(BiocManager)
library(rngtools)
library(mia)
library(miaViz)
library(ANCOMBC)
library(phyloseq) ;packageVersion("phyloseq")
library(Biostrings) ;packageVersion("Biostrings")
library(ggplot2) ; packageVersion("ggplot2") #est dans phyloseq ?
library(tidyverse) ;packageVersion("tidyverse")
library(dplyr) ;packageVersion("dplyr")
library(vegan) ;packageVersion("vegan")
library(BiocGenerics)
library(dendextend) ;packageVersion("dendextend")
library(ape)
library(DT)
library(dplyr)
library(SummarizedExperiment)
library(microbiome)
library(ggpubr)


# Script R2A_CO2 + PCS_CO2 vs the other treatment

count_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
transform <- t(count_tab)
df2 <- transform[-(19:36),]
df3 <- t(df2)
taxa_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/Essai ANCOM/tax_tab_filter_essai_ANCOMBC.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
filtered_df <- filter(sample_info, Etat == "Enrichissement") 
filtered_df$Echantillon <- c("B", "B", "B", "A", "B", "A", "B", "B", "B", "A", "B", "A", "B", "B", "B", "A", "B", "A")
filtered_df$Echantillon <- c("A", "A", "A", "B", "A", "B", "A", "A", "A", "B", "A", "B", "A", "A", "A", "B", "A", "B")


phylo_data <- phyloseq(otu_table(df3, taxa_are_rows=TRUE), 
                       sample_data(filtered_df), 
                       tax_table(taxa_tab))

# Aggregate to taxonomy level
phylum_data = aggregate_taxa(phylo_data, "Genus")
# The taxonomy table
tax_mat = as(tax_table(phylum_data), "matrix")

out = ancombc(data = phylum_data, formula = "Echantillon", tax_level = "Genus",
              p_adj_method = "holm", lib_cut = 1000,
              group = "Echantillon", struc_zero = FALSE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = TRUE, prv_cut = 0.01)
res = out$res
res
a <- as.data.frame(res)
a
str(a)

a$diffexpressed <- "NO"
a$diffexpressed[a$lfc.EchantillonB > 0 & a$q_val.EchantillonB < 0.05] <- "UP"
a$diffexpressed[a$lfc.EchantillonB < 0 & a$q_val.EchantillonB < 0.05] <- "DOWN"
a$delabel <- NA
a$delabel[a$diffexpressed != "NO"] <- a$lfc.taxon[a$diffexpressed != "NO"]
b <- subset(a, diffexpressed != "NO")


b$delabel  <- factor(b$delabel , levels = b$delabel[order(b$lfc.EchantillonB)])

library(writexl)

write_xlsx(b, "/Users/marie/Desktop/new_file_b.xlsx")


p <- ggplot(b, aes(x = lfc.EchantillonB, y = delabel, fill = diffexpressed)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(xmin = lfc.EchantillonB - se.EchantillonB, xmax = lfc.EchantillonB + se.EchantillonB), width = 0.2) +# Flip the coordinates to have the taxonomy names on the y-axis  # Set the colors you want
  labs(x = "Effect size (log fold change)", y = "Taxonomic Family", fill = "Significance") +  # Use a minimal theme
  theme(legend.position = "bottom") 
p + ggtitle("Differential Abundance (ANCOM-BC) between PCS and PCS_G") +
  xlab("PCS_G and R2A_G vs the others treamtments effect size (log fold change)") + ylab("Genus") 