library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ggtext) 

#----Abundance relativ plot 16S-----

count_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
taxa_tab <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
taxa_tab$Phylum <- str_replace(taxa_tab$Phylum, "Firmicutes", "Bacillota")
taxa_tab$Phylum <- str_replace(taxa_tab$Phylum, "Proteobacteria", "Pseudomonadota")
taxa_tab$Phylum <- str_replace(taxa_tab$Phylum, "Actinobacteriota", "Actinomycetota")
taxa_tab$Phylum <- str_replace(taxa_tab$Phylum, "Cyanobacteria", "Cyanobacteriota")
taxa_tab$Genus <- str_replace(taxa_tab$Genus, "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobium")

taxa_tab <- as.matrix(taxa_tab)
sample_info <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("~/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A" )
sample_info$etat <- c("Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate" )
sample_info$Essai <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM_c", "MSM_G_c", "PCS_c", "PCS_G_c", "R2A_c", "R2A_G_c", "MSM_c", "MSM_G_c", "PCS_c", "PCS_G_c", "R2A_c", "R2A_G_c", "MSM_c", "MSM_G_c", "PCS_c", "PCS_G_c", "R2A_c")
sample_info$Essai <- c("MSM_atm", "MSM_atm_C", "MSM_CO2_C", "MSM_CO2", "MSM_atm", "MSM_atm_C", "MSM_CO2_C", "MSM_CO2","MSM_atm", "MSM_atm_C", "MSM_CO2_C", "MSM_CO2", "PCS_atm", "PCS_atm_C", "PCS_CO2_C", "PCS_CO2", "PCS_atm", "PCS_atm_C", "PCS_CO2_C", "PCS_CO2", "PCS_atm", "PCS_atm_C", "PCS_CO2_C", "PCS_CO2", "R2A_atm", "R2A_atm_C", "R2A_CO2_C", "R2A_CO2", "R2A_atm", "R2A_atm_C", "R2A_CO2_C", "R2A_CO2", "R2A_atm", "R2A_atm_C", "R2A_CO2_C", "R2A_CO2")

phylo_data <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
                       sample_data(sample_info), 
                       tax_table(taxa_tab))

library(dplyr)
library(forcats)
library(ggplot2)

#--Phylum----
ps <- tax_glom(phylo_data, "Phylum")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Essai")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
sample_data_df <- as.data.frame(sample_data(ps2))

sample_data_df$etat <- c("Enrichment", "Substrate","Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate")
sample_data_df$Essai <- c("MSM_atm", "MSM_atm", "MSM_CO2", "MSM_CO2", "PCS_atm", "PCS_atm", "PCS_CO2", "PCS_CO2","R2A_atm", "R2A_atm", "R2A_CO2", "R2A_CO2")
sample_data(ps2) <- sample_data(sample_data_df)


plot_bar(ps2, fill="Phylum")
df <- psmelt(ps2)

top_phyla <- df %>%
  group_by(Essai, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_phyla

top20 <- top_phyla$Phylum[1:22]
df0 <- df %>%
  mutate(Phylum = fct_other(Phylum, top20))

a <- ggplot(df0, aes(Essai, Abundance, fill = Phylum)) +
  geom_col() + facet_grid( ~ etat, scales = "free", space = "free") +
  theme(legend.position = "top") +
  xlab("") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_get()
a + theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), strip.text.x = element_text(size = 11), axis.title = element_text(size = 12), axis.text.x = element_text(size = 12, angle = 90), legend.title = element_text(size=12), legend.text = element_text(size=10))



#--Genus----
ps <- tax_glom(phylo_data, "Genus")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Essai")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
sample_data_df <- as.data.frame(sample_data(ps2))

sample_data_df$etat <- c("Enrichment", "Substrate","Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate")
sample_data_df$Essai <- c("MSM_atm", "MSM_atm", "MSM_CO2", "MSM_CO2", "PCS_atm", "PCS_atm", "PCS_CO2", "PCS_CO2","R2A_atm", "R2A_atm", "R2A_CO2", "R2A_CO2")
sample_data(ps2) <- sample_data(sample_data_df)


plot_bar(ps2, fill="Genus")
df <- psmelt(ps2)

top_phyla <- df %>%
  group_by(Essai, Genus) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_phyla

top20 <- top_phyla$Genus[1:40]
df0 <- df %>%
  mutate(Genus = fct_other(Genus, top20))

a <- ggplot(df0, aes(Essai, Abundance, fill = Genus)) +
  geom_col() + facet_grid( ~ etat, scales = "free", space = "free") +
  theme(legend.position = "top") +
  xlab("") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_get() + scale_color_brewer(palette = "Dark2")
a + theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), strip.text.x = element_text(size = 11), axis.title = element_text(size = 12), axis.text.x = element_text(size = 12, angle = 90), legend.title = element_text(size=12), legend.text = element_text(size=10, face = "italic"))



#----Abundance relativ plot ITS-----

count_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
taxa_tab <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
taxa_tab$Phylum <- str_replace(taxa_tab$Phylum, "Fungi_phy_Incertae_sedis", "Unknown")
taxa_tab <- as.matrix(taxa_tab)

sample_info <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G" )
sample_info$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A")
sample_info$etat <- c("Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Enrichment", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate", "Substrate")
sample_info$Essai <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM_c", "MSM_G_c", "PCS_c", "PCS_G_c", "R2A_c", "R2A_G_c", "MSM_c", "MSM_G_c", "PCS_c", "PCS_G_c", "R2A_c", "R2A_G_c", "MSM_c", "MSM_G_c", "PCS_c", "PCS_G_c", "R2A_c", "R2A_G_c")

taxa_tab_3 <- taxa_tab %>%
  transform(Genus=str_replace(Genus,"g__",""))
taxa_tab_4 <- taxa_tab_3 %>%
  transform(Phylum=str_replace(Phylum,"p__",""))
taxa_tab <- as.matrix(taxa_tab_4)

phylo_data <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), 
                       sample_data(sample_info), 
                       tax_table(taxa_tab))


library(dplyr)
library(forcats)
library(ggplot2)

#--Phylum----
ps <- tax_glom(phylo_data, "Phylum")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Essai")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
sample_data_df <- as.data.frame(sample_data(ps2))

sample_data_df$etat <- c("Enrichment", "Substrate","Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate")
sample_data_df$Essai <- c("MSM_atm", "MSM_atm", "MSM_CO2", "MSM_CO2", "PCS_atm", "PCS_atm", "PCS_CO2", "PCS_CO2","R2A_atm", "R2A_atm", "R2A_CO2", "R2A_CO2")

sample_data(ps2) <- sample_data(sample_data_df)

plot_bar(ps2, fill="Phylum")
df <- psmelt(ps2)

top_phyla <- df %>%
  group_by(Essai, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_phyla

top20 <- top_phyla$Phylum[1:15]
df0 <- df %>%
  mutate(Phylum = fct_other(Phylum, top20))

a <- ggplot(df0, aes(Essai, Abundance, fill = Phylum)) +
  geom_col() + facet_grid( ~ etat, scales = "free", space = "free") +
  theme(legend.position = "top") +
  xlab("") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_get()
a + theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), strip.text.x = element_text(size = 11), axis.title = element_text(size = 12), axis.text.x = element_text(size = 12, angle = 90), legend.title = element_text(size=12), legend.text = element_text(size=10))



#--Genus----
ps <- tax_glom(phylo_data, "Genus")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "Essai")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
sample_data_df <- as.data.frame(sample_data(ps2))

sample_data_df$etat <- c("Enrichment", "Substrate","Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate", "Enrichment", "Substrate")
sample_data_df$Essai <- c("MSM_atm", "MSM_atm", "MSM_CO2", "MSM_CO2", "PCS_atm", "PCS_atm", "PCS_CO2", "PCS_CO2","R2A_atm", "R2A_atm", "R2A_CO2", "R2A_CO2")
sample_data(ps2) <- sample_data(sample_data_df)


plot_bar(ps2, fill="Genus")
df <- psmelt(ps2)

top_phyla <- df %>%
  group_by(Essai, Genus) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_phyla

top20 <- top_phyla$Genus[1:40]
df0 <- df %>%
  mutate(Genus = fct_other(Genus, top20))

a <- ggplot(df0, aes(Essai, Abundance, fill = Genus)) +
  geom_col() + facet_grid( ~ etat, scales = "free", space = "free") +
  theme(legend.position = "top") +
  xlab("") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_get() + scale_color_brewer(palette = "Dark2")
a + theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), strip.text.x = element_text(size = 11), axis.title = element_text(size = 12), axis.text.x = element_text(size = 12, angle = 90), legend.title = element_text(size=12), legend.text = element_text(size=10, face = "italic"))



#----PCOA Bray-Curtis - 16S - Enrichment---------

library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ggExtra)
library(datasets)
library(vegan)
library(phyloseq)
library(dplyr)
library(microViz)
library(phyloseq)
library(tidytree)
library(stringr)

count_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
count_tab_transform <- decostand(count_tab, method = "hellinger")
transform <- t(count_tab_transform)
df2 <- transform[-(19:35),]
df3 <- t(df2)
taxa_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info$CO2 <- str_replace(sample_info$CO2, "20%", "10%")
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A" )
sample_info$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A")
filtered_df <- filter(sample_info, Etat == "Enrichissement") # #Enlever les ligne avec Échantillons de colonisation du sample_info

phylo_data <- phyloseq(otu_table(df3, taxa_are_rows=TRUE), 
                       sample_data(filtered_df), 
                       tax_table(taxa_tab))

bray_dist = phyloseq::distance(phylo_data, method="bray", weighted=F)
ordination = ordinate(phylo_data, method="PCoA", distance=bray_dist)
b <- plot_ordination(phylo_data, ordination, color="Media", shape="CO2") + theme(aspect.ratio=1) + theme_classic() + geom_point(size=5)

b + theme(text = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.title = element_text(size=17), legend.text = element_text(size=17))

adonis2(df2 ~ Milieu*CO2, data=filtered_df, permutations=9999, method="bray") # Permanova sur Milieu et CO2 avec leur interaction

c <- b + annotate(geom = "text", size = 3, x = 0, y = -0.3, label = "PERMANOVA; Media, P = 0.0001")


d <- c + annotate(geom = "text", size = 3, x = 0, y = -0.34, label = "PERMANOVA; CO2, P = 0.0001")
e <- d + annotate(geom = "text", size = 3, x = 0, y = -0.38, label = "PERMANOVA; Media:CO2, P = 0.0001")
e



# PCOA Avec Unifrac - 16S - Enrichment---------

library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ggExtra)
library(datasets)
library(vegan)
library(phyloseq)
library(dplyr)
library(microViz)
library(phyloseq)
library(tidytree)


setwd("~/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce")
fas <- "16S_ASVs.fa" 
seqs <- readDNAStringSet(fas)
seqs
aligned <- AlignSeqs(seqs)
phang.align <- phyDat(as(aligned, "matrix"), type="DNA") # convert to phyDat format
dm <- dist.ml(phang.align) # calculate pairwise distance matrix 
treeNJ <- NJ(dm) # perform neighbor-joining tree method
fit = pml(treeNJ, data=phang.align) # compute intermal max likelihood

count_tab <- as.matrix(read.table("~/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
count_tab_transform <- decostand(count_tab, method = "hellinger")
transform <- t(count_tab_transform)
df2 <- transform[-(19:35),]
df3 <- t(df2)
taxa_tab <- as.matrix(read.table("~/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("~/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info$CO2 <- str_replace(sample_info$CO2, "20%", "10%")
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A" )
sample_info$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A")

filtered_df <- filter(sample_info, Etat == "Enrichissement")
TREE = phy_tree(fit$tree) 

phylo_data <- phyloseq(otu_table(df3, taxa_are_rows=TRUE), 
                       sample_data(filtered_df), 
                       tax_table(taxa_tab), phy_tree(TREE))


wunifrac = phyloseq::distance(phylo_data, method="unifrac", weighted=T)
ordination = ordinate(phylo_data, method="PCoA", distance=wunifrac, weighted=TRUE)
b <- plot_ordination(phylo_data, ordination, color="Media", shape="CO2") + theme(aspect.ratio=1) + theme_classic() + geom_point(size=5)
b
b + theme(text = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.title = element_text(size=17), legend.text = element_text(size=17))


wunifrac = phyloseq::distance(phylo_data, method="unifrac", weighted=T)
adonis2(wunifrac ~ Milieu*CO2, data = filtered_df, permutations=9999)

c <- b + annotate(geom = "text", size = 3, x = 0, y = -0.16, label = "PERMANOVA; Media, P = 0.0001")
d <- c + annotate(geom = "text", size = 3, x = 0, y = -0.17, label = "PERMANOVA; CO2, P = 0.0001")
e <- d + annotate(geom = "text", size = 3, x = 0, y = -0.18, label = "PERMANOVA; Media:CO2, P = 0.0002")
e


# PCOA Avec Bray-Curtis - ITS - Enrichment---------

count_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
count_tab_transform <- decostand(count_tab, method = "hellinger")
transform <- t(count_tab_transform)
df2 <- transform[-(19:36),]
df3 <- t(df2)
taxa_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info$CO2 <- str_replace(sample_info$CO2, "20%", "10%")
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G" )
sample_info$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A")
filtered_df <- filter(sample_info, Etat == "Enrichissement") # #Enlever les ligne avec Échantillons de colonisation du sample_info

phylo_data <- phyloseq(otu_table(df3, taxa_are_rows=TRUE), 
                       sample_data(filtered_df), 
                       tax_table(taxa_tab))

bray_dist = phyloseq::distance(phylo_data, method="bray", weighted=F)
ordination = ordinate(phylo_data, method="PCoA", distance=bray_dist)
b <- plot_ordination(phylo_data, ordination, color="Media", shape="CO2") + theme(aspect.ratio=1) + theme_classic() + geom_point(size=5)
b
b + theme(text = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.title = element_text(size=17), legend.text = element_text(size=17))

adonis2(df2 ~ Milieu*CO2, data=filtered_df, permutations=9999, method="bray") # Permanova sur Milieu et CO2 avec leur interaction


# PCOA Unifrac - its Enrichment---------

setwd("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats")
fas <- "ITS_ASVs.fa"    
seqs <- readDNAStringSet(fas)
seqs
aligned <- AlignSeqs(seqs)
phang.align <- phyDat(as(aligned, "matrix"), type="DNA") # convert to phyDat format
dm <- dist.ml(phang.align) # calculate pairwise distance matrix 
treeNJ <- NJ(dm) # perform neighbor-joining tree method
fit = pml(treeNJ, data=phang.align) # compute intermal max likelihood

count_tab <- as.matrix(read.table("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
count_tab_transform <- decostand(count_tab, method = "hellinger")
transform <- t(count_tab_transform)
df2 <- transform[-(19:36),]
df3 <- t(df2)
taxa_tab <- as.matrix(read.table("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
library(stringr)
sample_info$CO2 <- str_replace(sample_info$CO2, "20%", "10%")
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G" )
sample_info$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A")
filtered_df <- filter(sample_info, Etat == "Enrichissement")
TREE = phy_tree(fit$tree)
count_tab_transform <- decostand(count_tab, method = "hellinger")

phylo_data <- phyloseq(otu_table(df3, taxa_are_rows=TRUE), 
                       sample_data(filtered_df), 
                       tax_table(taxa_tab), phy_tree(TREE))


wunifrac = phyloseq::distance(phylo_data, method="unifrac", weighted=T)
ordination = ordinate(phylo_data, method="PCoA", distance=wunifrac, weighted=TRUE)
b <- plot_ordination(phylo_data, ordination, color="Media", shape="CO2") + theme(aspect.ratio=1) + theme_classic() + geom_point(size=5)
b

b + theme(text = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.title = element_text(size=17), legend.text = element_text(size=17))

wunifrac = phyloseq::distance(phylo_data, method="unifrac", weighted=T)
adonis2(wunifrac ~ Milieu*CO2, data = filtered_df, permutations=9999)
library(pairwiseAdonis)
pairwise.adonis(wunifrac, phyloseq::sample_data(phylo_data)$Milieu)

#----PCOA  Bray-Curtis - 16S - Substrate---------

count_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
count_tab_transform <- decostand(count_tab, method = "hellinger")
transform <- t(count_tab_transform)
df2 <- transform[-(1:18),]
df3 <- t(df2)
taxa_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info$CO2 <- str_replace(sample_info$CO2, "20%", "10%")
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A" )
sample_info$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A")
filtered_df <- filter(sample_info, Etat == "Colonisation") # #Enlever les ligne avec Échantillons de enrichissement du sample_info

phylo_data <- phyloseq(otu_table(df3, taxa_are_rows=TRUE), 
                       sample_data(filtered_df), 
                       tax_table(taxa_tab))

bray_dist = phyloseq::distance(phylo_data, method="bray", weighted=F)
ordination = ordinate(phylo_data, method="PCoA", distance=bray_dist)
c <- plot_ordination(phylo_data, ordination, color="Media", shape="CO2") + theme(aspect.ratio=1) + theme_classic() + geom_point(size=5)
c
c + theme(text = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.title = element_text(size=17), legend.text = element_text(size=17))

adonis2(df2 ~ Milieu*CO2, data=filtered_df, permutations=9999, method="bray") # Permanova sur Milieu et CO2 avec leur interaction


# PCOA Unifrac - 16S - Substrate---------

library(dada2)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ggExtra)
library(datasets)
library(vegan)
library(phyloseq)
library(dplyr)
library(microViz)
library(phyloseq)
library(tidytree)


setwd("~/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce")
fas <- "16S_ASVs.fa" 
seqs <- readDNAStringSet(fas)
seqs
aligned <- AlignSeqs(seqs)
phang.align <- phyDat(as(aligned, "matrix"), type="DNA") # convert to phyDat format
dm <- dist.ml(phang.align) # calculate pairwise distance matrix 
treeNJ <- NJ(dm) # perform neighbor-joining tree method
fit = pml(treeNJ, data=phang.align) # compute intermal max likelihood

count_tab <- as.matrix(read.table("~/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
count_tab_transform <- decostand(count_tab, method = "hellinger")
transform <- t(count_tab_transform)
df2 <- transform[-(1:18),]
df3 <- t(df2)
taxa_tab <- as.matrix(read.table("~/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("~/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info$CO2 <- str_replace(sample_info$CO2, "20%", "10%")
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A" )
sample_info$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A")

filtered_df <- filter(sample_info, Etat == "Colonization")
TREE = phy_tree(fit$tree) 

phylo_data <- phyloseq(otu_table(df3, taxa_are_rows=TRUE), 
                       sample_data(filtered_df), 
                       tax_table(taxa_tab), phy_tree(TREE))


wunifrac = phyloseq::distance(phylo_data, method="unifrac", weighted=T)
ordination = ordinate(phylo_data, method="PCoA", distance=wunifrac, weighted=TRUE)
b <- plot_ordination(phylo_data, ordination, color="Media", shape="CO2") + theme(aspect.ratio=1) + theme_classic() + geom_point(size=5)
b

b + theme(text = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.title = element_text(size=17), legend.text = element_text(size=17))


wunifrac = phyloseq::distance(phylo_data, method="unifrac", weighted=T)
adonis2(wunifrac ~ Milieu*CO2, data = filtered_df, permutations=9999)

# PCOA Bray-Curtis - ITS - substrate---------

count_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
count_tab_transform <- decostand(count_tab, method = "hellinger")
transform <- t(count_tab_transform)
df2 <- transform[-(1:18),]
df3 <- t(df2)
taxa_tab <- as.matrix(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info$CO2 <- str_replace(sample_info$CO2, "20%", "10%")
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G" )
sample_info$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A")
filtered_df <- filter(sample_info, Etat == "Colonisation") # #Enlever les ligne avec Échantillons de colonisation du sample_info

phylo_data <- phyloseq(otu_table(df3, taxa_are_rows=TRUE), 
                       sample_data(filtered_df), 
                       tax_table(taxa_tab))

bray_dist = phyloseq::distance(phylo_data, method="bray", weighted=F)
ordination = ordinate(phylo_data, method="PCoA", distance=bray_dist)
b <- plot_ordination(phylo_data, ordination, color="Media", shape="CO2") + theme(aspect.ratio=1) + theme_classic() + geom_point(size=5)
b
b + theme(text = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.title = element_text(size=17), legend.text = element_text(size=17))

adonis2(df2 ~ Milieu*CO2, data=filtered_df, permutations=9999, method="bray") # Permanova sur Milieu et CO2 avec leur interaction

# PCOA Unifrac - ITS - substrate---------

setwd("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats")
fas <- "ITS_ASVs.fa"    
seqs <- readDNAStringSet(fas)
seqs
aligned <- AlignSeqs(seqs)
phang.align <- phyDat(as(aligned, "matrix"), type="DNA") # convert to phyDat format
dm <- dist.ml(phang.align) # calculate pairwise distance matrix 
treeNJ <- NJ(dm) # perform neighbor-joining tree method
fit = pml(treeNJ, data=phang.align) # compute intermal max likelihood

count_tab <- as.matrix(read.table("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
count_tab_transform <- decostand(count_tab, method = "hellinger")
transform <- t(count_tab_transform)
df2 <- transform[-(1:18),]
df3 <- t(df2)
taxa_tab <- as.matrix(read.table("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/tax_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info <- as.data.frame(read.table("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t"))
sample_info$CO2 <- str_replace(sample_info$CO2, "20%", "10%")
sample_info$Echantillon <- c("MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G", "MSM", "MSM_G", "PCS", "PCS_G", "R2A", "R2A_G" )
sample_info$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A")
filtered_df <- filter(sample_info, Etat == "Colonisation")
TREE = phy_tree(fit$tree)
count_tab_transform <- decostand(count_tab, method = "hellinger")

phylo_data <- phyloseq(otu_table(df3, taxa_are_rows=TRUE), 
                       sample_data(filtered_df), 
                       tax_table(taxa_tab), phy_tree(TREE))


wunifrac = phyloseq::distance(phylo_data, method="unifrac", weighted=T)
ordination = ordinate(phylo_data, method="PCoA", distance=wunifrac, weighted=TRUE)
b <- plot_ordination(phylo_data, ordination, color="Media", shape="CO2") + theme(aspect.ratio=1) + theme_classic() + geom_point(size=5)
b

b + theme(text = element_text(size = 16), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), legend.title = element_text(size=17), legend.text = element_text(size=17))

wunifrac = phyloseq::distance(phylo_data, method="unifrac", weighted=T)
adonis2(wunifrac ~ Milieu*CO2, data = filtered_df, permutations=9999)



# db-rda---------

# https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html#8_Some_Other_Examples
library(vegan)
library(ggplot2)
library(phyloseq) ;packageVersion("phyloseq")
library(Biostrings) ;packageVersion("Biostrings")
library(tidyverse) ;packageVersion("tidyverse")
library(dplyr) ;packageVersion("dplyr")
library(vegan) ;packageVersion("vegan")
library(BiocGenerics)
library(dendextend) ;packageVersion("dendextend")
library(ape)

library(readxl)
library(ggsci)
library(ggrepel)
library(ggforce)
library(concaveman)
library(BiodiversityR)
library(ggordiplots)


setwd("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats")
c <- read.csv("Test_philippe2.csv", sep = ";") 
#view(c)
d <- c[,-(3)]
setwd("~/Desktop/Thèse/Séquençage/Essai_16S/Résultat_avec_espèce")
a <- read.csv("essai_rda_pourcentage_croissance.csv") 
d$Mycelium <- a[ , c(2)]
names(d)[names(d) == "Milieu"] <- "Media"

d$Mycelium  <- decostand(d$Mycelium, method = "standardize")
d$Proportion <- decostand(d$Proportion, method = "standardize")

str(d)


count_tab <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/count_tab_filter.tsv", header=T, row.names=1, check.names=F, sep="\t"))
names(count_tab)[names(count_tab) == "A1"] <- "MSM"
names(count_tab)[names(count_tab) == "A1G"] <- "MSM_G"
names(count_tab)[names(count_tab) == "B1"] <- "PCS"
names(count_tab)[names(count_tab) == "B1G"] <- "PCS_G"
names(count_tab)[names(count_tab) == "C1"] <- "R2A"
names(count_tab)[names(count_tab) == "C1G"] <- "R2A_G"

names(count_tab)[names(count_tab) == "A2"] <- "MSM"
names(count_tab)[names(count_tab) == "A2G"] <- "MSM_G"
names(count_tab)[names(count_tab) == "B2"] <- "PCS"
names(count_tab)[names(count_tab) == "B2G"] <- "PCS_G"
names(count_tab)[names(count_tab) == "C2"] <- "R2A"
names(count_tab)[names(count_tab) == "C2G"] <- "R2A_G"

names(count_tab)[names(count_tab) == "A3"] <- "MSM"
names(count_tab)[names(count_tab) == "A3G"] <- "MSM_G"
names(count_tab)[names(count_tab) == "B3"] <- "PCS"
names(count_tab)[names(count_tab) == "B3G"] <- "PCS_G"
names(count_tab)[names(count_tab) == "C3"] <- "R2A"
names(count_tab)[names(count_tab) == "C3G"] <- "R2A_G"

names(count_tab)[names(count_tab) == "A1C"] <- "MSM"
names(count_tab)[names(count_tab) == "A1CG"] <- "MSM_G"
names(count_tab)[names(count_tab) == "B1C"] <- "PCS"
names(count_tab)[names(count_tab) == "B1CG"] <- "PCS_G"
names(count_tab)[names(count_tab) == "C1C"] <- "R2A"
names(count_tab)[names(count_tab) == "C1CG"] <- "R2A_G"

names(count_tab)[names(count_tab) == "A2C"] <- "MSM"
names(count_tab)[names(count_tab) == "A2CG"] <- "MSM_G"
names(count_tab)[names(count_tab) == "B2C"] <- "PCS"
names(count_tab)[names(count_tab) == "B2CG"] <- "PCS_G"
names(count_tab)[names(count_tab) == "C2C"] <- "R2A"
names(count_tab)[names(count_tab) == "C2CG"] <- "R2A_G"

names(count_tab)[names(count_tab) == "A3C"] <- "MSM"
names(count_tab)[names(count_tab) == "A3CG"] <- "MSM_G"
names(count_tab)[names(count_tab) == "B3C"] <- "PCS"
names(count_tab)[names(count_tab) == "B3CG"] <- "PCS_G"
names(count_tab)[names(count_tab) == "C3C"] <- "R2A"
names(count_tab)[names(count_tab) == "C3G"] <- "R2A_G"

count_tab_2 <- t(count_tab)
colnames(count_tab_2) <- as.character(1:1051)
prefix <- "ASV_"  # Replace with your desired prefix
colnames(count_tab_2) <- paste(prefix, colnames(count_tab_2), sep = "")
count_tab_3 <- as.matrix(count_tab_2)
otu3 <- decostand(count_tab_3, method="hel")
otu4 <- otu3[-(1:18),] #enlève les données enrichissement


a = dbRDA=capscale(otu4 ~ + Mycelium + Proportion + MSM + PCS + R2A + Atm + Atm20, data=d, dist="bray")
adjR2.dbrda <- RsquareAdj (a)$adj.r.squared
adjR2.dbrda # 0.2315348 - adjusted R2 explained by all 6 variables is 23.1%

#plot(a)
anova.cca(dbRDA)

e <- d[,-(1)]
f <- e[,-(7)]

sel.fs <- adespatial::forward.sel(Y = otu4, X = f, adjR2thresh = adjR2.dbrda) 
n.tests <- ncol (f)  # number of tests equals to number of all variables from which is being selected
pval.adj <- p.adjust (sel.fs$pval, method = 'holm', n = n.tests)
sel.fs$pval.adj <- pval.adj
sel.fs
a = dbRDA=capscale(otu4 ~ Mycelium + PCS + Atm, data=f, dist="bray") 

#---Figure

plot1 <- ordiplot(a, scaling = 1, type = "text", choices=c(1,2))

## extract % explained by the first 2 axes

axis.long2 <- axis.long(a, choices=c(1, 2))

sc_si <- scores(a, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(a, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(a, display="bp", choices=c(1, 2), scaling=1)

lol <- as.data.frame(sc_bp)
lol$labels <- c("Mycelium", "PCS", "Atm")

lol1 <- as.data.frame(sc_si)
lol1$Media <- c("MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A", "R2A", "MSM", "MSM", "PCS", "PCS", "R2A")
lol1$CO2 <- c("ATM", "10%", "ATM", "10%", "ATM", "10%", "ATM", "10%", "ATM", "10%", "ATM", "10%", "ATM", "10%", "ATM", "10%", "ATM")

lol2 <- as.data.frame(sc_sp)
lol2$ASV <- NA
lol2$ASV[c(19, 92, 144, 186, 195)] <- c("ASV_19", "ASV_92", "ASV_144", "ASV_186", "ASV_195")

BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 16),
  axis.text = element_text(size = 16, colour = "gray25"),
  axis.title = element_text(size = 16, colour = "gray25"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.key = element_blank())

plotgg2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long2[1, "label"]) +
  ylab(axis.long2[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=lol1,
             aes(x=CAP1, y=CAP2, colour=Media, shape=CO2), 
             size=4) +
  geom_point(data=lol2, 
             aes(x=CAP1, y=CAP2), color = "grey50", shape=16, size=3) +
  geom_text_repel(data=lol2, 
                  aes(x=CAP1, y=CAP2, label=ASV),
                  colour="grey5", size=4) +
  geom_segment(data=lol, 
               aes(x=0, y=0, xend=CAP1*1.7, yend=CAP2*1.7), 
               colour="brown4",  lineend = "round", linejoin = "round", size = 0.3, arrow = arrow(length = unit(0.2, "inches"))) + 
  geom_text_repel(data=lol, 
                  aes(x=CAP1*1.85, y=CAP2*1.85, label=labels),
                  colour="brown4", size=5) +
  BioR.theme +
  ggsci::scale_colour_npg() +
  coord_fixed(ratio=1) 
plotgg2
plotgg2 + expand_limits(x=c(-1,1.5), y=c(-0.5, 1))


 # mycelium area figure + tukey test ---------

library(ggplot2)
library(dplyr)
library(forcats)
library(multcompView)

setwd("~/Desktop/THÈSE/Article 1 - Pleurotus - Brevundimonas/Photo_mycomatériaux_stage")
c <- read.csv("Essai_anova_Mycomatériaux.csv") 
str(c)
c$Essai[c$Essai == 'T'] <- 'Control'
c$couleur <- c("MSM", "MSM", "MSM", "MSM", "MSM", "MSM", "PCS", "PCS", "PCS", "PCS", "PCS", "PCS", "R2A", "R2A", "R2A", "R2A", "R2A", "control", "control", "control")

anova <- aov(data = c, X..de.mycélium ~ Essai)
summary(anova)
TUKEY <- TukeyHSD(anova, conf.level=.95)
TUKEY

plot(TukeyHSD(anova, conf.level=.95), las = 1)

cld <- multcompLetters4(anova, TUKEY)
print(cld)

c$Essai <- c("MSM_atm", "MSM_CO2", "MSM_atm", "MSM_CO2", "MSM_atm", "MSM_CO2", "PCS_atm", "PCS_CO2", "PCS_atm", "PCS_CO2", "PCS_atm", "PCS_CO2", "R2A_atm", "R2A_CO2", "R2A_atm", "R2A_CO2", "R2A_atm", "control", "control", "control")


#-plot


BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 16),
  axis.text = element_text(size = 16, angle = 90, colour = "gray25"),
  axis.title = element_text(size = 16, colour = "gray25"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.key = element_blank())


p <- ggplot(c, aes(x=Essai, y=X..de.mycélium, fill=couleur))+
  geom_boxplot()+ 
  theme_bw() +
  coord_fixed(ratio=0.2) + BioR.theme +
  scale_fill_manual(values = c("MSM" = "indianred1", "PCS" = "green3", "R2A" = "dodgerblue", "Control" = "white"))

p

# Modify legend titles
p2 <- p + xlab("Treatments") + ylab("Percentage of mycelial growth (%)") + annotate(geom = "text", x = 2, y = 45, label = "c") + 
  annotate(geom = "text", x = 3, y = 51, label = "bc") +
  annotate(geom = "text", x = 4, y = 70, label = "ab") +
  annotate(geom = "text", x = 5, y = 59, label = "abc") +
  annotate(geom = "text", x = 6, y = 64, label = "abc") +
  annotate(geom = "text", x = 7, y = 59, label = "abc") +
  annotate(geom = "text", x = 6, y = 40, label = "ANOVA; P = 0.00287", size = 4) + annotate(geom = "text", x = 1, y = 71, label = "a") + BioR.theme + theme(legend.position= "none")

p2

p2 

theme(text = element_text(size = 16))



# Plot enrichiment vs substrate - species richness- Bacteria---------

essai <- as.data.frame(read.table("/Users/marie/Desktop/THÈSE/Séquençage/Essai_16S/Résultat_avec_espèce/alpha_diversity_3_avecfacteurs.tsv", header=T, check.names=F, sep="\t"))
essai
essai$Essai <- c("MSM_atm", "MSM_atm_C", "MSM_CO2_C", "MSM_CO2", "MSM_atm", "MSM_atm_C", "MSM_CO2_C", "MSM_CO2","MSM_atm", "MSM_atm_C", "MSM_CO2_C", "MSM_CO2", "PCS_atm", "PCS_atm_C", "PCS_CO2_C", "PCS_CO2", "PCS_atm", "PCS_atm_C", "PCS_CO2_C", "PCS_CO2", "PCS_atm", "PCS_atm_C", "PCS_CO2_C", "PCS_CO2", "R2A_atm", "R2A_atm_C", "R2A_CO2_C", "R2A_CO2", "R2A_atm", "R2A_atm_C", "R2A_CO2_C", "R2A_CO2", "R2A_atm", "R2A_atm_C", "R2A_CO2")

anova <- aov(data = essai, Species_richness ~ Essai)
summary(anova)
TukeyHSD(anova, conf.level=.95)
plot(TukeyHSD(anova, conf.level=.95), las = 1)


p <- ggplot(essai, aes(Essai, Species_richness, fill=Stade), stat="identity", position="stack") + theme_bw()
p + geom_boxplot() + labs(x = "Treatments", y = "Species richness") + theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), strip.text.x = element_text(size = 16), axis.title = element_text(size = 16), axis.text.x = element_text(size = 16, angle = 90), legend.title = element_text(size=16), legend.text = element_text(size=16))



# Plot enrichiment vs substrate - species richness- Fungi---------

setwd("~/Desktop/THÈSE/Séquençage/Essai_ITS/Résultats/Résultats alpha diversité ITS")
essai <- as.data.frame(read.table("alpha_diversity_3_avecfacteur.tsv", header=T, check.names=F, sep="\t"))
essai
essai$Essai <- c("MSM_atm", "MSM_atm_C", "MSM_CO2_C", "MSM_CO2", "MSM_atm", "MSM_atm_C", "MSM_CO2_C", "MSM_CO2","MSM_atm", "MSM_atm_C", "MSM_CO2_C", "MSM_CO2", "PCS_atm", "PCS_atm_C", "PCS_CO2_C", "PCS_CO2", "PCS_atm", "PCS_atm_C", "PCS_CO2_C", "PCS_CO2", "PCS_atm", "PCS_atm_C", "PCS_CO2_C", "PCS_CO2", "R2A_atm", "R2A_atm_C", "R2A_CO2_C", "R2A_CO2", "R2A_atm", "R2A_atm_C", "R2A_CO2_C", "R2A_CO2", "R2A_atm", "R2A_atm_C", "R2A_CO2_C", "R2A_CO2")

anova <- aov(data = essai, Species_richness ~ Essai)
summary(anova)
TukeyHSD(anova, conf.level=.95)
plot(TukeyHSD(anova, conf.level=.95), las = 1)


p <- ggplot(essai, aes(Essai, Species_richness, fill=Stade), stat="identity", position="stack") + theme_bw()
p + geom_boxplot() + labs(x = "Treatments", y = "Species richness") + theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), strip.text.x = element_text(size = 16), axis.title = element_text(size = 16), axis.text.x = element_text(size = 16, angle = 90), legend.title = element_text(size=16), legend.text = element_text(size=16))



# Bateria abundance tukey test and plot for substrate based confrontation assay
library(ggplot2)
library(dplyr)
library(forcats)
library(multcompView)

setwd("~/Desktop")
J <- read.csv("J0-J5.csv")
str(J)

J0 <- rep("J0", 36)
J5 <- rep("J5", 54)
Jour <- c(J5, J0)
J_jour <- J 
J_jour$Jour <- Jour

anova <- aov(data = J, Concentration..copies.g.of.soil. ~ Traitement)
summary(anova)

TUKEY <- TukeyHSD(anova, conf.level=.95)
TUKEY

cld <- multcompLetters4(anova, TUKEY)
print(cld)

axis.text = element_text(size = 16),
axis.title = element_text(size = 16),
legend.title = element_text(size = 16),
legend.text = element_text(size = 16))

p <- ggplot(J_jour, aes(x=Traitement, y=Concentration..copies.g.of.soil., fill=Jour)) +
  geom_boxplot() + theme_classic()

p

p2 <- p + theme(text = element_text(size = 16),  axis.title = element_text(size = 14),
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 16))
p2 

p3 <- p2 + xlab("Treatments") + ylab("16S rRNA gene copie number per gram of lignocellulosic substrate")
p3
