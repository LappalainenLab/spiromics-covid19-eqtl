#!/usr/bin/env Rscript
#---------------------------------------------------
# PCA plot with the 1000G samples +
# eigenvectors projected onto MESA
# Author: Silva Kasela
#---------------------------------------------------

library(here)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(tidyr)
library(dplyr)

# Set theme
theme_set(theme_classic() +
            theme(plot.title = element_text(size = 14),
                  plot.subtitle = element_text(size = 13),
                  axis.title = element_text(size = 13),
                  axis.title.y = element_text(vjust = 2),
                  axis.title.x = element_text(vjust = -0.75),
                  axis.text = element_text(size = 12),
                  strip.text = element_text(size = 12),
                  legend.title = element_text(size = 13),
                  legend.text = element_text(size = 12)))

# Read in data ------------------------------------
# 1000G samples
ref <- read.table("~/lab/data/1kg/populations.txt", header = T, sep = "\t", stringsAsFactors = F)
x1kg_indiv <- read.table("~/projects_mesa/projects/genotype_data/pca/1kg/gazal_et_al_2019.table_S4.filtered_unrelated_outbred.txt", header = T, sep = "\t", stringsAsFactors = F)

#These populations have been divided into 5 super populations- https://www.internationalgenome.org/category/population/
#AFR, African
#AMR, Ad Mixed American
#EAS, East Asian
#EUR, European
#SAS, South Asian

# SPIROMICS RNA-seq samples samples
linking <- read.table(here("linking_files", "Topmed_rna_seq_link-3.present_in_added.txt"), header = T, sep = "\t", stringsAsFactors = F)
linking <- linking[linking$NWDID_from_master_present_in_freeze9 %in% 1, c("SUBJID", "NWDID_from_master")]

# Eigenvectors from smartpca
evec <- read.table(here("genotype", "pca", "1kg", "merged_1kg_spiromics.vcf.gz.filtered.ld_pruned.pca.evec"), header = F, stringsAsFactors = F)
colnames(evec) <- c("id", paste0("PC", 1:20), "group")
evec$id <- sapply(evec$id, function(x){unlist(strsplit(x, ":"))[1]})

# Define populations
evec$pop <- NA
evec[evec$group %in% "1000G", "pop"] <- with(evec[evec$group %in% "1000G", ], x1kg_indiv[match(id, x1kg_indiv$IID), "POP"])
evec[evec$group %in% "SPIROMICS", "pop"] <- "SPIROMICS"
evec$pop_long <- evec$pop
evec[evec$group %in% "1000G", "pop_long"] <- with(evec[evec$group %in% "1000G", ], ref[match(pop, ref$Population.Code), "Population.Description"])

# Define superpopulations
evec$superpop <- NA
evec[evec$group %in% "1000G", "superpop"] <- with(evec[evec$group %in% "1000G", ], x1kg_indiv[match(id, x1kg_indiv$IID), "SUPER_POP"])
evec[evec$group %in% "SPIROMICS", "superpop"] <- "SPIROMICS"
evec$superpop_long <- ifelse(evec$superpop %in% "AFR", "African",
                             ifelse(evec$superpop %in% "AMR", "Ad Mixed American",
                                    ifelse(evec$superpop %in% "EAS", "East Asian",
                                           ifelse(evec$superpop %in% "EUR", "European",
                                                  ifelse(evec$superpop %in% "SAS", "South Asian",
                                                         ifelse(evec$superpop %in% "SPIROMICS", "SPIROMICS", NA))))))

# Leave out SAS group
evec <- evec[!evec$superpop %in% "SAS", ]

addmargins(table(evec$superpop, useNA = "ifany"))
addmargins(table(evec$pop_long, useNA = "ifany"))

# PCA plots - All samples --------------------------------
col <- brewer.pal(n = 4, name = "Set1")
col <- c(col, "grey50") # yellow color out
evec %>%
  ggplot(aes(x = PC1, y = PC2, col = superpop_long)) +
  labs(title = "PCA plot",
       subtitle = "SPIROMICS samples projected onto PCs from 1000G",
       x = "PC1",
       y = "PC2",
       color = "Population",
       shape = "Group") +
  geom_point(aes(pch = group)) +
  scale_shape_discrete(solid = F) +
  scale_color_manual(values = col)
ggsave(here("genotype", "pca", "fig", "pca_plot.1kg_phase3_and_spiromics.superpopulations.pdf"), width = 5.5, height = 3.8)

# PCA plots - RNA-seq samples only -----------------------------------
spiromics_rnaseq <- evec %>%
  filter(id %in% linking$NWDID_from_master)

ggplot(data = evec[evec$group %in% "1000G",], aes(x = PC1, y = PC2, col = superpop_long)) +
  labs(title = "PCA plot",
       subtitle = "SPIROMICS samples projected onto PCs from 1000G",
       x = "PC1",
       y = "PC2",
       color = "Population",
       shape = "Group") +
  geom_point(alpha = 0.4) +
  scale_color_brewer(palette = "Set1") +
  geom_point(data = spiromics_rnaseq, col = "grey50", size = 1.5)
ggsave(here("genotype", "pca", "fig", "pca_plot.1kg_phase3_and_spiromics_rnaseq.superpopulations.pdf"), width = 5.5, height = 3.8)
