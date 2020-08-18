#!/usr/bin/env Rscript
#---------------------------------------------------------
# Fianl set of PCS to use
# Author: Silva Kasela
#---------------------------------------------------------

library(here)
library(ggplot2)

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

# Arguments (chosen after running this script)
group <- c("SPIROMICS")
k <- 4

# Read in data ------------------------------------

# SPIROMICS samples
linking <- read.table(here("linking_files", "Topmed_rna_seq_link-3.present_in_added.txt"), header = T, sep = "\t", stringsAsFactors = F)
linking <- linking[linking$NWDID_from_master_present_in_freeze9 %in% 1, c("SUBJID", "NWDID_from_master")]

# Eigenvectors from smartpca
evec <- read.table(here("genotype", "pca", "smartpca.pca.evec"), header = F, stringsAsFactors = F)
colnames(evec) <- c("id", paste0("PC", 1:20), "group")
evec$id <- sapply(evec$id, function(x){unlist(strsplit(x, ":"))[2]})
evec <- evec[evec$id %in% linking$NWDID_from_master,]

# Subpopulations using 1000 Genomes Project samples and k-nearest neighbors clustering
knn_pop <- read.table(here("genotype", "pca", "1kg", "spiromics.knn_populations_1kg.txt"), header = T, sep = "\t", stringsAsFactors = F)
stopifnot(evec$id == knn_pop$NWDID)

# PCs to choose ------------------

# 1) PCs vs KNN-populations - superpopulations
cat("PCs vs KNN-populations (superpopulations)", fill = TRUE)
res <- data.frame(group = "spiromics_rnaseq",
                  pc = 1:20,
                  adjr2 = NA,
                  pval = NA,
                  stringsAsFactors = F)

for (i in 1:nrow(res)) {
  m <- summary(lm(evec[, (i+1)] ~ as.factor(knn_pop$inferred_superpop)))
  fstat <- m$fstatistic
  pval <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  res[i, 3:4] <- c(m$adj.r.squared, pval)
}

cat("Significant PCs (P < 0.05/20)", fill = T)
print(res[res$pval < 0.05/20,])

# Figure
res$adjr2 <- ifelse(res$adjr2 < 0, 0, res$adjr2)
ggplot(res, aes(x = pc, y = adjr2, col = -log10(pval))) +
  labs(title = "Association with superpopulations from 1KG",
       x = "Genotype PC",
       y = "Adj. R2") +
  geom_point() +
  # geom_point(data = res, aes(x = pc, y = -log10(pval)/120), col = "indianred")
  scale_color_binned(type = "viridis") +
  scale_x_continuous(breaks = 1:20) +
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(here("genotype", "pca", "fig", "adjr2.associations_with_superpopulatsions_from_1kg.pdf"), width = 5.5, height = 2.75)

# 2) PCs vs KNN-populations - subpopulations

# Remove subpopulations with just one individual assaigned to it
count <- table(knn_pop$inferred_pop)
count <- count[count == 1]
knn_pop_1 <- knn_pop[!knn_pop$inferred_pop %in% names(count), ]
evec_1 <- evec[match(knn_pop_1$NWDID, evec$id),]

cat("PCs vs KNN-populations (subpopulatiosn)", fill = TRUE)
res <- data.frame(group = "spiromics_rnaseq",
                  pc = 1:20,
                  adjr2 = NA,
                  pval = NA,
                  stringsAsFactors = F)

for (i in 1:nrow(res)) {
  m <- summary(lm(evec[, (i+1)] ~ as.factor(knn_pop$inferred_pop)))
  fstat <- m$fstatistic
  pval <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  res[i, 3:4] <- c(m$adj.r.squared, pval)
}

cat("Significant PCs (P < 0.05/20)", fill = T)
print(res[res$pval < 0.05/20,])

# Figure
res$adjr2 <- ifelse(res$adjr2 < 0, 0, res$adjr2)
ggplot(res, aes(x = pc, y = adjr2, col = -log10(pval))) +
  labs(title = "Association with subpopulations from 1KG",
       x = "Genotype PC",
       y = "Adj. R2") +
  geom_point() +
  #geom_point(data = res, aes(x = pc, y = -log10(pval)/120), col = "indianred")
  scale_color_binned(type = "viridis") +
  scale_x_continuous(breaks = 1:20) +
  scale_y_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(here("genotype", "pca", "fig", "adjr2.associations_with_subpopulatsions_from_1kg.pdf"), width = 5.5, height = 2.75)

# Write out files with correct number of PCs + gender ---------

# Get sex of the samples (imputed sex from plink)
impute_sex <- read.table(here("genotype", "sex_imputed", "sex_imputed.sexcheck"), header = T, stringsAsFactors = F)
impute_sex <- impute_sex[match(evec$id, impute_sex$FID),]
# By default, F estimates smaller than 0.2 yield female calls, and values larger than 0.8 yield male

cat(paste0("Select top ", k, " genotype PCs"), fill = TRUE)
evec <- evec[, c(1:(k + 1))]
colnames(evec)[1] <- "ID"
evec$gender <- ifelse(impute_sex$SNPSEX %in% 1, 1, 0)
cat(nrow(evec), "individuals", fill = T)
write.table(evec, here("genotype", "pca", paste0("spiromics.144_samples.covariates_pcs.txt")), col.names = T, row.names = F, sep = "\t", quote = F)
