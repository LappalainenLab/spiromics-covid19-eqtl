#!/usr/bin/env Rscript
#---------------------------------------------------------------
# K-nearest Neighbors clustering for 1KG and SPIROMICS-RNAseq samples
# Author: Silva Kasela
#--------------------------------------------------------------

library(here)
library(class)

# Read in data ------------------------------------
# SPIROMICS samples
linking <- read.table(here("linking_files", "Topmed_rna_seq_link-3.present_in_added.txt"), header = T, sep = "\t", stringsAsFactors = F)
linking <- linking[linking$NWDID_from_master_present_in_freeze9 %in% 1, c("SUBJID", "NWDID_from_master")]

# 1000G samples
ref <- read.table("~/data/1kg/populations.txt", header = T, sep = "\t", stringsAsFactors = F)
x1kg_indiv <- read.table("~/data/1kg/gazal_et_al_2019.table_S4.filtered_unrelated_outbred.txt", header = T, sep = "\t", stringsAsFactors = F)

# Eigenvectors from smartpca
evec <- read.table(here("genotype", "pca", "1kg", "merged_1kg_spiromics.vcf.gz.filtered.ld_pruned.pca.evec"), header = F, stringsAsFactors = F)
colnames(evec) <- c("id", paste0("PC", 1:20), "group")
evec$id <- sapply(evec$id, function(x){unlist(strsplit(x, ":"))[1]})

# Define superpopulations
evec$superpop <- NA
evec[evec$group %in% "1000G", "superpop"] <- with(evec[evec$group %in% "1000G", ], x1kg_indiv[match(id, x1kg_indiv$IID), "SUPER_POP"])
evec[evec$group %in% "SPIROMICS", "superpop"] <- "SPIROMICS"

# Define populations
evec$pop <- NA
evec[evec$group %in% "1000G", "pop"] <- with(evec[evec$group %in% "1000G", ], x1kg_indiv[match(id, x1kg_indiv$IID), "POP"])
evec[evec$group %in% "SPIROMICS", "pop"] <- "SPIROMICS"

# Leave out SAS group
evec <- evec[!evec$superpop %in% "SAS", ]

# Choose SPIROMICS eQTL samples
evec <- evec[evec$group %in% "1000G" | evec$id %in% linking$NWDID_from_master,]

# KNN ------------------------------------
# Train and test on 1KG dataset, then apply on the SPIROMICS samples
set.seed(4120)

# Normalize data - same scale
normalize <- function(x) {(x - min(x))/(max(x) - min(x))}
evec_norm <- as.data.frame(apply(evec[, 2:21], 2, normalize))
evec_norm_1kg <- evec_norm[evec$group %in% "1000G",]
evec_norm_spiromics <- evec_norm[evec$group %in% "SPIROMICS",]

# Split 1KG samples into training and test dataset, 2/3 and 1/3 of the samples respectively
idx <- sample(c("train", "test"), nrow(evec_norm_1kg), replace = TRUE, prob = c(0.67, 0.33))
train <- evec_norm_1kg[idx == "train", ]
test <- evec_norm_1kg[idx == "test", ]

# True group labels
labels <- evec[evec$group %in% "1000G", "pop"]
train_labels <- labels[idx == "train"]
test_labels <- labels[idx == "test"]

# Find optimal k (optimal = best accuracy)
accuracy <- c()
knn_model <- list()
for (i in 1:10) {
  m <- knn(train = train, test = test, cl = train_labels, k = i)
  knn_model[[i]] <- m
  accuracy[i] <- mean(test_labels == m)
}
k <- which.max(round(accuracy, 2))
cat("Using", k, "nearest neighbors", fill = T)

pdf(here("genotype", "pca", "fig", "knn_prediction_accuracy.pdf"), width = 4, height = 4)
plot(1:10, accuracy, type = "b",
     main = paste0("Choosing k for KNN\n"), xlab = "k",
     col = ifelse(1:10 == k, "indianred", "black"),
     pch = ifelse(1:10 == k, 19, 1))
dev.off()

# Use on SPIROMICS samples
knn_spiromics <- knn(train = train, test = evec_norm_spiromics, cl = train_labels, k = k)
addmargins(table(knn_spiromics))
knn_spiromics <- as.character(knn_spiromics)
addmargins(table(knn_spiromics))

# Write out results
res <- data.frame("NWDID" = evec[evec$group %in% "SPIROMICS", "id"],
                  "inferred_pop" = knn_spiromics,
                  stringsAsFactors = F)
res$inferred_pop_long <- ref[match(res$inferred_pop, ref$Population.Code), "Population.Description"]
addmargins(table(res$inferred_pop_long))
res$inferred_superpop <- ref[match(res$inferred_pop, ref$Population.Code), "Super.Population.Code"]
addmargins(table(res$inferred_superpop))
write.table(res, file = here("genotype", "pca", "1kg", paste0("spiromics.knn_populations_1kg.txt")), col.names = T, row.names = F, sep = "\t", quote = F)
