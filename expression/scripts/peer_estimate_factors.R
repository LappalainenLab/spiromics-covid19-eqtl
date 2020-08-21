#!/usr/bin/env Rscript
#-------------------------------------------------------------------------------------------------
# PEER factors for normalized expression data:
#   * based on inverse transformed data matrix,
#   * only chr1-22
# Author: Silva Kasela
#-------------------------------------------------------------------------------------------------

## PEER tutorial:
## https://github.com/PMBio/peer/wiki/Tutorial

library(here)
library(data.table)
library(peer)

###########################################################################

# Define arguments -----------------------------------

args <- commandArgs(trailingOnly = TRUE)
print(args)
input <- args[1] # path to normalized expression file
cov <- args[2] # "no_cov" "cov_genopc_sex"

cat("Estimating PEER factors", fill = T)

# Load data -------------------------------------------

# Normalized and filtered data matrix - keep genes on chr1-22
data <- fread(input, header = T, sep = "\t", stringsAsFactors = F, data.table = F)
data <- data[!data[,1] %in% "chrX",]
rownames(data) <- data$gene_id
data <- data[, -c(1:4)]
data <- t(data)
cat("Dimensions of the data matrix w/o probes on chrX:", dim(data), fill = T)

# Optional covariates for PEER: -------------------------------
if (cov %in% "cov_genopc_sex") {
  # Could be genotype PC and sex, for example
  cov_data <- read.table("/path/to/additional_covariates.txt", header = T, sep = "\t", stringsAsFactors = F)
}

###########################################################################

# Estimate PEERs ---------------------------------------

if (cov %in% "cov_genopc_sex") {
  cat("PEERs: with covariates", fill = T)
} else {
  cat("PEERs: no covariates", fill = T)
}

# Number of hidden confounders to infer
k <- 50

# Create the model object
model <- PEER()

# Infer k hidden confounders
PEER_setNk(model, k)
PEER_getNk(model)

# set the observed data
PEER_setPhenoMean(model, data)
dim(PEER_getPhenoMean(model))

# set observed covariates
if (cov %in% "cov_genopc_sex") {
  PEER_setCovariates(model, as.matrix(cov_data))
  head(PEER_getCovariates(model))
}

# set max iterations
PEER_setNmax_iterations(model, 10000)

# using the default priors (uninformative priors)
# PEER uses uninformative priors on weight precision and noise precision by default (Alpha a = 0.001, Alpha b = 0.1, Eps a = 0.1, Eps b = 10)

# perform the inference
time <- system.time(PEER_update(model))
print(time)

# Results -------------------------------------------------

X <- PEER_getX(model)  # samples x PEER factors: posterior mean of the inferred confounders (NxK matrix)
dim(X)

# add relevant row/column names
if (cov %in% "cov_genopc_sex") {
  c <- c(colnames(cov_data), paste0("InferredCov", 1:k))
} else {
  c <- paste0("InferredCov", 1:k)
}
colnames(X) <- c
rownames(X) <- rownames(data)

A <- PEER_getAlpha(model)  # PEER factors x 1: precision (inverse variance) of the weights (Kx1 matrix)
dim(A)
colnames(A) <- "Alpha"
A <- as.data.frame(A)
A$Relevance <- 1.0 / A$Alpha

# plot variance of factors - should be similar to elbow plot
pdf(here("peer", "factor_relevance", paste0("peer_factor_relevance.", cov, ".pdf")), width = 9, height = 8)
plot(A$Relevance, xlab = "Factors", ylab = "Factor relevance", main = "Relevance of PEER factors", type = "b")
dev.off()

# Write out results --------------------------------------

cat("PEER: writing results ... ", fill = T)
write.table(t(X), file = here("peer", paste0("peer_factors.", cov, ".txt")), sep = "\t", col.names = T, row.names = T, quote = F)
write.table(A, file = here("peer", "factor_relevance", paste0("peer_factor_relevance.", cov, ".txt")), sep = "\t", col.names = T, row.names = T, quote = F)

cat("Done!", fill = T)
