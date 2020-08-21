#!/usr/bin/env Rscript
#--------------------------------------------------------------------------
# Preparing covariates for eQTL mapping
#   * Optimanl number of PEER factors
#   * Genotype PCs + sex
# Author: Silva Kasela
#--------------------------------------------------------------------------

library(here)

#################################################################

# Define arguments ---------------

args <- commandArgs(trailingOnly = TRUE)
print(args)
k <- args[1]
peer_file <- args[2]
geno_pc_file <- args[3]

cat(paste0("Preparing covariates matrix with ", k, " PEERs, genotype PCs = sex"), fill = T)

# Read in data ---------------

## Genotype PCs + sex
geno_pc <- read.table(geno_pc_file, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
geno_pc <- t(geno_pc)

## PEERs
peer <- read.table(peer_file, header = T, sep = "\t", stringsAsFactors = F)
stopifnot(colnames(peer) == colnames(geno_pc))

# Write out cov file with PEER factors + genotype PCs and sex ------
out <- cbind("id" = c(rownames(geno_pc), rownames(peer)[c(1:k)]),
             rbind(geno_pc, peer[1:k,]))
write.table(out, file = "spiromics_covariates.txt", col.names = T, row.names = F, sep = "\t", quote = F)
