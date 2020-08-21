#!/usr/bin/env Rscript
#----------------------------------
# Script to get sdY for coloc
# Author: Silva Kasela
#----------------------------------

library(here)
library(data.table)
library(doParallel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Need to specify four arguments")

pheno_file <- args[1]
cov_file <- args[2]
cores <- as.numeric(args[3])
out_file <- args[4]

cat("Read in phenotype file: ", pheno_file, fill = T)
pheno <- fread(here(pheno_file), header = T, sep = "\t", stringsAsFactors = F, data.table = FALSE)
rownames(pheno) <- pheno$gene_id
pheno <- pheno[,-c(1:4)]

cat("Covariates file: ", cov_file, fill = T)
cov <- fread(here(cov_file), header = T, sep = "\t", stringsAsFactors = F, data.table = FALSE)
rownames(cov) <- cov$id
cov$id <- NULL
cov <- t(cov)
stopifnot(rownames(cov) == colnames(pheno))

cat("Regressing the phenotypes on the covariates", fill = T)

run_regress <- function(pheno_id) {
  y <- as.numeric(pheno[pheno_id,])
  m <- lm(y ~ cov)
  sdY <- sd(m$residuals)
  names(sdY) <- pheno_id
  return(sdY)
}

if (cores > 2) {
  # Parallel computing
  num_cores <- cores - 1
  cl <- makeCluster(num_cores, outfile = "")
  registerDoParallel(cl)
  result <- foreach(i = rownames(pheno), .combine = c, .export = "cov") %dopar% run_regress(i)
  stopCluster(cl)
} else{
  result <- foreach(i = rownames(pheno), .combine = c, .export = "cov") %do% run_regress(i)
}

write.table(result, file = here(out_file), col.names = F, row.names = T, sep = "\t", quote = F)

cat("Done!", fill = T)
