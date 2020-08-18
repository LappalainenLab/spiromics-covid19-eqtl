#!/usr/bin/env Rscript
#---------------------------------------------------------
# Annotate variants lookup table with rsIDs (b151_GRCh38p7)
#   * if multiple rsIDs per variant, keep the first instance,
#     and save all the rsIDs into another variable
# Author: Silva Kasela
#--------------------------------------------------------

library(here)
library(data.table)
library(seqminer)

# Read in data ---------

lookup <- fread(here("genotype", "freeze9.pass_only.spiromics_144samples.maf05.biallelic.lookup_table.txt.gz"), header = F, sep = "\t", stringsAsFactors = F, data.table = F)
colnames(lookup) <- c("chr", "variant_pos", "variant_id", "ref", "alt")
lookup$rs_id_dbSNP151_GRCh38p7 <- "."
lookup$rs_id_dbSNP151_GRCh38p7.all_ids <- "."

# Process chromosome by chromosome

genofile <- "~/data/dbsnp/b151_GRCh38p7/00-All.biallelic_new_ids.vcf.gz"
chr <- c(1:22, "X")

for (i in chr) {
  cat("Processing chr", i, fill = TRUE)
  idx <- which(lookup$chr %in% paste0("chr", i))

  # process 50,000 lines once
  k1 <- 1
  k2 <- 50000
  while (k1 < length(idx)) {
    k2 <- min(k2, length(idx))
    range <- paste0(i, ":", min(lookup$variant_pos[idx[c(k1:k2)]]), "-", max(lookup$variant_pos[idx[c(k1:k2)]]))
    sel_vcf <- tabix.read.table(tabixFile = genofile, tabixRange = range, stringsAsFactors = FALSE)
    overlap <- intersect(lookup$variant_id[idx[c(k1:k2)]], sel_vcf$ID)
    sel_vcf <- sel_vcf[sel_vcf$ID %in% overlap, ]

    # get rs id from INFO field
    sel_vcf$rs <- sapply(sel_vcf$INFO, function(x){
      rs <- unlist(strsplit(x, ";"))[1]
      rs <- paste0("rs", unlist(strsplit(rs, "="))[2])
      return(rs)
    })
    lookup[match(overlap, lookup$variant_id), "rs_id_dbSNP151_GRCh38p7"] <- sel_vcf[match(overlap, sel_vcf$ID), "rs"]

    # add all the other IDs too
    lookup[match(overlap, lookup$variant_id), "rs_id_dbSNP151_GRCh38p7.all_ids"] <- lookup[match(overlap, lookup$variant_id), "rs_id_dbSNP151_GRCh38p7"]
    count <- table(sel_vcf$ID)
    # if count > 1
    if (sum(count > 1) > 0) {
      ids_to_add <- names(count[count > 1])
      rs <- sapply(ids_to_add, function(x) {
        rows <- which(sel_vcf$ID %in% x)
        paste(sel_vcf$rs[rows], collapse = ",")
      })
      lookup[match(ids_to_add, lookup$variant_id), "rs_id_dbSNP151_GRCh38p7.all_ids"] <- rs
    }
    k1 <- k2 + 1
    k2 <- k2 + 50000
  }
}

# Write out ------------------
fwrite(lookup, file = here("genotype", "freeze9.pass_only.spiromics_144samples.maf05.biallelic.lookup_table.all_ids.txt.gz"), col.names = T, row.names = F, sep = "\t", quote = F)
lookup$rs_id_dbSNP151_GRCh38p7.all_ids <- NULL
fwrite(lookup, file = here("genotype", "freeze9.pass_only.spiromics_144samples.maf05.biallelic.lookup_table.txt.gz"), col.names = T, row.names = F, sep = "\t", quote = F)

cat("Done!", fill = TRUE)
