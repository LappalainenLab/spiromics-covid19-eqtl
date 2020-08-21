#!/usr/bin/env Rscript
#----------------------------------------------
# Add GRCh positions to the UKBB GWAS varaints file
# Author: Silva Kasela
#----------------------------------------------

library(here)
library(data.table)

# Read in data
data <- fread(cmd = paste0("gunzip -c ", here("coloc/input/gwas/variants.tsv.bgz")), header = T, sep = "\t", stringsAsFactors = F, data.table = FALSE)
data$hg19 <- paste0("chr", data$chr, "_", data$pos)

input <- fread(here("coloc", "input", "gwas", "ukbb_variants.input.bed"), header = F, sep = "\t", stringsAsFactors = F, data.table = FALSE)
input$id <- paste0(input$V1, "_", input$V3)
stopifnot(length(unique(input$id)) == nrow(input))
output <- fread(here("coloc", "input", "gwas", "ukbb_variants.output.bed"), header = F, sep = "\t", stringsAsFactors = F, data.table = FALSE)
output$id <- paste0(output$V1, "_", output$V3)
addmargins(table(output$V1))
unlifted <- read.table(here("coloc", "input", "gwas", "ukbb_variants.unlifted.bed"), header = F, sep = "\t", stringsAsFactors = F)
unlifted$id <- paste0(unlifted$V1, "_", unlifted$V3)

cat("No unlifted positions: ", nrow(unlifted), fill = T)

# Remove unlifted read
input <- input[!input$id %in% unlifted$id, ]
rm(unlifted)

output$id.hg19 <- input$id
# Keep onlc chr1-22, chrX
output <- output[output$V1 %in% paste0("chr", c(1:22, "X")), ]

# Add GRCh id to the GWAS file
data$GRCh38_id <- output[match(data$hg19, output$id.hg19), "id"]
data$GRCh38_chr <- output[match(data$hg19, output$id.hg19), "V1"]
data$GRCh38_pos <- output[match(data$hg19, output$id.hg19), "V3"]

# Remove positions that do not map to GRCh38 or where new chr is different from the original
cat("Do not map:", sum(is.na(data$GRCh38_id)), fill = T)
data <- data[!is.na(data$GRCh38_id),]
cat("Maps on other chromosome:", sum(paste0("chr", data$chr) != data$GRCh38_chr), fill = TRUE)
data <- data[paste0("chr", data$chr) == data$GRCh38_chr, ]

# Write out
fwrite(data, file = here("coloc", "input", "gwas", "ukbb_variants.grch38.txt.gz"), col.names = T, row.names = F, sep = "\t", quote = F)

cat("Done!", fill = TRUE)
