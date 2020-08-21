#!/usr/bin/env Rscript
#-------------------------------------------------------
# Process RNA-seq data
#   1) according to the GTEx pipeline using TMM and inverse normal transformation
#   2) DESeq2 size factors for aFC calculations
# Author: Silva Kasela
#------------------------------------------------------

# GTEx pipeline (https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl):
# 1) Read counts are normalized between samples using TMM (Robinson & Oshlack, Genome Biology, 2010)
# 2) Genes are selected based on the following expression thresholds:
#     >=0.1 TPM in >=20% samples AND
#     >=6 reads (unnormalized) in >=20% samples
# 3) Each gene is inverse normal transformed across samples.

# DESeq2 size factors:
# Size factor is a measure of sequence depth, estimated by the median-of-ratios method
# Count values can be brought to a common scale by dividing by the corresponding size factor

library(here)
library(data.table)
library(edgeR)
library(DESeq2)

# Arguments ------------

args <- commandArgs(trailingOnly = TRUE)
print(args)

normalization_method <- args[1] # "tmm" or "deseq2"
tpm_file <- args[2] # path to tmp
counts_file <- args[3] # path to counts
linking_file <- ifelse(args[4] %in% "NA", NA, args[4]) # linking file: expression_id (first column), genotype_id (second column) (if needed)
annot_file <- args[5]
out_file <- args[6]

tpm_threshold <- 0.1
count_threshold <- 6
sample_frac_threshold <- 0.2

# Read in data -----------

# Annotation file
annot <- read.table(annot_file, header = F, sep = "\t", stringsAsFactors = F)
## Keep only genes from chr1-22, chrX
annot <- annot[annot$V1 %in% c(paste0("chr", 1:22), "chrX"), ]
annot <- annot[annot$V3 %in% "gene", ]
annot$gene_id <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub("gene_id ", "", unlist(strsplit(x, ";"))[1]), "[.]"))[1]
})
## Add start (TSS -1) and end (TSS)
## Note: if strand == +, then (start - 1, start)
## Note: if strand == -, then (end -1, end)
annot$start <- ifelse(annot$V7 %in% "+", annot$V4 - 1, annot$V5 - 1)
annot$end <- ifelse(annot$V7 %in% "+", annot$V4, annot$V5)
# sort annot file by chr and start
annot$chr_number <- as.numeric(ifelse(annot$V1 %in% "chrX", 23, sub("chr", "", annot$V1)))
annot <- annot[order(annot$chr_number, annot$start),]

# Linking file
if (!is.na(linking_file)) {
  cat("Using linking file to use new IDs in the output file", fill = TRUE)
  linking <- read.table(linking_file, header = T, sep = "\t", stringsAsFactors = F)
  linking <- linking[linking$NWDID_from_master_present_in_freeze9 %in% 1, c("SUBJID", "NWDID_from_master")]
  # Sort by the second column
  linking <- linking[order(linking[,2]), ]
}

# TPM matrix
tpm <- read.table(tpm_file, header = T, sep = "\t", stringsAsFactors = F)
rownames(tpm) <- tpm$gene_id
tpm$gene_id <- NULL
tpm[1:5, 1:5]

# Counts
counts <- read.table(counts_file, header = T, sep = "\t", stringsAsFactors = F)
counts <- t(counts)
counts[1:5, 1:5]

# Double checks
stopifnot(colnames(tpm) == colnames(counts))
stopifnot(rownames(tpm) %in% rownames(counts))

# Select overlapping gene_ids
gene_id <- sort(intersect(rownames(tpm), rownames(counts)))
length(gene_id)
tpm <- tpm[gene_id,]
counts <- counts[gene_id,]

# Process data --------------------

# Mask lowly expressed genes
ns <- ncol(tpm)
keep_tpm <- which(apply(tpm, 1, function(x) {sum(as.numeric(x) >= tpm_threshold)}) >= sample_frac_threshold*ns)
keep_counts <- which(apply(counts, 1, function(x){sum(as.numeric(x) >= count_threshold)}) >= sample_frac_threshold*ns)
keep <- gene_id[intersect(keep_tpm, keep_counts)]
cat("No. of genes retrieved after filtering: ", length(keep), fill = T)

# Normalize
if (normalization_method == "tmm") {
  cat("Normalize read counts between samples using TMM", fill = T)
  # Discard genes with all-zero counts
  allzero <- apply(counts, 1, function(x){sum(x == 0)}) == ncol(counts)
  d <- DGEList(counts = counts[!allzero,])
  d <- calcNormFactors(d, method = "TMM", logratio_trim = 0.3, sum_trim = 0.05, Acutoff = -1e10)
  head(d$samples)
  # normalized/rescaled CPM (counts per million): counts_df / (lib_size * tmm) * 1e6
  cpm <- cpm(d, normalized.lib.sizes = TRUE, log = FALSE)
  # Expression measurements for each gene transformed to the quantiles of the standard normal distribution
  # Inverse normal transformation to the rows that passed filtering
  cpm <- cpm[keep, ]
  norm <- apply(cpm, 1, function(x){
    qnorm((rank(x, na.last = 'keep') - 0.5) / sum(!is.na(x)))
  })
  norm <- t(norm)
  dim(norm)

} else if (normalization_method == "deseq2") {
  cat("Normalize read counts using DESeq2 sife factors", fill = T)
  colData <- data.frame("ID" = colnames(counts),
                        stringsAsFactors = F)
  dds <- DESeqDataSetFromMatrix(countData = counts[keep, ],
                                colData = colData,
                                design = ~ 1)
  # Estimate size factor
  dds <- estimateSizeFactors(dds)
  print(head(colData(dds)))
  # Normalized counts - count values divided by size factors
  norm <- counts(dds, normalized = TRUE) # same as counts(dds)[1:6,1]/colData(dds)$sizeFactor[1]
  # log2 transformation
  norm <- apply(norm, 1, function(x) (log2(x + 1)))
  norm <- t(norm)
} else {
  stop("Normalization method must be tmm or deseq2")
}

# Use new IDs sorted by the second column, if linking file is provided
if (!is.na(linking_file)) {
  stopifnot(colnames(norm) %in% linking[,1])
  linking <- linking[linking[,1] %in% colnames(norm), ]
  norm <- norm[, linking[,1]]
  colnames(norm) <- linking[,2]
}

# Result in BED format -------------------
# chr    start (TSS -1)  end (TSS)    gene_id            ID1
if (!all(rownames(norm) %in% annot$gene_id)) {
  cat("Don't find match for all the genes from the reference file - ", sum(!rownames(norm) %in% annot$gene_id), "missing", fill = T)
  idx <- which(rownames(norm) %in% annot$gene_id)
  norm <- norm[idx,]
}

stopifnot(rownames(norm) %in% annot$gene_id)
annot <- annot[annot$gene_id %in% rownames(norm), ]
if (length(unique(annot$gene_id)) != nrow(annot)) {stop("Duplicated gene IDs in the annotation file!")}
addmargins(table(annot$V1))

out <- cbind("#chr" = annot$V1,
             "start" = annot$start,
             "end" = annot$end,
             "gene_id" = annot$gene_id,
             norm[annot$gene_id,])

# Write out -------------------
write.table(out, here("expression", "normalized", out_file), col.names = T, row.names = F, sep = "\t", quote = F)

cat("Done!", fill = T)
