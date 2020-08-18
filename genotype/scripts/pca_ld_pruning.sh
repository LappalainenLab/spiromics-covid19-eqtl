#!/bin/bash

# LD-pruning of the common autosomal SNPs in not long-range LD regions with PLINK, create PED file and .pedind file for smartpca
## Regions of long-range LD in European populations from Table 1, Price et al. 2008 AJHG
# Preparing data for smartpca

module load plink/1.90-b3.29

PLINK_BINARY="$1"
[[ -z "${PLINK_BINARY}" ]] && (echo "Please provide PLINK_BINARY file" >&2; exit 1)

INDIV="$2"
[[ -z "${INDIV}" ]] && (echo "Please provide list of indivudals to keep" >&2; exit 1)

OUT_DIR="$3"
[[ -z "${OUT_DIR}" ]] && (echo "Please provide output directory file" >&2; exit 1)

# LD-pruning
## Note: using --snps-only option
plink --bfile ${PLINK_BINARY} \
    --geno 0.01 \
    --mind 0.05 \
    --maf 0.05 \
    --keep ${INDIV} \
    --not-chr x \
    --snps-only \
    --exclude range references/regions_with_long_range_ld.hg38.txt \
    --indep-pairwise 200 100 0.1 \
    --out ${OUT_DIR}/ld_pruned_variants

# Save in ped format
plink --bfile ${PLINK_BINARY} \
    --extract ${OUT_DIR}/ld_pruned_variants.prune.in \
    --keep ${INDIV} \
    --recode \
    --out ${OUT_DIR}/ld_pruned
