#!/bin/bash

# Liftover GRCh37 (hg19) positions to GRCh38

## UKBB GWAS data from the Neale Lab in GRCh37
## Download variants file:
### wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

module load liftover

# 1) Prepare files for liftover
## Provide BED format file (e.g. input.bed)
## BED files are 0 based

echo "Prepare bed file"
zcat variants.tsv.bgz | cut -f 2,3 | grep -v "chr" | sed 's/\t/_/g' | uniq | awk 'BEGIN{OFS="\t";FS="_"} {print("chr"$1,$2-1,$2)}' > ukbb_variants_input.bed

# 2) Run liftover

echo "Run Liftover"
liftOver ukbb_variants_input.bed ~/data/chain_files/hg19ToHg38.over.chain.gz ukbb_variants_output.bed ukbb_variants_unlifted.bed
