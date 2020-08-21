#!/bin/bash

# Calculate aFC and 95% confidence intervals for top eQTLs

# Calculate aFC for cis-eQTLs from:
## 1) Gene count data that was normalized with DESeq size factors and log2 transformed, with the aFC arguments `--min_samps 2` and `--min_alleles 1` and including the same covariates that were used in the cis-eQTL mapping;

mkdir -p eqtl/cis/afc/files

# 1) Get a set of top eQTLs
zcat eqtl/cis/spiromics.cis_qtl.txt.gz | cut -f 1,7 | sed 's/\t/:/g' | grep -v 'gene_id:variant_id' | sort -u | sed 's/:/\t/g' | awk 'BEGIN{OFS="\t";FS="\t"} {if(NR==1){print("pid","sid","sid_chr","sid_pos")}; split($2,a,"_");print($1,$2,a[1],a[2])}' > eqtl/cis/afc/files/top_eqtls.txt

# 2) Run aFC
# Chromosome by chromosome
DIR='~/projects/spiromics_eqtl'
for i in {1..22} X; do
    python ~/programs/aFC/aFC.py --vcf ${DIR}/genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic.vcf.gz --pheno ${DIR}/expression/normalized/spiromics.deseq_log2_expression.bed.gz --qtl ${DIR}/eqtl/cis/afc/files/top_eqtls.txt --cov ${DIR}/cov/spiromics.covariates.txt --chr chr${i} --log_xform 1 --log_base 2 --min_samps 2 --min_alleles 1 --o ${DIR}/eqtl/cis/afc/tmp_${i}.txt --boot 100
done

# 3) Merge aFC across chromosomes
cat ${DIR}/eqtl/cis/afc/tmp_1.txt > ${DIR}/eqtl/cis/afc/spiromics.afc_top_eqtl.txt
for i in {2..22} X; do
    cat ${DIR}/eqtl/cis/afc/tmp_${i}.txt | grep -v "pid" >> ${DIR}/eqtl/cis/afc/spiromics.afc_top_eqtl.txt
done
gzip ${DIR}/eqtl/cis/afc/spiromics.afc_top_eqtl.txt

# 4) Clean up
rm ${DIR}/eqtl/cis/afc/tmp_*
