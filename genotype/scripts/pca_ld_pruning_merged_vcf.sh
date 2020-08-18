#!/bin/bash

# Merged vcf to binary PLINK format - only autosomal SNPs, CR 99% and MAF 5% across all the samples, exclude regions of long-range LD in European populations (Table 1 from Price et al. 2008 AJHG)
# LD-pruning
# Prepare data for smartpca

module load plink/1.90-b3.29

VCF="$1"
[[ -z "${VCF}" ]] && (echo "Please provide VCF file" >&2; exit 1)
echo ${VCF}

echo "Process merged VCF"
## Note: using --snps-only option because eigenstrat does not allow very long SNP IDs
plink --vcf ${VCF} \
    --make-bed \
    --geno 0.01 \
    --maf 0.05 \
    --not-chr x \
    --snps-only \
    --output-chr chrM \
    --exclude range references/regions_with_long_range_ld.hg38.txt \
    --vcf-half-call missing \
    --out ${VCF}.filtered

echo "LD-pruning"
plink --bfile ${VCF}.filtered \
    --indep-pairwise 200 100 0.1 \
    --out ${VCF}.filtered.ld_pruned_variants

echo "Save in ped format"
plink --bfile ${VCF}.filtered \
    --extract ${VCF}.filtered.ld_pruned_variants.prune.in \
    --recode \
    --out ${VCF}.filtered.ld_pruned

echo "Generate the ‘.pedind’ file for smartpca"
cat ${VCF}.filtered.ld_pruned.ped | cut -d ' ' -f 1-5 | awk 'BEGIN{OFS=" ";FS=" "} {if(substr($1,1,3) != "NWD") print($1,$2,$3,$4,$5,"1000G"); else print($1,$2,$3,$4,$5,"SPIROMICS")}' > ${VCF}.filtered.ld_pruned.pedind
