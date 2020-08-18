# eQTL mapping in SPIROMICS bronchial samples for genes related to COVID-19 response: processing genotype data

This repository documents the workflow and code used to process genotype data for eQTL mapping in SPIROMICS.

Using the following software:

* bcftools/1.9
* plink/1.90-b3.29
* eigensoft/6.1.3
* R/3.4.1

## a) Pull out genotypes, sort out PCs, MAF filter

Using genotype data from the [TOPMed](https://doi.org/10.1101/563866) Freeze 9 (_freeze.9.readme.txt_):

> Freeze 9: After variant filtering, this produces a dataset with approximately 781 M SNP records and 62 M short indel records which pass variant QC, and an additional 104 M records where support vector machine (SVM) filtering suggests that genotypes are less reliable. /.../ Multiallelic variants are represented as biallelic records on consecutive rows.

* minDP10 = genotypes for all individuals at all sites, with genotypes shown as missing ("./.") when an individual has less than 10x sequencing depth at this site.
* Spiromics - 2,710 samples from SPIROMICS (study PI Meyers)

### __Process unphased genotypes data__

* Select only SNPs and Indels that PASS all the filters & select only 2,710 SPIROMICS samples
* Select common variants with MAF >= 1%
* Define variant IDs - chrXX_pos_ref_alt_b38

```bash
scripts/genotype.process_chr.sh
```

* Merge individual files by chromosome into one

```bash
## Create file list to be merged
echo "genotype/tmp_chr1.freeze9.pass_only.spiromics_2710samples.maf01.bcf" > genotype/file_list.txt
for x in {2..22} X; do
    echo "genotype/tmp_chr${x}.freeze9.pass_only.spiromics_2710samples.maf01.bcf" >> genotype/file_list.txt
done

# Merge chromosomes into one vcf
bcftools concat --threads 12 -Ob --file-list genotype/file_list.txt -o genotype/freeze9.pass_only.spiromics_2710samples.maf01.bcf
```

* Keep only biallelic sites
	* "True multiallelic sites are not observed very frequently unless you look at very large cohorts, so they are often taken as a sign of a noisy region where artifacts are likely." [GATK Glossary]("https://gatk.broadinstitute.org/hc/en-us/articles/360035890771-Biallelic-vs-Multiallelic-sites)

```bash
# Firstly, join biallelic sites into multiallelic records
# Secondly, include only these variants that have exactly two alleles listed in REF and ALT columns
bcftools norm --threads 10 -Ou --multiallelics + genotype/freeze9.pass_only.spiromics_2710samples.maf01.bcf | bcftools view --threads 10 -Ob -m2 -M2 -o genotype/freeze9.pass_only.spiromics_2710samples.maf01.biallelic.bcf
```

* Convert to VCF and index with tabix

```bash
bcftools view --threads 20 -Oz -o genotype/freeze9.pass_only.spiromics_2710samples.maf01.biallelic.vcf.gz genotype/freeze9.pass_only.spiromics_2710samples.maf01.biallelic.bcf && bcftools index --tbi genotype/freeze9.pass_only.spiromics_2710samples.maf01.biallelic.vcf.gz
```

* Convert VCF to Plink binary files

```bash
PLINK_BINARY='genotype/freeze9.pass_only.spiromics_2710samples.maf01.biallelic'
plink --make-bed --output-chr chrM --vcf ${PLINK_BINARY}.vcf.gz --out ${PLINK_BINARY}

# Delete tmp files
rm genotype/tmp_chr*
```

* VCF file for individuals with RNA-seq data + filtering of MAF >= 0.05 - __144 indiv__

```bash
# Hack to exclude monomorphic variants: AF[0]==0.5 && AC[0]==N && COUNT(GT="het")=N (happens when all samples are hets)
bcftools view -Ou --samples-file linking_files/Topmed_rna_seq_link-3.NWDID_present_in_freeze9_144indiv.txt genotype/freeze9.pass_only.spiromics_2710samples.maf01.biallelic.bcf | bcftools view -Ou -i 'MAF>=0.05' | bcftools view -e 'AF[0]==0.5 && AC[0]==144 && COUNT(GT=\"het\")=144' -Ob -o genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic.bcf

# To VCF and index
bcftools view --threads 20 -Oz -o genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic.vcf.gz genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic.bcf && bcftools index --tbi genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic.vcf.gz

# To plink binary (genotype CR >= 0.9)
PLINK_BINARY='genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic'
plink --make-bed --output-chr chrM --geno 0.1 --vcf ${PLINK_BINARY}.vcf.gz --out ${PLINK_BINARY}
```

* Create lookup table for WGS variants

Two files: 1) if variant is present multiple times, then take first rsID from the annotation file (ID from the oldest dbSNP version, usually); 2) list all rsIDs

```bash
# Lookup table with the following columns: chr, variant_pos, variant_id, ref, alt, rs_id_dbSNP151_GRCh38p7
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic.bcf > genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic.lookup_table.txt
gzip genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic.lookup_table.txt

# Annotate lookup table
scripts/create_lookup_table.R
```

### __Principal Component Analysis (PCA)__

Perform PCA on a set of LD-independent autosomal SNPs from not long-range LD regions with a call rate >= 99% and MAF >= 0.05 using [EIGENSTRAT](https://github.com/argriffing/eigensoft/blob/master/EIGENSTRAT).

* Perform LD pruning on the post-sample and variant QC WGS VCF of unrelated individuals

```bash
PLINK_BINARY='genotype/freeze9.pass_only.spiromics_2710samples.maf01.biallelic'
INDIV='linking_files/master_with_resends_and_topmedids.cleaned_formatted.unrelated_2684indiv_plink.txt'
OUT_DIR='genotype/pca'
mkdir -p ${OUT_DIR}/fig

# Run LD pruning of autosomal common SNPs with plink,
## Keep only unrelated individuals
## SNPs with a 99% genotyping rate (1% missing), MAF >= 0.05
## Individuals with 95% genotyping rate (5% missing)
scripts/pca_ld_pruning.sh ${PLINK_BINARY} ${INDIV} ${OUT_DIR}
```

* Check identity-by-decent (kinship >= 1/32 = 0.03125)

```bash
OUT_DIR='genotype/pca'
plink --file ${OUT_DIR}/ld_pruned --genome --min 0.03 --out ${OUT_DIR}/ld_pruned
```

* Run `smartpca` to calculate PCs for a given race/ethnic group

```bash
OUT_DIR='genotype/pca'

# Generate the ‘.pedind’ file for smartpca
cat ${OUT_DIR}/ld_pruned.ped | cut -d ' ' -f 1-6 | awk -v group='spiromics' 'BEGIN{OFS=" ";FS=" "} {print($1,$2,$3,$4,$5,group)}' > ${OUT_DIR}/ld_pruned.pedind

# Run smartpca from eigenstrat (same options as in GTEx v8)
smartpca.perl -i ${OUT_DIR}/ld_pruned.ped -a ${OUT_DIR}/ld_pruned.map -b ${OUT_DIR}/ld_pruned.pedind -k 20 -m 0 -o ${OUT_DIR}/smartpca.pca -e ${OUT_DIR}/smartpca.eval -p ${OUT_DIR}/smartpca.plot -l ${OUT_DIR}/smartpca.log
```

* PCA plot

```bash
mkdir -p 'genotype/pca/fig'
PREFIX='SPIROMICS'
HIGHLIGHT_SAMPLES='linking_files/Topmed_rna_seq_link-3.present_in_added.txt'

scripts/pca_plot.R ${PREFIX} ${HIGHLIGHT_SAMPLES}
```

### __Project SPIROMICS samples onto 1KG populations__

To understand the clusters on PCA plots, and to choose the number of PCs to be used in eQTL mapping

Using processed 1KG VCF file:

* Common SNPs and MAF > 5% ('MAF[0]>0.05 & TYPE=\"snp\"')
* Selected only unrelated and outbread individuals from Gazal et al. 2015

Steps for projecting SPIROMICS individuals onto 1KG populations:

* Merge 1000G and SPIROMICS VCF files
* Merged VCF to binary PLINK format and additional filtering before LD-pruning - only autosomal SNPs, CR 99% and MAF 5% across all the samples, exclude regions of long-range LD (Table 1 of Price et al. 2008 AJHG)
* LD-pruning to create list of SNPs to be used in smartpca, output in ped format for smartpca
* Estimate PCs from 1000G samples, and project SPIROMICS samples onto those eigenvectors

```bash
# Merge 1KG VCF with SPIROMICS VCF file
X1KG_VCF='~/data/1kg/merged_genotypes.1kg_phase3_grch38.vcf.gz'
SPIROMICS_VCF='genotype/freeze9.pass_only.spiromics_2710samples.maf01.biallelic.vcf.gz'
mkdir -p genotype/pca/1kg

bcftools merge --threads 16 -m id ${X1KG_VCF} ${SPIROMICS_VCF} -Oz -o genotype/pca/1kg/merged_1kg_spiromics.vcf.gz

# Process merged VCF for smartpca
VCF='genotype/pca/1kg/merged_1kg_spiromics.vcf.gz'
scripts/pca_ld_pruning_merged_vcf.sh ${VCF}

# Run smartpca from eigenstrat (same options as in GTEx v8)
# Using option `-w poplistname` to infer eigenvectors using only individuals from a subset of populations, and then project individuals from all populations onto those eigenvectors.
FILE='genotype/pca/1kg/merged_1kg_spiromics.vcf.gz.filtered.ld_pruned'
smartpca.perl -i ${FILE}.ped -a ${FILE}.map -b ${FILE}.pedind -w genotype/pca/1kg/poplist.txt -k 20 -m 0 -o ${FILE}.pca -e ${FILE}.eval -p ${FILE}.plot -l ${FILE}.log

# Use KNN to inferre the population of the SPIROMICS samples
scripts/pca_predict_1kg_pop_knn.R
```

Make figures

```bash
scripts/pca_plot_projected_onto_1kg.R
```

### __Final number of PCs to use__

Number of PCs to be used: __4__:

* Top 4 PCs proportionally capture the most variance across SPIROMICS samples, PVE > 0.1.
* Top 4 PCs significantly correlate with subpopulations inferred for SPIROMICS samples, using 1000 Genomes Project samples and k-nearest neighbors clustering (P < 1 × 10−6 from F test; adjusted R2 = 0.36 − 0.98)

```bash
scripts/pca_select_pcs.R
```
