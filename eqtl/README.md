# eQTL mapping in SPIROMICS bronchial samples for genes related to COVID-19 response: eQTL mapping and downstream analyses

This repository documents the workflow and code used to map eQTLs for COVID-19-related genes in SPIROMICS bronchial epithelium samples and perform downstream analyses.

Using the following software/packages:

* R/3.4.1
* python/3.6.4
* [tensorQTL](https://github.com/broadinstitute/tensorqtl)
* [aFC](https://github.com/secastel/aFC)
* [phenoscanner](https://github.com/phenoscanner/phenoscanner) R package
* [coloc](https://github.com/chr1swallace/coloc/tree/condmask) R package (`condmask` branch)

## e) Map eQTLs with tensorQTL

Mapping eQTLs using the tensorQTl software ([Taylor-Weiner et al. 2019](https://doi.org/10.1186/s13059-019-1836-7))

### Cis-eQTL mapping

-**_Mapping of cis-eQTLs with tensorQTL_**

* 15 PEERs + 4 genotype PCs, sex
* variants with minor allele frequency >= 0.05
* window-size +/- 1Mb of TSS
* 10,000 permutations to control for multiple testing
* gene-level q-values with a fixed P-value interval for the estimation of pi0 (the ‘lambda’ parameter was set to 0.85)
* False discovery rate (FDR) threshold of < 0.05 for significant eGenes

```bash
# Input data
GENO_PATH='/path/to/plink_genotype_file'
PHENO_FILE='/path/to/normalized_expression.bed.gz'
COV_FILE="/path/to/covariates.txt"
PREFIX="spiromics"
OUT="/path/to/output_folder"

# Phenotype-level summary statistics of cis-QTLs with empirical p-values
python3 -m tensorqtl ${GENO_PATH} ${PHENO_FILE} ${PREFIX} \
    --mode cis \
    --covariates ${COV_FILE} \
    --permutations 10000 \
    --window 1000000 \
    --fdr 0.05 \
    --qvalue_lambda 0.85 \
    --seed 124456677 \
    --output_dir ${OUT}

# Summary statistics of cis-QTLs for all gene-variants pairs
python3 -m tensorqtl ${GENO_PATH} ${PHENO_FILE} ${PREFIX} \
    --mode cis_nominal \
    --covariates ${COV_FILE} \
    --window 1000000 \
    --output_dir ${OUT}
```

-**_Mapping of conditionally independent cis-eQTLs with tensorQTL_**

Independent eQTLs at FDR 5% were identified by forward stepwise regression followed by a backwards selection step. The probe-level significance threshold was set to be the maximum beta-adjusted P-value (correcting for multiple-testing across the variants) over all eGenes.

* Forward stage: performing a scan for cis-eQTLs, correcting for all previously discovered variants and all covariates used in regular cis-eQTL mapping.
   * If the beta adjusted P-value for the lead variant was significant at the probe-level threshold, the lead variant was added to the list of discovered cis-eQTLs.
* Backwards stage: testing each variant separately, controlling for all other discovered variants.
   * If no variant was significant at the probe-level threshold the variant in question was dropped, otherwise the lead variant from this scan, which controls for all other signals (except one) found in the forward stage, was chosen as the variant that represents the signal best in the full model.

```bash
# Conditionally independent cis-eQTLs
CIS_FILE="/path/to/spiromics.cis_qtl.txt.gz"

python3 -m tensorqtl ${GENO_PATH} ${PHENO_FILE} ${PREFIX} \
    --mode cis_independent \
    --covariates ${COV_FILE} \
    --cis_output ${CIS_FILE} \
    --permutations 10000 \
    --window 1000000 \
    --seed 124456677 \
    --output_dir ${OUT}
```

### aFC estimates

Using allelic fold change (aFC) as the measure of eQTL effect size ([Mohammadi et al. 2017](https://doi.org/10.1101/gr.216747.116))

* gene count data that was normalized with DESeq size factors and log2 transformed
* including the same covariates that were used in the cis-eQTL mapping

Note: aFC estimates are capped at log2(100) = 6.64 (_it's suggested to remove eQTLs where the aFC has it the cap value_)

```bash
scripts/calculate_afc.sh
```

### Summary of eQTL mapping

Overview of cis-eQTLs for COVID-19-related genes in bronchial epithelium

```bash
summary_of_eqtl_mapping.Rmd
```

Regional association plots for high interest genes - _ACE2_, _TMPRSS2_, and the six genes from the 3p21.31 locus associated with COVID-19 ([Ellinghaus et al. 2020](https://doi.org/10.1056/NEJMoa2020283) and [COVID-19 Host Genetics Initiative](https://www.covid19hg.org/))

```bash
summary_eqtl_regional_association_plots.Rmd
```

## f) Downstream analyses

Using PhenoScanner v2 ([Staley et al. 2016](https://doi.org/10.1093/bioinformatics/btw373), [Kamat et al. 2019](https://doi.org/10.1093/bioinformatics/btz469)) to search if the regulatory variants for COVID-19-related genes have phenotype associations and coloc ([Giambartolomei et al. 2014](https://doi.org/10.1371/journal.pgen.1004383), [Wallace 2020](https://doi.org/10.1371/journal.pgen.1008720)) to run colocalization analysis

### PhenoScanner

Querying the PhenoScanner database of genotype-phenotype associations and overview of the results

```bash
summary_phenoscanner_lookup.Rmd
```

### Colocalization

Test for colocalization if the regulatory variant for COVID-19-related gene has also association with blood cell traits (_EFO parent category = Hematological measurement_), pulmonary function traits (_EFO parent category = Pulmonary function measurement_), or respiratory system disease (_EFO parent category = Respiratory disease_)

Using both the standard method and coloc-cond/mask (development version):

  * When using coloc-cond/mask, then using `method = mask` for the GWAS dataset and LD from the corresponding 1000G population to match the ancestry of the discovery population (_for example, 1000G CEU population if the discovery population is of European ancestry_)
  * When using coloc-cond/mask, then using `method = single` for the eQTL dataset if the eGene does not have multiple independent eVariants, otherwise using conditional p-values (_all the eGenes tested in coloc had only one independent eVariant_)
   * Run coloc in 500-kb region centered on the lead eQTL (+/-250 kb from the lead variant)
   * Priors p1 = 1e-4, p2 = 1e-4, p3 = 5e-6

-**_Processing data_**

Get full summary statistics for COVID-19-related genes that have pheWAS hits

```bash
# Using here a python script for subsetting parquet file
# Alternatively one could use the `arrow::read_parquet()` function to read parquet files into R
cat phewas_phenoscanner.covid19_related_egenes_v11.txt | cut -f 25 | grep -v "gene_id" | sort -u > covid19_related_egenes.txt

for CHR in {1..22} X; do
   scripts/lookup_phenotype_from_parquet.py --parquet /path/to/spiromics.cis_qtl_pairs.chr${CHR}.parquet --str1-list covid19_related_egenes.txt --output spiromics.cis_eqtl.covid19_related_egenes.allpairs.${CHR}.txt.gz
done
```

GWAS summary statistics from the GWAS catalog and by the Neale Lab using UKBB data

* Using the REST API to get GWAS summary statistics from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads/summary-statistics)
* UKBB GWAS data from the [Neale Lab](http://www.nealelab.is/uk-biobank), using the GWAS round 2 results
   * Note that UKBB GWAS data from the Neale Lab was in GRCh37, used `liftover` tool to lift the positions over to GRCh38

```bash
# https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291
# Download variants file
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

# Lift over positions from hg19 to GRCh38
scripts/ukbb_liftover_hg19_to_grch38.sh

# Add GRCh position to the variants file
scripts/ukbb_add_grch38pos.R

# Download phenotypes file to get # of cases
wget https://www.dropbox.com/s/d4mlq9ly93yhjyt/phenotypes.both_sexes.tsv.bgz -O phenotypes.both_sexes.tsv.gz
```

Calculate sdY for each of the eGenes

* sdY = standard deviation of residual expression (adjusted for the covariates used in eQTL mapping)

```bash
# Get sdY for expression
pheno_file="/path/to/normalized_expression.bed.gz"
cov_file="/path/to/covariates.txt"
cores="5"
out_file="spiromics_sdY.txt"
scripts/get_sdY.R ${pheno_file} ${cov_file} ${cores} ${out_file}
```

-**_coloc-standard_**

Assuming one causal variant per trait

```bash
# Example to run coloc-standard
METHOD='standard'
G="ERMP1"
scripts/run_coloc_for_phenoscanner_hits.R ${G} ${METHOD}
```

-**_coloc-cond/mask_**

Using masking in GWAS dataset, LD data based on the corresponding 1000G population

```bash
# Example to run coloc-cond/mask: masking in GWAS dataset
METHOD='mask'
G="ERMP1"
scripts/run_coloc_for_phenoscanner_hits.R ${G} ${METHOD}
```

-**_Summary of colocalization analysis_**

Overview of colocalization analysis with coloc

```bash
summary_coloc.Rmd
```
