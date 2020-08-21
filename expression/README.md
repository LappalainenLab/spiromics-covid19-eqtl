# eQTL mapping in SPIROMICS bronchial samples for genes related to COVID-19 response: processing gene expression data

This repository documents the workflow and code used to process gene expression data and calculate PEER factors for eQTL mapping.

Using the following software/packages:

* R/3.4.1
* python/3.6.4
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) R package
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) R package
* [peer](https://github.com/PMBio/peer/wiki) R package

## c) Normalized gene expression data

Processing gene expression data according to the [GTEx v8 pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl) ([pre-print](https://doi.org/10.1101/787903) in bioRxiv):

### __Generate normalized expression data in BED format__

The expression data are normalized as follows:

* Read counts are normalized between samples using TMM ([Robinson & Oshlack, 2010](https://doi.org/10.1186/gb-2010-11-3-r25)) from the `edgeR` package
* Genes are selected based on the following expression thresholds:
    * \>= 0.1 TPM in >= 20% samples, and
    * \>=6 reads (unnormalized) in >=20% samples
* Each gene is inverse normal transformed across samples

```bash
NORMALIZATION_METHOD="tmm"
TPM_FILE="/path/to/spiromics_epibrush_tpm.txt"
COUNTS_FILE="/path/to/spiromics_epibrush_counts.txt"
LINKING_FILE="/path/to/linking_file.txt" # linking file: expression_id (first column), genotype_id (second column) (if needed)
ANNOT_FILE="/path/to/gencode_v33_annotation.gtf"
OUT_FILE="normalized_expression.bed"

scripts/process_expression_data.R ${NORMALIZATION_METHOD} ${TPM_FILE} ${COUNTS_FILE} ${LINKING_FILE} ${ANNOT_FILE} ${OUT_FILE}
```

`DESeq2` size factors ([Love et al. 2014](https://doi.org/10.1186/s13059-014-0550-8)) normalized file for aFC calculations:

```bash
NORMALIZATION_METHOD="deseq2"
OUT_FILE="deseq_log2_normalized_expression.bed"

scripts/process_expression_data.R ${NORMALIZATION_METHOD} ${TPM_FILE} ${COUNTS_FILE} ${LINKING_FILE} ${ANNOT_FILE} ${OUT_FILE}
```

Tabix index .bed files

```bash
bgzip ${BED_FILE}
tabix -p bed ${BED_FILE}.gz
```

## d) PEER factors

Probabilistic Estimation of Expression Residuals (PEER, [Stegle et al. 2010](https://doi.org/10.1371/journal.pcbi.1000770), [Stegle et al. 2012](https://doi.org/10.1038/nprot.2011.457)) is a Bayesian framework that uses factor analysis methods to infer hidden determinants and their effects from gene expression profiles.

Thus, PEER factors are latent factors that capture and correct for technical variation (e.g., batch effects) and unwanted biological variation (e.g., cellular heterogeneity) in molecular datasets.

### __Estimate PEER factors__

Using the `peer` R package

```bash
INPUT='/path/to/normalized_expression.bed'
COV='no_cov' # one could also input known covariates to take into account

scripts/peer_estimate_factors.R ${INPUT} ${COV}
```

### __Optimal number of PEERs__

According to the optimization analysis for selection of PEERs by sample size to maximize cis-eGene discovery done in [GTEx v8](https://doi.org/10.1101/787903), 15 PEERs were chosen to be used as covariates in eQTL mapping together with 4 genotype PCs and sex.

```bash
# 15 PEERs + genotype PCs and sex
K='15'
PEER_FILE='/path/to/peer_factors.txt'
GENOPC_FILE='/path/to/genotype_pcs_and_sex.txt'

scripts/prepare_cov.R ${K} ${PEER_FILE} ${GENOPC}
```

Alternatively, one could run the cis-eQTL discovery pipeline with increments of 5 PEERs and with reduced number of permutations (100 instead of 10,000) to assess the optimal number of PEERs.
