# spiromics-covid19-eqtl

Mapping of cis-eQTLs for COVID-19-related genes in the SPIROMICS dataset

This repository contains code to conduct eQTL mapping in bronchial epithelium in the SPIROMICS cohort as presented in the article "_Genetic and non-genetic factors affecting the expression of COVID-19-relevant genes in the large airway epithelium_" ([Kasela et al., 2021, Genome Med](https://doi.org/10.1186/s13073-021-00866-2), see the preprint in [medRxiv](https://www.medrxiv.org/content/10.1101/2020.10.01.20202820v1)).

## Overview of this repository

* [genotype](https://github.com/LappalainenLab/spiromics-covid19-eqtl/tree/master/genotype) - scripts to process genotype data for eQTL mapping in SPIROMICS
* [expression](https://github.com/LappalainenLab/spiromics-covid19-eqtl/tree/master/expression) - scripts to process gene expression data and calculate PEER factors for eQTL mapping
* [candidate_genes](https://github.com/LappalainenLab/spiromics-covid19-eqtl/tree/master/candidate_genes) - selection of COVID-19 candidate genes
* [eqtl](https://github.com/LappalainenLab/spiromics-covid19-eqtl/tree/master/eqtl) - scripts to perform eQTL mapping and downstream analyses

Full eQTL summary statistics for the 496 COVID-19-related genes can be downloaded from [here](https://github.com/LappalainenLab/spiromics-covid19-eqtl/tree/master/eqtl/summary_stats).
