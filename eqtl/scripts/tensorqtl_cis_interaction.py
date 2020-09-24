#!/usr/bin/env python3

"""Run interaction-QTL mapping with EigenMT correction"""

# Script modified from https://github.com/broadinstitute/tensorqtl

import os
import sys
import argparse

import numpy as np
import pandas as pd

import tensorqtl
from tensorqtl import genotypeio, cis, eigenmt

# def main():
parser = argparse.ArgumentParser()
parser.add_argument(
    '-i',
    '--interaction',
    dest='interaction',
    type=str,
    required=True,
    default=None,
    metavar='interaction',
    help="Provide interaction term"
)
if not sys.argv[1:]:
    sys.exit(parser.print_help())

args = vars(parser.parse_args()) # type: Dict[str, Any]

# Input data
geno_path = 'genotype/freeze9.pass_only.spiromics_144samples.maf05.biallelic'
pheno_file = 'expression/normalized/spiromics.normalized_expression.bed.gz'
cov_file = 'cov/spiromics.covariates.txt'
interaction_file = 'cov/cibersort_interaction/spiromics.%(interaction)s_cibersort.txt' % args
prefix = 'spiromics'
output_dir = 'eqtl/cis_interaction/cibersort/results/%(interaction)s' % args

# Load phenotypes and covariates:
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(pheno_file)
covariates_df = pd.read_csv(cov_file, sep='\t', index_col=0).T # samples x covariates
assert np.all(phenotype_df.columns==covariates_df.index)

# Load interaction data
interaction_s = pd.read_csv(interaction_file, sep='\t', index_col=0, header=None, squeeze=True)
## Select individuals that are in the interaction dataset
phenotype_df = phenotype_df.iloc[:, phenotype_df.columns.isin(interaction_s.index)]
covariates_df = covariates_df[covariates_df.index.isin(interaction_s.index)]
assert np.all(phenotype_df.columns==covariates_df.index)
assert covariates_df.index.isin(interaction_s.index).all()
interaction_s = interaction_s.loc[covariates_df.index].astype(np.float32)

# Load genotypes (for VCFs with hard GT calls only, specify type as np.int8 to save memory)
pr = genotypeio.PlinkReader(geno_path, select_samples=phenotype_df.columns, dtype=np.int8)

# Load genotypes for each chromosome separately
top_df = []
for chrom in pr.chrs:
    g, pos_s = pr.get_region(chrom)
    genotype_df = pd.DataFrame(g, index=pos_s.index, columns=pr.fam['iid'])[phenotype_df.columns]
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    # Map cis_nominal with intercation term and eigenMT correction
    chr_df = cis.map_nominal(genotype_df, variant_df[variant_df['chrom']==chrom], phenotype_df[phenotype_pos_df['chr']==chrom], phenotype_pos_df[phenotype_pos_df['chr']==chrom], covariates_df, prefix, interaction_s=interaction_s, maf_threshold_interaction=0.1, window=1000000, output_dir=output_dir, write_top=False, run_eigenmt=True)
    top_df.append(chr_df)

top_df = pd.concat(top_df)
top_df.to_csv(os.path.join(output_dir, '{}.cis_qtl_top_assoc.txt.gz'.format(prefix)), sep='\t', float_format='%.6g')

# if __name__ == '__main__':
#     main()
