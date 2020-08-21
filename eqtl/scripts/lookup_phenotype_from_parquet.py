#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Lookup all rows containing given phenotype ID from a parquet file."""

import os
import sys
import argparse
import pandas as pd

def _set_args(): # type: (None) -> argparse.ArgumentParser
    parser = argparse.ArgumentParser() # type: argparse.ArgumentParser
    parser.add_argument(
        '-p',
        '--parquet',
        dest='parquet',
        required=True,
        type=str,
        metavar='parquet',
        help="Input file as parquet file"
    )
    parser.add_argument(
        #'-g',
        '--str1-list',
        dest='str1_list',
        required=True,
        type=str,
        metavar='str1 list',
        help="String 1 list"
    )
    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        type=str,
        metavar='output',
        help="Output file"
    )
    return parser


def read_phenotypes(phenotypes): # type: (str) -> Set[str]
    """Read a list of phenotypes"""
    pheno_list = list() # type: List
    with open(phenotypes) as gfile:
        for line in gfile: # type: str
            line = line.strip() # type: str
            pheno_list.append(line)
    return set(pheno_list)


def main(): # type: (None) -> None
    """Run the program"""
    parser = _set_args() # type: argparse.ArgumentParser
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args()) # type: Dict[str, Any]
    # Read in list of phenotypes to be selected
    str1_list = read_phenotypes(args['str1_list']) # type: Set[str]
    # Read in parquet file
    df = pd.read_parquet(args['parquet'])
    # Select rows where phenotype_id is present in the list
    df = df.loc[(df['phenotype_id'].isin(str1_list))]
    # Add chr and pos information (needed for coloc to make ID's)
    df['chr'] = df.variant_id.str.split("_").str[0]
    df['pos'] = df.variant_id.str.split("_").str[1]
    # Write file out
    df.to_csv(args['output'], sep='\t', compression='gzip', index=False)


if __name__ == '__main__':
    main()
