import os
import sys
import numpy as np
import pandas as pd


def tsv2euf():
    tsv = pd.read_csv("../flat2euf/m6aSACseq/GSE198246/GSE198246_2ng_sites.tsv.gz", delimiter="\t", header=None,
                      names=['chrom', 'pos', 'strand', 'motif', 'frac'])
    tsv['score'] = tsv['frac'] * 100

    # Set thickStart and thickEnd to pos and pos+1, respectively
    tsv['thickStart'] = tsv['pos']
    tsv['thickEnd'] = tsv['pos'] + 1

    # Drop the original frac column
    tsv = tsv.drop(columns=['frac'])

    # Write output file in BED format
    with open('output_file.bed', 'w') as f:
        for _, row in tsv.iterrows():
            chrom = row['chrom']
            start = row['pos']
            end = start + 1
            name = row['motif']
            score = row['score']
            strand = row['strand']
            thick_start = row['thickStart']
            thick_end = row['thickEnd']
            item_rgb = '0,0,0'
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score:.2f}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}\n')

