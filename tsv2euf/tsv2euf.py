import os
import sys
import numpy as np
import pandas as pd
import yaml

EUF_VERSION = "bedModv1.2"


def write_header(file):
    config = yaml.load(open("config.yaml"), Loader=yaml.FullLoader)
    euf_header_keys = [
        "fileformat",
        "organism",
        "modification_type",
        "assembly",
        "annotation_source",
        "annotation_version",
        "sequencing_platform",
        "basecalling",
        "bioinformatics_workflow",
        "experiment",
    ]

    # metadata
    euf_header = dict()
    for key in euf_header_keys:
        euf_header[key] = config["options"].get(key, None)
    euf_header["fileformat"] = EUF_VERSION
    for k, v in euf_header.items():
        file.write(f"#{k}={v}\n")


def tsv2euf():
    tsv = pd.read_csv("../flat2euf/m6aSACseq/GSE198246/GSE198246_2ng_sites.tsv.gz", delimiter="\t")#,
                      #names=['chrom', 'pos', 'strand', 'motif', 'frac'])
    tsv['score'] = tsv['frac'] * 100

    # Set thickStart and thickEnd to pos and pos+1, respectively
    tsv['thickStart'] = tsv['pos']
    # convert dtype of columns
    # print(tsv.convert_dtypes().info())
    tsv["pos"] = pd.to_numeric(tsv["pos"], errors="coerce")
    # print(tsv["pos"])
    tsv['thickEnd'] = tsv["pos"] + 1

    # Drop the original frac column
    tsv = tsv.drop(columns=['frac'])

    # Write output file in BED format
    with open('../flat2euf/m6aSACseq/euf/output_file.bed', 'w') as f:
        write_header(f)
        f.write("chrom\tchromStart\tchromEnd\tname\tscore\tthickStart\tthickEnd\titemRgb\tcoverage\tfrequency"
                "\trefBase\n")
        for _, row in tsv.iterrows():
            #print(row)
            chrom = row['chr']
            start = row['pos']
            end = start + 1
            # name = row['motif']
            name = "."
            score = 0
            strand = row['strand']
            thick_start = row['thickStart']
            thick_end = row['thickEnd']
            item_rgb = '0,0,0'
            coverage = "-1"
            frequency = row['score']
            refBase = ""
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                    f'\t{coverage}\t{frequency:.5f}\t{refBase}\n')


if __name__ == "__main__":
    tsv2euf()
    # write_header()
