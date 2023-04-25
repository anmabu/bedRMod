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

    # build the header from metadata
    euf_header = dict()
    for key in euf_header_keys:
        euf_header[key] = config["options"].get(key, None)
    euf_header["fileformat"] = EUF_VERSION

    # check for additional keys and append them to the header
    additional_keys = []
    for key in config["options"].keys():
        if key not in euf_header_keys:
            additional_keys.append(key)
    # append additional keys
    if len(additional_keys) > 0:
        for key in additional_keys:
            # if there are nested dictionaries, they get appended here
            if isinstance(config["options"].get(key, None), dict):
                npairs = ""
                for nkey, nvalue in config["options"].get(key, None).items():
                    npairs += f"{nkey}:{nvalue};"
                npairs = npairs[:-1]
                euf_header[key] = npairs
            else:
                euf_header[key] = config["options"].get(key, None)
    for k, v in euf_header.items():
        file.write(f"#{k}={v}\n")


def tsv2euf():
    tsv = pd.read_csv("../flat2euf/m6aSACseq/GSE198246/GSE198246_2ng_sites.tsv.gz", delimiter="\t")
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
    # with open("../flat2euf/m6aSACseq/euf/output_header_file.bed", "w") as output:
    #     write_header(output)
