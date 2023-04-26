import os
import sys
import numpy as np
import pandas as pd
import yaml

EUF_VERSION = "bedModv1.2"


def write_header(config_yaml, output_file):
    """
    reads information from the config yaml and writes it to the header of the bedMod file.
    the structure of the config file is quite rigid as of now.
    :param config_yaml: this .yaml file contains the options to write to the header.
    :param output_file: this is the file where the header is written into. File already has to be open for this to work! T
    :return:
    """
    config = yaml.load(open(config_yaml), Loader=yaml.FullLoader)
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
        output_file.write(f"#{k}={v}\n")


def tsv2euf(input_file, config_yaml, output_file):
    """
    converts tab-seperated files into bedMod format. only works with special columns as of now.
    These columns are: "chr", "pos", "strand", "motif", "frac"
    :param config_yaml: path/to/config.yaml, containing the header information for the new bedMod file.
    :param input_file: path/to/input_file.tsv(.gz)
    :param output_file: path/to/output_file.bed
    :return:
    """
    tsv = pd.read_csv(input_file, delimiter="\t")
    tsv['score'] = tsv['frac'] * 100

    # Set thickStart and thickEnd to pos and pos+1, respectively
    tsv['thickStart'] = tsv['pos']
    # convert dtype of columns
    tsv["pos"] = pd.to_numeric(tsv["pos"], errors="coerce")
    tsv['thickEnd'] = tsv["pos"] + 1

    # Drop the original frac column
    tsv = tsv.drop(columns=['frac'])

    # Write output file in BED format
    with open(output_file, 'w') as f:
        write_header(config_yaml, f)
        f.write("chrom\tchromStart\tchromEnd\tname\tscore\tthickStart\tthickEnd\titemRgb\tcoverage\tfrequency"
                "\trefBase\n")
        for _, row in tsv.iterrows():
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
    tsv2euf("../flat2euf/m6aSACseq/GSE198246/GSE198246_2ng_sites.tsv.gz", "config.yaml",
            "../flat2euf/m6aSACseq/euf/output_file.bed")
    # with open("../flat2euf/m6aSACseq/euf/output_header_file.bed", "w") as output:
    #     write_header(output)
