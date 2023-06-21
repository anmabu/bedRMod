import os

import pandas as pd
import yaml

from helper import write_header


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
        config = yaml.load(open(config_yaml), Loader=yaml.FullLoader)
        write_header(config, f)
        f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage\tfrequency"
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
                    f'\t{coverage}\t{frequency}\t{refBase}\n')


def proEUF2euf(input_file, config_yaml, output_file):
    """
    converts profile into EUF. needed second file to show at which positions are which modifications
    :param input_file:
    :param config_yaml:
    :param output_file:
    :return:
    """

    proEUF = pd.read_csv(input_file, delimiter="\t")

    # Set thickStart and thickEnd to pos and pos+1, respectively
    proEUF['thickStart'] = proEUF['pos']
    # convert dtype of columns
    proEUF["pos"] = pd.to_numeric(proEUF["pos"])
    proEUF['thickEnd'] = proEUF["pos"] + 1

    config = yaml.load(open(config_yaml), Loader=yaml.FullLoader)
    if config["modifications_file"]:
        mod_file = pd.read_csv(config["modifications_file"])
        with open(output_file, 'w') as f:
            write_header(config, f)
            f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage\tfrequency"
                    "\trefBase\n")
            for _, row in proEUF.iterrows():
                chrom = row["ref_seg"]
                start = row['pos']
                end = start + 1
                score = 0
                strand = row['strand']
                selected_row = mod_file[(mod_file["ref_seg"] == chrom) & (mod_file["mod_index"] == start)]
                if selected_row.empty:
                    name = "."
                    frequency = 0
                else:
                    name = selected_row.iloc[0]["mod_type"]
                    frequency = 100
                thick_start = row['thickStart']
                thick_end = row['thickEnd']
                item_rgb = '0,0,0'
                coverage = row["cov"]
                refBase = "U" if row["ref_base"].upper() == "T" else row["ref_base"].upper()
                f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                        f'\t{coverage}\t{frequency}\t{refBase}\n')

    else:
        # Write output file in BED format without modifications
        with open(output_file, 'w') as f:
            write_header(config, f)
            f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                    "\tfrequency\trefBase\n")
            for _, row in proEUF.iterrows():
                chrom = row["ref_seg"]
                start = row['pos']
                end = start + 1
                name = "."
                score = 0
                strand = row['strand']
                thick_start = row['thickStart']
                thick_end = row['thickEnd']
                item_rgb = '0,0,0'
                coverage = row["cov"]
                frequency = 0
                refBase = "U" if row["ref_base"].upper() == "T" else row["ref_base"].upper()
                f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                        f'\t{coverage}\t{frequency}\t{refBase}\n')


if __name__ == "__main__":
    print(os.getcwd())
    proEUF2euf("test_files/MH1601_both_GCF_ref_localN1L10nofwD20R3k1.proEUF", "config.yaml",
               "example_files/test_frankenstein.bed")
