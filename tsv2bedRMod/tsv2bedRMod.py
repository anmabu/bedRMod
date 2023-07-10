import os

import pandas as pd
import yaml

from helper import write_header
from helper import get_modification_color


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

    config = yaml.safe_load(open(config_yaml, "r"))

    # Write output file in BED format
    with open(output_file, 'w') as f:
        write_header(config, f)
        f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage\tfrequency"
                "\trefBase\tcustom\n")
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
            custom = ""
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                    f'\t{coverage}\t{frequency}\t{refBase}\t{custom}\n')


def proEUF2euf(input_file, config_yaml, output_file):
    """
    converts profile into EUF. needed second file to show at which positions are which modifications
    :param input_file: (path to) input file in proEUF format.
    :param config_yaml: (path to) config file containing the information on the metadata
    :param output_file: (path to) output file.
    :return:
    """

    proEUF = pd.read_csv(input_file, delimiter="\t")

    # Set thickStart and thickEnd to pos and pos+1, respectively
    proEUF['thickStart'] = proEUF['pos']
    # convert dtype of columns
    proEUF["pos"] = pd.to_numeric(proEUF["pos"])
    proEUF['thickEnd'] = proEUF["pos"] + 1

    path, ending = os.path.splitext(output_file)
    if not ending == ".euf.bed":
        output_file = path + ".euf.bed"
        print(f"filename changed to {output_file}")

    directory, file = os.path.split(output_file)
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"the given path does not lead to a directory: {directory}")

    config = yaml.safe_load(open(config_yaml, "r"))
    if config["modifications_file"]:
        mod_file = pd.read_csv(config["modifications_file"])
        with open(output_file, 'w') as f:
            write_header(config, f)
            f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                    "\tfrequency\trefBase\tcustom\n")
            for _, row in proEUF.iterrows():
                chrom = row["ref_seg"]
                start = row['pos']
                end = start + 1
                score = 0
                strand = row['strand']
                # print(mod_file.loc[(mod_file["ref_seg"] == chrom) & (mod_file["mod_index"] == start)])
                selected_row = mod_file[(mod_file["ref_seg"] == chrom) & (mod_file["mod_index"] == start)]
                if selected_row.empty:
                    name = "."
                    frequency = 0
                    item_rgb = '0,0,0'
                    custom = ""
                else:
                    selected_row = selected_row.iloc[0]
                    name = selected_row["mod_type"]
                    score = 954
                    frequency = 100
                    item_rgb = get_modification_color(name)
                    custom = f"p-value={selected_row['p-value']};experimenter=John;"
                thick_start = row['thickStart']
                thick_end = row['thickEnd']
                coverage = row["cov"]
                refBase = "U" if row["ref_base"].upper() == "T" else row["ref_base"].upper()
                f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                        f'\t{coverage}\t{frequency}\t{refBase}\t{custom}\n')

    else:
        # Write output file in BED format without modifications
        with open(output_file, 'w') as f:
            write_header(config, f)
            f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                    "\tfrequency\trefBase\tcustom\n")
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
                custom = ""
                f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                        f'\t{coverage}\t{frequency}\t{refBase}\t{custom}\n')


if __name__ == "__main__":
    proEUF2euf("test_files/MH1601_both_GCF_ref_localN1L10nofwD20R3k1.proEUF", "config.yaml",
               "example_files/test_frankenstein.txt")
