import os

import pandas as pd
import yaml

from helper import write_header
from helper import get_modification_color
from helper import parse_excel


def tsv2bedRMod(input_file, config_yaml, output_file):
    """
    converts tab-seperated files into bedMod format. only works with special columns as of now.
    These columns are: "chr", "pos", "strand", "motif", "frac"
    :param config_yaml: path/to/config.yaml, containing the header information for the new bedRMod file.
    :param input_file: path/to/input_file.tsv(.gz)
    :param output_file: path/to/output_file.bedrmod
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


def proEUF2bedRMod(input_file, config_yaml, output_file):
    """
    converts proEUF into bedRMod. needed second file to show at which positions are which modifications. This file is linked to in the config file
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
                    "\tfrequency\trefBase\n")
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
                    # custom = ""
                else:
                    selected_row = selected_row.iloc[0]
                    name = selected_row["mod_type"]
                    score = 954
                    frequency = 100
                    item_rgb = get_modification_color(name)
                    # custom = f"p-value={selected_row['p-value']};experimenter=John;"
                thick_start = row['thickStart']
                thick_end = row['thickEnd']
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
                # custom = ""
                f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                        f'\t{coverage}\t{frequency}\t{refBase}\n')


def bid_mouse2bedRMod(input_file, config_yaml, output_file=None, sheet_name=0):
    """
    converts bid-seq files into EUF.
    :param input_file: (path to) input file.
    :param config_yaml: (path to) config file containing the information on the metadata
    :param output_file: (path to) output file.
    :param sheet_name: if multiple sheets within an excel file are present, they can be passed by sheet name.
    :return:
    """
    bid = pd.read_excel(input_file, sheet_name=sheet_name, header=3)
    if output_file is None:
        output_file = input_file
    path, ending = os.path.splitext(output_file)

    if sheet_name != 0:
        sheet_name = sheet_name.replace(" ", "")
        path += f"_{sheet_name}"

    if not ending == ".bedrmod":
        output_file = path + ".bedrmod"
        print(f"filename changed to {output_file}")
    else:  # this is the case if a sheet name is added to the filename
        output_file = path + ending

    directory, file = os.path.split(output_file)
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"the given path does not lead to a directory: {directory}")

    config = yaml.safe_load(open(config_yaml, "r"))
    with open(output_file, 'w') as f:
        write_header(config, f)
        f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                "\tfrequency\trefBase\tcustom\n")
        for _, row in bid.iterrows():
            chrom = row["chr"]
            start = pd.to_numeric(row['pos'])
            end = start + 1
            name = "psi"
            score = row["Frac_Ave %"] * 100
            strand = row['strand']
            thick_start = start
            thick_end = end
            item_rgb = get_modification_color(name)
            # coverage is not directly indicated, but can be reconstructed from the given values
            coverage = ((row["Deletion_count_rep1"] / row["Deletion_rep1"])
                        + (row["Deletion_count_rep2"] / row["Deletion_rep2"]))
            frequency = row["Frac_Ave %"]
            refBase = "U" if row["Motif_1"][2].upper() == "T" else row["Motif_1"][2].upper()
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                    f'\t{coverage}\t{frequency}\t{refBase}\n')


def bid_human2bedRMod(input_file, config_yaml, output_file=None, sheet_name=0):
    """
    converts bid-seq files into EUF.
    :param input_file: (path to) input file.
    :param config_yaml: (path to) config file containing the information on the metadata
    :param output_file: (path to) output file.
    :param sheet_name: if multiple sheets within an excel file are present, they can be passed by sheet name.
    :return:
    """
    bid = pd.read_excel(input_file, sheet_name=sheet_name, header=3)
    if output_file is None:
        output_file = input_file
    path, ending = os.path.splitext(output_file)

    if sheet_name != 0:
        sheet_name = sheet_name.replace(" ", "")
        path += f"_{sheet_name}"

    if not ending == ".bedrmod":
        output_file = path + ".bedrmod"
        print(f"filename changed to {output_file}")
    else:  # this is the case if a sheet name is added to the filename
        output_file = path + ending

    directory, file = os.path.split(output_file)
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"the given path does not lead to a directory: {directory}")

    config = yaml.safe_load(open(config_yaml, "r"))
    with open(output_file, 'w') as f:
        write_header(config, f)
        f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                "\tfrequency\trefBase\tcustom\n")
        for _, row in bid.iterrows():
            chrom = row["chr"]
            start = pd.to_numeric(row['pos'])
            end = start + 1
            name = "psi"
            score = row["Frac_Ave %"] * 100
            strand = row['strand']
            thick_start = start
            thick_end = end
            item_rgb = get_modification_color(name)
            # coverage is not directly indicated, but can be reconstructed from the given values
            coverage = ((row["Deletion count_rep1"] / row["Deletion_rep1"])
                        + (row["Deletion count_rep2"] / row["Deletion_rep2"])
                        + (row["Deletion count_rep3"] / row["Deletion_rep3"]))
            frequency = row["Frac_Ave %"]
            refBase = "U" if row["Motif_1"][2].upper() == "T" else row["Motif_1"][2].upper()
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                    f'\t{coverage}\t{frequency}\t{refBase}\n')


def etam2bedRMod(input_file, config_yaml, output_file=None):
    """
    converts etam-seq files into EUF.
    :param input_file: (path to) input file.
    :param config_yaml: (path to) config file containing the information on the metadata
    :param output_file: (path to) output file.
    :return:
    """
    etam = pd.read_csv(input_file, delimiter="\t")

    if output_file is None:
        output_file = input_file

    path, ending = os.path.splitext(output_file)
    if not ending == ".bedrmod":
        output_file = path + ".bedrmod"
        print(f"output file: {output_file}")

    directory, file = os.path.split(output_file)
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"the given path does not lead to a directory: {directory}")

    config = yaml.safe_load(open(config_yaml, "r"))

    with open(output_file, 'w') as f:
        write_header(config, f)
        f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                "\tfrequency\trefBase\n")
        for _, row in etam.iterrows():
            ref_seg, pos, strand = row["pos"].split("_")
            chrom = ref_seg
            start = int(pos)
            end = start + 1
            name = "m6A"
            score = int(1000 - (row["FDR"] * 1000))
            strand = strand
            thick_start = start
            thick_end = end
            item_rgb = get_modification_color(name)
            # coverage is not directly indicated, but can be reconstructed from the given values
            coverage = (row["FTOm_total_count"] * row["accessbility"])
            frequency = row["methylation"]
            refBase = "U" if row["motif"][2].upper() == "T" else row["motif"][2].upper()
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                    f'\t{coverage}\t{frequency}\t{refBase}\n')


def csv2bedRMod(input_file, config_yaml, delimiter=None, ref_seg="ref_seg", start="pos", modi="m1A", modi_column=False,
                score=None, score_function=None, strand="strand", coverage=None, coverage_function=None, frequency=None,
                frequency_function=None, refBase=None, refbase_function=None):
    """
    converts arbitrary csv files into bedRMod format.
    The parameters usually pass the column name of the csv which contains the respective information.
    The name of the output file is infered from the input file and put in the same directory as the input file.
    :param input_file:(path to) input csv file.
    :param config_yaml: (path to) config file containing the information on the metadata
    :param delimiter: delimiter of the passed csv file. If "None" is it infered by pandas.
    :param ref_seg: column name of the column containing the reference sequence. i.e. the chromosome
    :param start: column name of the column that contains the positions of the modification
    :param modi: contains the column name of the column containing the modification or the name of the
    modification for the whole file.
    :param modi_column: indicates whether the value passed to "modi" contains the column name containing the
    modification (True) or denominates the modifiation itself (False)
    :param score: can either be a fixed value e.g. 0 or 1000 if score is totally unknown, or can be a calculation of the
    score e.g. int(1000 - (row["FDR"] * 1000)) or indicate a column name containing the score.
    :param score_function:
    :param strand: indicates the column name of the column containing the strandedness. can be "+" or "-" to indicate
    same strandedness for whole file.
    :param coverage: column name of column containing the coverage at this position.
    :param coverage_function:
    :param frequency:
    :param frequency_function:
    :param refBase:
    :param refbase_function:
    :return:
    """
    file = pd.read_csv(input_file, delimiter=delimiter)
    # file = pd.read_excel(input_file, header=3)

    output_file = input_file
    path, ending = os.path.splitext(output_file)
    if not ending == ".bedrmod":
        output_file = path + ".bedrmod"
        print(f"output file: {output_file}")

    config = yaml.safe_load(open(config_yaml, "r"))

    with open(output_file, 'w') as f:
        write_header(config, f)
        f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                "\tfrequency\trefBase\n")
        for _, row in file.iterrows():
            chrom = row[ref_seg]
            # chrom = ref_seg
            # start_col = int(pos)
            start_col = int(row[start])
            end = start_col + 1
            name = row[modi] if modi_column else modi
            # does this make sense?
            if score_function is not None:
                if type(score) == list:
                    params = [row[col] for col in score]
                elif isinstance(score, str):
                    params = row[score]
                else:
                    params = score
                score_column = score_function(params)
            else:
                if isinstance(score, str):
                    score_column = row[score]
                else:
                    score_column = score
            if strand == "+":
                strandedness = "+"
            elif strand == "-":
                strandedness = "-"
            else:
                strandedness = row[strand]
            thick_start = start_col
            thick_end = end
            item_rgb = get_modification_color(name)
            if coverage_function is not None:
                if isinstance(coverage, list):
                    params = [row[col] for col in coverage]
                elif isinstance(coverage, str):
                    params = row[coverage]
                coverage_col = coverage_function(params)
            else:
                if coverage in file.columns:
                    coverage_col = row[coverage_col]
                else:
                    coverage_col = coverage
            if frequency_function is not None:
                if isinstance(frequency, list):
                    params = [row[col] for col in frequency]
                elif isinstance(frequency, str):
                    params = row[frequency]
                frequency_col = frequency_function(params)
            else:
                if isinstance(frequency, str):
                    frequency_col = row[frequency]
                elif isinstance(frequency, (int, float)):
                    frequency_col = frequency
            if refbase_function is not None:
                if isinstance(refBase, list):
                    params = [row[col] for col in refBase]
                elif isinstance(refBase, str):
                    params = row[refBase]
                refBase_col = refbase_function(params)
            else:
                if refBase in file.columns:
                    refBase_col = row[refBase]
                else:
                    refBase_col = refBase
            f.write(f'{chrom}\t{start_col}\t{end}\t{name}\t{score_column}\t{strandedness}\t{thick_start}\t{thick_end}'
                    f'\t{item_rgb}\t{coverage_col}\t{frequency_col}\t{refBase_col}\n')


if __name__ == "__main__":
    proEUF2bedRMod("test_files/MH1601_both_GCF_ref_localN1L10nofwD20R3k1.proEUF", "config.yaml",
                   "example_files/test_frankenstein.bedrmod")
