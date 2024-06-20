import os
import pandas as pd
import yaml

# from tsv2bedRMod import bid_mouse2bedRMod, bid_human2bedRMod, etam2bedRMod, csv2bedRMod
from tsv2bedRMod import csv2bedRMod

from helper import write_header
from helper import get_modification_color
from helper import parse_excel


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
                "\tfrequency\n")
        for _, row in bid.iterrows():
            chrom = row["chr"]
            start = pd.to_numeric(row['pos'])
            end = start + 1
            name = "Y"
            score = round(row["Frac_Ave %"] * 10)
            strand = row['strand']
            thick_start = start
            thick_end = end
            item_rgb = get_modification_color(name)
            # coverage is not directly indicated, but can be reconstructed from the given values
            coverage = round(((row["Deletion_count_rep1"] / row["Deletion_rep1"])
                        + (row["Deletion_count_rep2"] / row["Deletion_rep2"])))
            frequency = round(row["Frac_Ave %"])
            # refBase = "U" if row["Motif_1"][2].upper() == "T" else row["Motif_1"][2].upper()
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                    f'\t{coverage}\t{frequency}\n')


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
                "\tfrequency\n")
        for _, row in bid.iterrows():
            chrom = row["chr"]
            start = pd.to_numeric(row['pos'])
            end = start + 1
            name = "Y"
            score = round(row["Frac_Ave %"] * 10)
            strand = row['strand']
            thick_start = start
            thick_end = end
            item_rgb = get_modification_color(name)
            # coverage is not directly indicated, but can be reconstructed from the given values
            coverage = round(((row["Deletion count_rep1"] / row["Deletion_rep1"])
                        + (row["Deletion count_rep2"] / row["Deletion_rep2"])
                        + (row["Deletion count_rep3"] / row["Deletion_rep3"])))
            frequency = round(row["Frac_Ave %"])
            # refBase = "U" if row["Motif_1"][2].upper() == "T" else row["Motif_1"][2].upper()
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                    f'\t{coverage}\t{frequency}\n')


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
                "\tfrequency\n")
        for _, row in etam.iterrows():
            ref_seg, pos, strand = row["pos"].split("_")
            chrom = ref_seg
            start = int(pos)
            end = start + 1
            name = "m6A"
            score = round(1000 - (row["FDR"] * 1000))
            strand = strand
            thick_start = start
            thick_end = end
            item_rgb = get_modification_color(name)
            # coverage is not directly indicated, but can be reconstructed from the given values
            coverage = round((row["FTOm_total_count"] * row["accessbility"] / 100))
            frequency = round(row["methylation"])
            #refBase = "U" if row["motif"][2].upper() == "T" else row["motif"][2].upper()
            f.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}'
                    f'\t{coverage}\t{frequency}\n')


def convert_bid_mouse():
    dirpath = "/home/annebusch/anne02/euf-data/bid-seq/"
    for file in os.listdir(dirpath):
        if "Mouse" in file and file.endswith(".xlsx"):
            sheet_names = parse_excel(dirpath + file)
            for name in sheet_names:
                bid_mouse2bedRMod(dirpath + file, "/home/annebusch/anne02/euf-data/bid-seq/mouse_config.yaml", sheet_name=name)


def convert_bid_human():
    dirpath = "/home/annebusch/anne02/euf-data/bid-seq/"
    for file in os.listdir(dirpath):
        if "Mouse" not in file and "shControl" not in file and file.endswith(".xlsx"):
            bid_human2bedRMod(dirpath + file, dirpath + "human_config.yaml")


def convert_etam():
    def score_func(param):
        return int(1000 - (param * 1000))

    def cov_function(params):
        one, two = params
        return one * two

    def refbase_func(param):
        return "U" if param[2].upper() == "T" else param[2].upper()
    dirpath = "/home/annebusch/anne02/euf-data/etam/"
    for file in os.listdir(dirpath):
        if "mesc" in file and file.endswith(".txt"):
            etam2bedRMod(dirpath + file, dirpath + "mouse_config.yaml")
        elif "hela" in file and file.endswith(".txt"):
            etam2bedRMod(dirpath + file, dirpath + "hela_config.yaml")


def rename_glori():
    """
    rename downloaded glori files, as the orignal filenames do not contain mouse/hela information
    :return:
    """
    dirpath = "/home/annebusch/anne02/euf-data/glori/"
    os.chdir(dirpath)
    rename_dict = {}
    with open(dirpath + "rename.txt", 'r') as r:
        for line in r:
            id, name = line.strip().split()
            rename_dict.__setitem__(id, name)
    for file in os.listdir(dirpath):
        if file.endswith(".csv"):
            id, rest = file.split()
            new_filename = id + rename_dict[id] + ".csv"
            os.rename(file, new_filename)


def convert_glori():
    def score_func(p_value):
        return round(1000 - (p_value * 1000))

    def frequency_func(ratio):
        return round(ratio * 100)

    dirpath = "/home/annebusch/anne02/euf-data/glori/"
    os.chdir(dirpath)
    for file in os.listdir():
        if file.endswith(".csv"):
            if "mESC" in file or "MEF" in file:
                conf = "mouse_config.yaml"
            elif "HeLa_si" in file or "STM" in file:
                conf = "hela-si_hela-stm_hek-stm_config.yaml"
            else:
                conf = "hek_hela-hypoxia_config.yaml"
            csv2bedRMod(file, conf,
                        delimiter="\t",
                        ref_seg="Chr", start="Sites", strand="Strand",
                        modi="m6A",
                        coverage="Acov",
                        score="Pvalue",
                        score_function=score_func,
                        frequency="Ratio",
                        frequency_function=frequency_func)

# convert_glori()
# convert_bid_human()
# convert_bid_mouse()
# convert_etam_human()
# convert_etam()
def score_func(param):
    return int(1000 - (param * 1000))

def cov_function(params):
    one, two = params
    return one * two
def score_func(p_value):
    return round(1000 - (p_value * 1000))

def frequency_func(ratio):
    return round(ratio * 100)

csv2bedRMod("/home/annebusch/anne02/euf-data/etam/GSE211303_hela.polya.wt.ftom.ftop.rep1.deep.hits.txt",
            "/home/annebusch/anne02/euf-data/etam/hela_config.yaml",
            delimiter="\t",
            modi="m6A",
            coverage=["FTOm_total_count", "accessbility"],
            coverage_function=cov_function,
            score="FDR",
            score_function=score_func,
            frequency="methylation",
            frequency_function=frequency_func)