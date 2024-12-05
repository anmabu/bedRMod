import pandas as pd


def read_header(file):
    """
    reads the header of a bedRMod file and returns a dictionary
    :param file: (path to) bedRMod file
    :return: header of bedRMod file
    """
    header_dict = {}
    with open(file, "r") as f:
        for line in f:
            if line.startswith("#chrom"):
                break  # stop at last line of header because these are the column names of the data section
            if line.startswith("#"):
                line = line[1:].rstrip()
                k, v = line.split("=")
                header_dict[k] = v
    return header_dict


def read_data(file):
    """
    reads the data section of a bedRMod file and returns a pd.DataFrame with its contents
    :param file: (path to) bedRMod file
    :return: pandas dataframe of data in bedRMod file
    """
    bedrmod = pd.read_csv(file, sep="\t", comment="#", names=["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "coverage", "frequency"])
    return bedrmod


def read_bedRMod(file):
    """
    reads a bedRMod file and returns a (header, data) tuple
    :param file: (path to) bedRMod file
    :return: tuple of header and pandas dataframe of bedRMod file
    """
    return read_header(file), read_data(file)
