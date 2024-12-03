import pandas as pd

def read_header(file):
    """
    reads the header of a bedRMod file and returns a dictionary
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
    """
    bedrmod = pd.read_csv(file, sep="\t", comment="#", names=["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "coverage", "frequency"])
    return bedrmod

def read_bedRMod(file):
    """reads a bedRMod file an returns a (header, data) tuple"""
    return read_header(file), read_data(file)

if __name__ ==  "__main__":
    test_file = "/home/anne/Link to Dokumente/PhD/bedRMod/test/test_df2bedrmod.bedrmod"
    print(read_header("/home/anne/Link to Dokumente/PhD/bedRMod/test/test_df2bedrmod.bedrmod"))
    print(read_data(test_file))
    print(read_bedRMod(test_file))