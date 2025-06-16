import pandas as pd
from ruamel.yaml import YAML

yaml = YAML()
yaml.sort_base_mapping_type_on_output = False  # disable sorting of keys


from helper import EUF_VERSION, get_modification_color


def write_header_from_config(config_yaml, output_file):
    """
    reads information from the config yaml and writes it to the header of the bedMod file.
    the structure of the config file is quite rigid as of now.
    :param config_yaml: the contents of the config.yaml file that contains the options to write to the header.
    :param output_file: this is the file where the header is written into. File already has to be open for this to work!
    """

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
        "external_source"
    ]

    config = yaml.load(open(config_yaml, "r"))

    # build the header from metadata
    euf_header = dict()
    for key in euf_header_keys:
        euf_header[key] = config["options"].get(key, "")
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
            if isinstance(config["options"].get(key, ""), dict):
                npairs = ""
                for nkey, nvalue in config["options"].get(key, "").items():
                    npairs += f"{nkey}:{nvalue};"
                npairs = npairs[:-1]  # remove last ;
                euf_header[key] = npairs
            else:
                euf_header[key] = config["options"].get(key, "")

    with open(output_file, "w") as f:
        for k, v in euf_header.items():
            if isinstance(v, dict):
                npairs = ""
                for ke, va in v.items():
                    npairs += f"{ke}:{va};"
                npairs = npairs[:-1]  # remove last ;
                f.write(f"#{k}={npairs}\n")
                continue  # don't write it twice
            if v is not None:
                f.write(f"#{k}={v}\n")
            else:
                value = ""
                # these are required fields in the config file
                if k in ["fileformat", "organism", "modification_type", "assembly", "annotation_source",
                         "annotation_version"]:
                    print(f"There is a problem with the config.yaml file. {k} is required to have a value. Please correct "
                          f"this and convert again!")
                    return False
                else:
                    f.write(f"#{k}={value}\n")
    return True


def write_header_from_dict(header_dict, file):
    """
    Write header dict in correct format to new bedRMod file
    :param header_dict: Dictionary containing the values of the header to write
    :param file: output bedrmod file

    :return:
    """
    with open(file, 'w') as f:
        for key, value in header_dict.items():
            f.write('#' + key + '=' + value + '\n')
    return

def write_data_from_df(data_df, file):
    """
    append pandas data frame to bedRMod file, that already contains a header!
    :param data_df: Dataframe with values to write to bedrmod file
    :param file: output bedrmod file
    :return:
    """
    column_headers = True
    with open(file, 'r') as f:
        last = f.readlines()[-1]
        if not last.startswith("#chrom"):
            column_headers = False

    with open(file, 'a') as f:
        if not column_headers:
            f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                    "\tfrequency\n")
        data_df.to_csv(f, sep='\t', index=False, header=False, mode='a')

    return


def write_single_data_row(file, chrom, start, name, score, strand, coverage, frequency,
                          end=None, thickstart=None, thickend=None, itemRgb=None):
    """
    Appends contents of a single row to an already existing file.
    The positional arguments of this function are required for output. The keyword arugments can be inferred using the
    data from the positional arguments.
    :param file: output bedrmod file
    :param chrom: chromosome/reference segment information
    :param start: start position (0-based)
    :param name: modification name
    :param score: score (0 - 1000)
    :param strand: strand (+, -, .)
    :param coverage: coverage
    :param frequency: frequency
    :param end: end position, If not passed to the function, it will be calculated from the start position.
    :param thickstart: display start position. If not passed to the function, it will be calculated from the start position.
    :param thickend: display end position. If not passed to the function, it will be calculated from the end position.
    :param itemRgb: display color of the modification. If not passed, it will be calculated.
    :return:
    """

    if end is not None and end != start + 1:
        raise ValueError("end must be either None or start position + 1")
    if end is None:
        end = start + 1

    if thickstart is not None and thickstart != start:
        raise ValueError("thickstart must be either None or start position")
    if thickstart is None:
        thickstart = start

    if thickend is not None and ((thickend != end) or (thickend != thickstart + 1)):
        raise ValueError("thickend must be either None or (thick)start position + 1")
    if thickend is None:
        thickend = end

    # don't throw an error if the itemRGB is not the color from the list. This can be customized.
    if itemRgb is None:
        itemRgb = get_modification_color(name)


    # Do checks whether last header line is written correctly
    column_headers = True
    with open(file, 'r') as f:
        last = f.readlines()[-1]
        if not last.startswith("#chrom"):
            column_headers = False

    with open(file, 'a') as f:
        if not column_headers:
            f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                    "\tfrequency\n")
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thickstart}\t{thickend}\t{itemRgb}\t{coverage}\t{frequency}\n")

    return


def write_bedRMod(file, header_dict, data_df):
    """
    write header dict and pandas dataframe to new bedRMod file
    :param file: output bedrmod file
    :param header_dict: Dictionary containing the values of the header to write
    :param data_df: Dataframe with values to write to bedrmod file
    :return:
    """
    write_header_from_dict(file, header_dict)
    write_data_from_df(file, data_df)
    return

