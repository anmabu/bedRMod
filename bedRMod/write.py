import pandas as pd
from ruamel.yaml import YAML

yaml = YAML()
yaml.sort_base_mapping_type_on_output = False  # disable sorting of keys


from helper import EUF_VERSION


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
    for k, v in euf_header.items():
        if isinstance(v, dict):
            npairs = ""
            for ke, va in v.items():
                npairs += f"{ke}:{va};"
            npairs = npairs[:-1]  # remove last ;
            output_file.write(f"#{k}={npairs}\n")
            continue  # don't write it twice
        if v is not None:
            output_file.write(f"#{k}={v}\n")
        else:
            value = ""
            # these are required fields in the config file
            if k in ["fileformat", "organism", "modification_type", "assembly", "annotation_source",
                     "annotation_version"]:
                print(f"There is a problem with the config.yaml file. {k} is required to have a value. Please correct "
                      f"this and convert again!")
                return False
            else:
                output_file.write(f"#{k}={value}\n")
    return True


def write_header_from_dict(file, header_dict):
    """
    Write header dict in correct format to new bedRMod file
    :param file: output bedrmod file
    :param header_dict: Dictionary containing the values of the header to write
    :return:
    """
    with open(file, 'w') as f:
        for key, value in header_dict.items():
            f.write('#' + key + '=' + value + '\n')
    return

def write_data_from_df(file, data_df):
    """
    append pandas data frame to bedRMod file, that already contains a header!
    :param file: output bedrmod file
    :param data_df: Dataframe with values to write to bedrmod file
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

