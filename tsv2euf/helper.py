import yaml

EUF_VERSION = "bedModv1.2"


def write_header(config, output_file):
    """
    reads information from the config yaml and writes it to the header of the bedMod file.
    the structure of the config file is quite rigid as of now.
    :param config: this .yaml file contains the options to write to the header.
    :param output_file: this is the file where the header is written into. File already has to be open for this to work! T
    :return:
    """
    # config = yaml.load(open(config_yaml), Loader=yaml.FullLoader)
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
