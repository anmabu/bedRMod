import yaml

EUF_VERSION = "bedRModv1.4"


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


def get_modification_color(modi):
    """
    looks up the color of the modification in the rgb dictionary and returns the assiciated rgb value
    :param modi:
    :return:
    """
    rgb_colors = {'?A': "255,0,0",
                  'm1A': "0,255,0",
                  'm2A': "0,0,255",
                  'i6A': "255,255,0",
                  'ms2i6A': "255,0,255",
                  'm6A': "0,255,255",
                  't6A': "128,0,0",
                  'm6t6A': "0,128,0",
                  'ms2t6A': "0,0,128",
                  'Am': "128,128,0",
                  'I': "128,0,128",
                  'm1I': "0,128,128",
                  'Ar(p)': "255,128,0",
                  'io6A': "255,0,128",
                  '?C': "128,255,0",
                  's2C': "0,255,128",
                  'Cm': "128,0,255",
                  'ac4C': "0,128,255",
                  'm5C': "255,128,128",
                  'm3C': "128,255,128",
                  'k2C': "128,128,255",
                  'f5C': "192,192,192",
                  'f5Cm': "192,0,0",
                  '?G': "192,192,0",
                  'm1G': "0,192,0",
                  'm2G': "192,0,192",
                  'Gm': "0,192,192",
                  'm22G': "0,0,192",
                  'm22Gm': "0,0,64",
                  'm7G': "0,64,0",
                  'fa7d7G': "64,0,0",
                  'Q': "0,205,139",
                  'manQ': "64,64,64",
                  'galQ': "64,64,0",
                  'yW': "64,0,64",
                  'o2yW': "0,64,64",
                  '?U': "255,64,0",
                  'mnm5U': "255,0,64",
                  's2U': "64,255,0",
                  'Um': "0,255,64",
                  's4U': "64,0,255",
                  'ncm5U': "0,64,255",
                  'mcm5U': "255,64,64",
                  'mnm5s2U': "64,255,64",
                  'mcm5s2U': "64,64,255",
                  'cmo5U': "205,205,205",
                  'mo5U': "139,0,0",
                  'cmnm5U': "139,139,0",
                  'cmnm5s2U': "0,139,0",
                  'acp3U': "139,0,139",
                  'mchm5U': "0,139,139",
                  'cmnm5Um': "0,0,139",
                  'ncm5Um': "0,0,70",
                  'D': "0,70,0",
                  'psi': "205,139,0",
                  'm1psi': "139,205,0",
                  'm5U': "205,0,139",
                  'm5s2U': "0,205,0",
                  'm5Um': "0,139,205"}
    return rgb_colors.get(modi)
