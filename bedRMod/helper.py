import math
import pandas as pd

import ruamel
from ruamel.yaml import YAML

yaml = YAML()  # don't enable safe mode as this erases the comments!
yaml.default_flow_style = False
yaml.sort_base_mapping_type_on_output = False  # disable sorting of keys


EUF_VERSION = "bedRModv1.8"


def write_config_header(config, output_file):
    """
    reads information from the config yaml and writes it to the header of the bedMod file.
    the structure of the config file is quite rigid as of now.
    :param config: the contents of the config.yaml file that contains the options to write to the header.
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


def funcify(expression):
    """
    Takes a string of an expression as an input and converts it into a python function.
    :return: function of passed expression string
    """
    eval_dict = {
        "log10": math.log10
    }
    func = eval(expression, eval_dict)
    return func


def write_bioinformatics_keys(config_yaml, workflow=None, coverage_function=None, frequency_function=None, score_function=None):
    """
    check if customizable functions are in the config and adds them to the config file if not included.
    :param config_yaml: (path to) config file
    :param workflow: string of the workflow that has been used
    :param score_function: string of score function
    :param coverage_function: string of coverage function
    :param frequency_function: string of frequency function
    :return:
    """

    # change representation of None in the output file to "", so that nothing gets written in bedRMod
    #class EmptyStringDumper(yaml.SafeDumper):
    #    def represent_none(self, _):
    #        return self.represent_scalar('tag:yaml.org,2002:str', '')
    def represent_none(self, _):
        return self.represent_scalar('tag:yaml.org,2002:str', '')
    # Add the custom representer to Yaml
    yaml.Representer.add_representer(type(None), represent_none)

    # EmptyStringDumper.add_representer(type(None), EmptyStringDumper.represent_none)

    config = yaml.load(open(config_yaml, "r"))
    with open(config_yaml + ".backup", 'w') as file:
        yaml.dump(config, file)
    try:
        if workflow is not None or score_function is not None or coverage_function is not None or frequency_function is not None:
            if type(config["options"]["bioinformatics_workflow"]) is not ruamel.yaml.comments.CommentedMap:
                config["options"]["bioinformatics_workflow"] = {"workflow": config["options"]["bioinformatics_workflow"]}
            if workflow is not None and "bioinformatics_workflow" not in config["options"].keys():
                config["options"]["bioinformatics_workflow"] = {"workflow": workflow}
            elif workflow is not None and "bioinformatics_workflow" in config["options"].keys():
                if "score_function" not in config["options"]["bioinformatics_workflow"].keys():
                    config["options"]["bioinformatics_workflow"]["workflow"] = workflow
                else:
                    if workflow != config["options"]["bioinformatics_workflow"]["workflow"]:
                        print(f"The workflow from the config file "
                              f"({config['options']['bioinformatics_workflow']['workflow']}) "
                              f"does not match the newly given workflow ({workflow}). "
                              f"Proceeding with the given workflow {workflow} "
                              f"and overwriting the config file.")
                        config["options"]["bioinformatics_workflow"]["workflow"] = workflow
            if coverage_function is not None and "bioinformatics_workflow" not in config["options"].keys():
                config["options"]["bioinformatics_workflow"] = {"coverage_function": coverage_function}
            elif coverage_function is not None and "bioinformatics_workflow" in config["options"].keys():
                if "coverage_function" not in config["options"]["bioinformatics_workflow"].keys():
                    config["options"]["bioinformatics_workflow"]["coverage_function"] = coverage_function
                else:
                    if coverage_function != config["options"]["bioinformatics_workflow"]["coverage_function"]:
                        print(f"The coverage function from the config file "
                              f"({config['options']['bioinformatics_workflow']['coverage_function']}) "
                              f"does not match the newly given coverage function ({coverage_function}). "
                              f"Proceeding with the given coverage function {coverage_function} "
                              f"and overwriting the config file.")
                        config["options"]["bioinformatics_workflow"]["coverage_function"] = coverage_function
            if frequency_function is not None and "bioinformatics_workflow" not in config["options"].keys():
                config["options"]["bioinformatics_workflow"] = {"frequency_function": frequency_function}
            elif frequency_function is not None and "bioinformatics_workflow" in config["options"].keys():
                if "frequency_function" not in config["options"]["bioinformatics_workflow"].keys():
                    config["options"]["bioinformatics_workflow"]["frequency_function"] = frequency_function
                else:
                    if frequency_function != config["options"]["bioinformatics_workflow"]["frequency_function"]:
                        print(f"The frequency function from the config file "
                              f"({config['options']['bioinformatics_workflow']['frequency_function']}) "
                              f"does not match the newly given frequency function ({frequency_function}). "
                              f"Proceeding with the given frequency function {frequency_function} "
                              f"and overwriting the config file.")
                        config["options"]["bioinformatics_workflow"]["frequency_function"] = frequency_function
            if score_function is not None and "bioinformatics_workflow" not in config["options"].keys():
                config["options"]["bioinformatics_workflow"] = {"score_function": score_function}
            elif score_function is not None and "bioinformatics_workflow" in config["options"].keys():
                if "score_function" not in config["options"]["bioinformatics_workflow"].keys():
                    config["options"]["bioinformatics_workflow"]["score_function"] = score_function
                else:
                    if score_function != config["options"]["bioinformatics_workflow"]["score_function"]:
                        print(f"The score function from the config file "
                              f"({config['options']['bioinformatics_workflow']['score_function']}) "
                              f"does not match the newly given score function ({score_function}). "
                              f"Proceeding with the given score function {score_function} "
                              f"and overwriting the config file.")
                        config["options"]["bioinformatics_workflow"]["score_function"] = score_function
        with open(config_yaml, 'w') as file:
            yaml.dump(config, file)
    except Exception as e:
        print("An exception occurred while trying to write the config file")
        with open(config_yaml + ".backup", 'r') as file:
            original_config = yaml.load(file)
        with open(config_yaml, 'w') as file:
            yaml.dump(original_config, file)


def read_bioinformatics_keys(config_yaml):
    """
    If the field "bioinformatics_workflow" is in the config file and contains an additional dictionary,
    the fields in this dictionary are read and returned.
    :param config_yaml: (path to) config file
    :return: values for workflow, score_function, coverage_function, frequency_function. Will return None if empty
    """
    workflow = None
    score_function = None
    coverage_function = None
    frequency_function = None
    config = yaml.load(open(config_yaml, "r"))
    if "bioinformatics_workflow" in config["options"].keys():
        if isinstance(config["options"].get("bioinformatics_workflow", ""), dict):
            for key in config["options"]["bioinformatics_workflow"].keys():
                if key == "workflow":
                    workflow = config["options"]["bioinformatics_workflow"]["workflow"]
                if key == "coverage_function":
                    coverage_function = config["options"]["bioinformatics_workflow"]["coverage_function"]
                if key == "frequency_function":
                    frequency_function = config["options"]["bioinformatics_workflow"]["frequency_function"]
                if key == "score_function":
                    score_function = config["options"]["bioinformatics_workflow"]["score_function"]

    return workflow, coverage_function, frequency_function, score_function


def check_value_range(result):
    """
    check whether returned values are in the allowed range
    :param result:
    :return:
    """
    chrom, start_col, end, name, score_column, strandedness, thick_start, thick_end, item_rgb, \
        coverage_col, frequency_col = result

    if not 0 <= score_column <= 1000:
        print(f"The score value ({score_column}) is not in the allowed range. Please check and try again.")
        return False

    if not 1 <= frequency_col <= 100:
        print(f"The frequency value ({frequency_col}) is not in the allowed range. Please check and try again.")
        return False
    return True

def get_modification_color(modi):
    """
    looks up the color of the modification in the rgb dictionary and returns the associated rgb value
    :param modi: short name of modification in Modomics
    :return: RGB value for modification
    """
    rgb_colors = {'pmnm5U': '255,0,0',
                     'm1Am': '0,255,0',
                     'pm1Am': '0,0,255',
                     'm1Gm': '255,255,0',
                     'pm1Gm': '255,0,255',
                     'm1Im': '0,255,255',
                     'pm1Im': '128,0,0',
                     'pse2U': '0,128,0',
                     'ps2U': '0,0,128',
                     'm1acp3Y': '128,128,0',
                     'pm1acp3Y': '128,0,128',
                     'm1A': '0,128,128',
                     'm1G': '255,128,0',
                     'm1I': '255,0,128',
                     'pm1I': '128,255,0',
                     'm1Y': '0,255,128',
                     'pm1Y': '128,0,255',
                     'pm1G': '0,128,255',
                     'pm5Um': '255,128,128',
                     'pAr(p)': '128,255,128',
                     'hm5Cm': '128,128,255',
                     'phm5Cm': '192,192,192',
                     'Am': '192,0,0',
                     'pAm': '192,192,0',
                     'Cm': '0,192,0',
                     'Gm': '192,0,192',
                     'Im': '0,192,192',
                     'pIm': '0,0,192',
                     'Ym': '0,0,64',
                     'pYm': '0,64,0',
                     'Um': '64,0,0',
                     'mcmo5Um': '0,205,139',
                     'pmcmo5Um': '64,64,64',
                     'Ar(p)': '64,64,0',
                     'Gr(p)': '64,0,64',
                     'pGr(p)': '0,64,64',
                     "N2'3'cp": '255,64,0',
                     'm2,8A': '255,0,64',
                     'pm2,8A': '64,255,0',
                     'msms2i6A': '0,255,64',
                     'pmsms2i6A': '64,0,255',
                     'ges2U': '0,64,255',
                     'pges2U': '255,64,64',
                     'k2C': '64,255,64',
                     'pk2C': '64,64,255',
                     'm2A': '205,205,205',
                     'pm2A': '139,0,0',
                     'ms2ct6A': '139,139,0',
                     'pms2ct6A': '0,139,0',
                     'ms2io6A': '139,0,139',
                     'pms2io6A': '0,139,139',
                     'ms2hn6A': '0,0,139',
                     'pms2hn6A': '0,0,70',
                     'pms2i6A': '0,70,0',
                     'ms2i6A': '205,0,0',
                     'ms2m6A': '139,205,0',
                     'pms2m6A': '205,139,0',
                     'ms2t6A': '205,0,139',
                     'pms2t6A': '0,205,0',
                     'se2U': '0,139,205',
                     's2Um': '255,0,0',
                     'ps2Um': '0,255,0',
                     's2C': '0,0,255',
                     'ps2C': '255,255,0',
                     's2U': '255,0,255',
                     'pm2G': '0,255,255',
                     'm3Um': '128,0,0',
                     'pm3Um': '0,128,0',
                     'acp3D': '0,0,128',
                     'pacp3D': '128,128,0',
                     'acp3Y': '128,0,128',
                     'pacp3Y': '0,128,128',
                     'acp3U': '255,128,0',
                     'pacp3U': '255,0,128',
                     'm3C': '128,255,0',
                     'pm3C': '0,255,128',
                     'm3Y': '128,0,255',
                     'pm3Y': '0,128,255',
                     'm3U': '255,128,128',
                     'pm3U': '128,255,128',
                     'imG-14': '128,128,255',
                     'pimG-14': '192,192,192',
                     'pm4C': '192,0,0',
                     's4U': '192,192,0',
                     'ps4U': '0,192,0',
                     'pm4Cm': '192,0,192',
                     'CoApN': '0,192,192',
                     'acCoApN': '0,0,192',
                     'malonyl-CoApN': '0,0,64',
                     'succinyl-CoApN': '0,64,0',
                     'ppN': '64,0,0',
                     "5'-OH-N": '0,205,139',
                     'NADpN': '64,64,64',
                     'pppN': '64,64,0',
                     'm5Cm': '64,0,64',
                     'pm5Cm': '0,64,64',
                     'm5Um': '255,64,0',
                     'pD': '255,0,64',
                     'mchm5Um': '64,255,0',
                     'pmchm5Um': '0,255,64',
                     'mchm5U': '64,0,255',
                     'pmchm5U': '0,64,255',
                     'pcmo5U': '255,64,64',
                     'phm5C': '64,255,64',
                     'inm5Um': '64,64,255',
                     'pinm5Um': '205,205,205',
                     'inm5s2U': '139,0,0',
                     'pinm5s2U': '139,139,0',
                     'inm5U': '0,139,0',
                     'pinm5U': '139,0,139',
                     'pmcm5s2U': '0,139,139',
                     'nm5ges2U': '0,0,139',
                     'pnm5ges2U': '0,0,70',
                     'nm5se2U': '0,70,0',
                     'pnm5se2U': '205,0,0',
                     'nm5s2U': '139,205,0',
                     'pnm5s2U': '205,139,0',
                     'nm5U': '205,0,139',
                     'pnm5U': '0,205,0',
                     'nchm5U': '0,139,205',
                     'pnchm5U': '255,0,0',
                     'ncm5Um': '0,255,0',
                     'pncm5Um': '0,0,255',
                     'ncm5s2U': '255,255,0',
                     'pncm5s2U': '255,0,255',
                     'ncm5U': '0,255,255',
                     'pncm5U': '128,0,0',
                     'chm5U': '0,128,0',
                     'pchm5U': '0,0,128',
                     'cm5s2U': '128,128,0',
                     'pcm5s2U': '128,0,128',
                     'cmnm5Um': '0,128,128',
                     'pcmnm5Um': '255,128,0',
                     'cmnm5ges2U': '255,0,128',
                     'pcmnm5ges2U': '128,255,0',
                     'cmnm5se2U': '0,255,128',
                     'pcmnm5se2U': '128,0,255',
                     'cmnm5s2U': '0,128,255',
                     'pcmnm5s2U': '255,128,128',
                     'cmnm5U': '128,255,128',
                     'pcmnm5U': '128,128,255',
                     'cm5U': '192,192,192',
                     'pcm5U': '192,0,0',
                     'cnm5U': '192,192,0',
                     'pcnm5U': '0,192,0',
                     'f5Cm': '192,0,192',
                     'pf5Cm': '0,192,192',
                     'f5C': '0,0,192',
                     'pf5C': '0,0,64',
                     'ho5C': '0,64,0',
                     'pho5C': '64,0,0',
                     'hm5C': '0,205,139',
                     'ho5U': '64,64,64',
                     'pho5U': '64,64,0',
                     'mcm5Um': '64,0,64',
                     'pmcm5Um': '0,64,64',
                     'mcm5s2U': '255,64,0',
                     'mcm5U': '255,0,64',
                     'pmcm5U': '64,255,0',
                     'mo5U': '0,255,64',
                     'pmo5U': '64,0,255',
                     'm5s2U': '0,64,255',
                     'pm5s2U': '255,64,64',
                     'mnm5ges2U': '64,255,64',
                     'pmnm5ges2U': '64,64,255',
                     'mnm5se2U': '205,205,205',
                     'pmnm5se2U': '139,0,0',
                     'mnm5s2U': '139,139,0',
                     'pmnm5s2U': '0,139,0',
                     'mnm5U': '139,0,139',
                     'm5C': '0,139,139',
                     'pm5C': '0,0,139',
                     'm5D': '0,0,70',
                     'pm5D': '0,70,0',
                     'm5U': '205,0,0',
                     'pm5U': '139,205,0',
                     'tm5s2U': '205,139,0',
                     'ptm5s2U': '205,0,139',
                     'tm5U': '0,205,0',
                     'ptm5U': '0,139,205',
                     'pm6,6A': '255,0,0',
                     'yW-86': '0,255,0',
                     'pyW-86': '0,0,255',
                     'yW-72': '255,255,0',
                     'yW-58': '255,0,255',
                     'pyW-58': '0,255,255',
                     'pyW-72': '128,0,0',
                     'preQ1base': '0,128,0',
                     'preQ1': '0,0,128',
                     'ppreQ1': '128,128,0',
                     'preQ0base': '128,0,128',
                     'preQ0': '0,128,128',
                     'ppreQ0': '255,128,0',
                     'm7G': '255,0,128',
                     'pm7G': '128,255,0',
                     'm8A': '0,255,128',
                     'pm8A': '128,0,255',
                     'pac4Cm': '0,128,255',
                     'pm4,4C': '255,128,128',
                     'A': '128,255,128',
                     'ApppppN': '128,128,255',
                     'AppppN': '192,192,192',
                     'ApppN': '192,0,0',
                     'pAp': '192,192,0',
                     'ppA': '0,192,0',
                     'pA': '192,0,192',
                     "pA2'3'cp": '0,192,192',
                     'pppA': '0,0,192',
                     'C+': '0,0,64',
                     'pC+': '0,64,0',
                     'mmpN': '64,0,0',
                     'mpN': '0,205,139',
                     'N': '64,64,64',
                     'G+': '64,64,0',
                     'pG+': '64,0,64',
                     'ct6A': '0,64,64',
                     'pct6A': '255,64,0',
                     'C': '255,0,64',
                     'pC': '64,255,0',
                     "pC2'3'cp": '0,255,64',
                     'D': '64,0,255',
                     'oQ': '0,64,255',
                     'poQ': '255,64,64',
                     'galQ': '64,255,64',
                     'pgalQ': '64,64,255',
                     'mpppN': '205,205,205',
                     'gluQ': '139,0,0',
                     'pgluQ': '139,139,0',
                     "pG2'3'cp": '0,139,0',
                     'G': '139,0,139',
                     'pG(pN)': '0,139,139',
                     'GpppN': '0,0,139',
                     'pGp': '0,0,70',
                     'ppG': '0,70,0',
                     'pG': '205,0,0',
                     'pppG': '139,205,0',
                     'ht6A': '205,139,0',
                     'pht6A': '205,0,139',
                     'OHyW': '0,205,0',
                     'pOHyW': '0,139,205',
                     'I': '255,0,0',
                     'pI': '0,255,0',
                     'imG2': '0,0,255',
                     'pimG2': '255,255,0',
                     'manQ': '255,0,255',
                     'pmanQ': '0,255,255',
                     'OHyWy': '128,0,0',
                     'pOHyWy': '0,128,0',
                     'mimG': '0,0,128',
                     'pmimG': '128,128,0',
                     'pm1A': '128,0,128',
                     'pac4C': '0,128,128',
                     'm2Gm': '255,128,0',
                     'pm2Gm': '255,0,128',
                     'm2,7Gm': '128,255,0',
                     'pm2,7Gm': '0,255,128',
                     'm2,7G': '128,0,255',
                     'm2,7GpppN': '0,128,255',
                     'pm2,7G': '255,128,128',
                     'm2,2Gm': '128,255,128',
                     'pm2,2Gm': '128,128,255',
                     'm2,2,7G': '192,192,192',
                     'm2,2,7GpppN': '192,0,0',
                     'pm2,2,7G': '192,192,0',
                     'm2,2G': '0,192,0',
                     'pm2,2G': '192,0,192',
                     'm2G': '0,192,192',
                     'm4Cm': '0,0,192',
                     'm4,4Cm': '0,0,64',
                     'pm4,4Cm': '0,64,0',
                     'm4,4C': '64,0,0',
                     'ac4Cm': '0,205,139',
                     'ac4C': '64,64,64',
                     'm4C': '64,64,0',
                     'm6Am': '64,0,64',
                     'pm6Am': '0,64,64',
                     'm6,6Am': '255,64,0',
                     'pm6,6Am': '255,0,64',
                     'm6,6A': '64,255,0',
                     'io6A': '0,255,64',
                     'pio6A': '64,0,255',
                     'ac6A': '0,64,255',
                     'pac6A': '255,64,64',
                     'f6A': '64,255,64',
                     'pf6A': '64,64,255',
                     'g6A': '205,205,205',
                     'pg6A': '139,0,0',
                     'hm6A': '139,139,0',
                     'phm6A': '0,139,0',
                     'hn6A': '139,0,139',
                     'phn6A': '0,139,139',
                     'pi6A': '0,0,139',
                     'i6A': '0,0,70',
                     'm6ApppppN': '0,70,0',
                     'm6AppppN': '205,0,0',
                     'm6ApppN': '139,205,0',
                     'm6t6A': '205,139,0',
                     'pm6t6A': '205,0,139',
                     'm6A': '0,205,0',
                     'pm6A': '0,139,205',
                     't6A': '255,0,0',
                     'pt6A': '0,255,0',
                     'm7GppppN': '0,0,255',
                     'm7GpppN': '255,255,0',
                     'pCm': '255,0,255',
                     'pGm': '0,255,255',
                     'pUm': '128,0,0',
                     'o2yW': '0,128,0',
                     'po2yW': '0,0,128',
                     'Y': '128,128,0',
                     'pY': '128,0,128',
                     'Qbase': '0,128,128',
                     'Q': '255,128,0',
                     'pQ': '255,0,128',
                     'OHyWx': '128,255,0',
                     'pOHyWx': '0,255,128',
                     'pN': '128,0,255',
                     'xX': '0,128,255',
                     'xA': '255,128,128',
                     'xC': '128,255,128',
                     'xG': '128,128,255',
                     'xU': '192,192,192',
                     'Xm': '192,0,0',
                     'U': '192,192,0',
                     'cmo5U': '0,192,0',
                     'mcmo5U': '192,0,192',
                     'pmcmo5U': '0,192,192',
                     'pU': '0,0,192',
                     "pU2'3'cp": '0,0,64',
                     'yW': '0,64,0',
                     'pyW': '64,0,0',
                     'pyyW': '0,205,139',
                     'imG': '64,64,64',
                     'pimG': '64,64,0'}
    if not rgb_colors.get(modi):
        print("Please check your modification name. It does not seem to be a valid MODOMICS shortname."
              "No RGB values were created.")
        return None
    else:
        return rgb_colors.get(modi)


def parse_excel_sheetnames(input_file):
    """
    parses the input excel file and returns a list of sheetnames.  
    This is useful if the Excel file contains multiple sheets with information to be converted into bedrmod.
    Not really efficient, if the input file is huge, though. 
    """
    file = pd.read_excel(input_file, None)
    return file.keys()