#! /usr/bin/env python3

"""Data wrangling wrapper to convert generic text/flat files
to EU format for input into Sci-ModoM. This script also generates
track hub components.
"""

import sys
import logging
import argparse
import json
import yaml

import numpy as np
import pandas as pd

from pandas.api.types import is_numeric_dtype

from pathlib import Path

# TODO: temporary
from paths import PROJECT_ROOT
UTILS_PATH = PROJECT_ROOT # / 'utils'

import utils as utils

logger = logging.getLogger(__name__)

# TODO: define in package
EUF_VERSION = "bedModv1.2"

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

db_json_keys = [
    "id",
    "title",
    "experiment", # TO BE REMOVED
    "summary", # TO BE REMOVED
    "description",
    "organism", 
    "assembly",
    "annotation_source",
    "annotation_version",
    "methods",
    "cell_tissue",
    "references",
    "tags",
]

euf_fields = [
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "coverage",
    "frequency",
    "refBase",
]
    
fmt_keys = [
    "delimiter",
    "excel",
    "comment",
    "skiprows",
    "freq_pc",
    "coords_idx",
    "score_transformation",
]
     

mod_colours = {
    "m6A": "255,0,0",
    "m5C": "0,255,0",
    "Y": "0,0,255",
}


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Data wrangling wrapper to convert generic 
                                     text/flat files to EU format for input into Sci-ModoM.""")

    parser.add_argument('config', help="""The yml configuration file.""")
    
    parser.add_argument('--header', help="""A list of fields to use as header if no header
                        is present. Fields must match those in the config.""", nargs="+", type=str)
    
    parser.add_argument('--split', help="""A list of two items 'field sep' to split the 
                        column 'field' formatted as chrom/sep/pos/sep/strand. The resulting
                        fields are taken from the config.""", nargs=2, type=str)
    
    utils.add_logging_options(parser)
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    utils.update_logging(args)
    
    msg = f'[flat2euf]: {" ".join(sys.argv)}'
    logger.info(msg)

    # read configuration file
    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    
    required_keys = ["destination_output", "input_files", 
                     "meta_options", "fmt_options", 
                     "file_meta_options", "file_fmt_options"]
    utils.check_keys_exist(config, required_keys)
    
    # create output directory
    euf_dir = Path(config["destination_output"], 'euf')
    euf_dir.mkdir(parents=True, exist_ok=False)
    
    # metadata
    euf_header = dict()
    for key in euf_header_keys:
        euf_header[key] = config["meta_options"].get(key, None)
    euf_header["fileformat"] = EUF_VERSION
    
    db_json = dict()
    for key in db_json_keys:
        db_json[key] = config["meta_options"].get(key, None)

    # formatting
    fmt_options = config["fmt_options"]
    
    # process
    for file_id, file_path in config["input_files"].items():
        
        msg = f"Converting {file_path} to EUF..."
        logger.info(msg)
        
        # check that all EUF fields and formatting options are present
        # if not, will throw an error
        file_fmt_options = dict()
        if config["file_fmt_options"] is not None:
            file_fmt_options = config["file_fmt_options"].get(file_id, dict())
        fmt = {k:file_fmt_options.get(k, fmt_options[k]) for k in fmt_keys + euf_fields}
        # currently default to chromStart/End 
        for k in ["thickStart", "thickEnd"]:
            fmt[k] = None
        
        kwargs = {
            "comment": fmt["comment"],
            "skiprows": fmt["skiprows"]
        }
        if fmt["excel"]:
            reader = pd.read_excel
        else:
            reader = pd.read_csv
            kwargs["sep"] = fmt['delimiter']
        if args.header:
            kwargs["header"] = None
            kwargs["names"] = args.header
        df = reader(file_path, **kwargs)
        
        # split if required
        if args.split is not None:
            df[[fmt["chrom"], fmt["chromStart"], fmt["strand"]]] = df[args.split[0]].str.split(
                args.split[1], 
                expand=True
            )
            # make sure chromStart is treated as number
            df[fmt["chromStart"]] = df[fmt["chromStart"]].astype(int)
        
        # we need to perform some checks before defining the fields
        columns = df.columns
        columns = [fmt[f] if fmt[f] in columns else i for i, f in enumerate(euf_fields)]
        # if chromStart == chromEnd 
        if columns[1] == columns[2]:
            columns[2] = 2
        
        euf = df[[c for c in columns if not isinstance(c, int)]].copy()
        # adjust for indexing of chromStart/End 
        # assume that valid fields are open at the End, else 
        # adjust for single-site (chromStart == chromEnd)
        values = euf[fmt["chromStart"]].values
        if fmt["coords_idx"] == 1:
            euf[fmt["chromStart"]] = values - 1
            
            msg = "Subtracting 1 from chromStart..."
            logger.debug(msg)
            
        if fmt["chromStart"] == fmt["chromEnd"]:
            fmt["chromEnd"] = "chromEnd"
            columns[2] = "chromEnd"
            if fmt["coords_idx"] == 1:
                euf["chromEnd"] = values
            else:
                euf["chromEnd"] = values + 1
                
            msg = "Adding chromEnd..."
            logger.debug(msg)
            
        for i, k1, k2 in zip([6,7], ["thickStart", "thickEnd"], [fmt["chromStart"], fmt["chromEnd"]]):
            fmt[k1] = k1
            columns[i] = k1
            euf[k1] = euf[k2].values
        
        # if itemRGB is undefined, add colour 
        if isinstance(columns[8], int):
            columns[8] = "itemRGB"
            fmt["itemRGB"] = "itemRGB"
            colour = mod_colours[fmt["name"]]
            euf["itemRGB"] = colour
            
            msg = f"Setting itemRGB to {colour}..."
            logger.debug(msg)
            
        # for remaining fields, process as is, report "non-valid" headers
        for c in columns:
            if not isinstance(c, int):
                continue
            h = euf_fields[c]
            columns[c] = h
            euf[h] = fmt[h]
            fmt[h] = h
            
            msg = f"Filling {h} with value from config: {fmt[h]}..."
            logger.debug(msg)
            
        if any([isinstance(c, int) for c in columns]):
            msg = f"Some header fields {','.join([str(c) for c in columns])} " \
                  f"remained unprocessed..."
            logger.critical(msg)
            
        euf = euf[columns]
        euf.columns = euf_fields
        
        # now adjust for score and frequency
        if fmt["score_transformation"] == "pvalue":
            # convert p-values to -log10 (use maximum real for p-values -> 0)
            # then scale to [0-1000] - transformation is data-dependent
            values = -np.log10(euf['score'])
            values[values == np.inf] = np.ma.masked_invalid(values).max()
            euf['score'] = np.round(1000*(values-min(values))/(max(values)-min(values)))
            euf['score'] = euf['score'].astype(int)
        elif fmt["score_transformation"] == "probability":
            euf['score'] = np.round(1000*euf['score'].values)
            euf['score'] = euf['score'].astype(int)
        elif fmt["score_transformation"] == "no":
            pass
        else:
            msg = f"score_transformation: {fmt['score_transformation']} " \
                  f"undefined..."
            logger.critical(msg)
        
        if not fmt["freq_pc"] and is_numeric_dtype(euf['frequency']):
            euf['frequency'] = euf['frequency']*100
            
        # sort on the chrom field, and then on the chromStart field.
        euf.sort_values(["chrom", "chromStart"], ascending=[True, True], inplace=True)
        
        # EUF header and metadata
        file_meta_options = dict()
        if config["file_meta_options"] is not None:
            file_meta_options = config["file_meta_options"].get(file_id, dict())
            
        header = {k:file_meta_options.get(k, euf_header[k]) for k in euf_header_keys}
        filen = Path(euf_dir, f"{file_id}.bed")
        if filen.is_file():
            msg = f"File {filen} already exists... removing..."
            logger.warning(msg)
            filen.unlink()
        
        msg = f"... writing to {filen} (incl. metadata)"
        logger.info(msg)
        
        with open(filen, 'a') as fh:
            for k, v in header.items():
                fh.write(f"#{k}={v}\n")
            header = "\t".join(euf.columns)
            fh.write(f"#{header}\n")
            euf.to_csv(
                fh, 
                sep="\t",
                header=False,
                index=False)
        
        metadata = {k:file_meta_options.get(k, db_json[k]) for k in db_json_keys}
        metadata["id"] = file_id
        
        # TODO: these keys will be removed
        metadata["experiment"] = metadata["title"]
        metadata["summary"] = metadata["title"]
        
        filen = Path(euf_dir, f"{file_id}.json")
        with open(filen, "w") as fh:
            # add [ ] to conform to current rmap_web version...
            json.dump([metadata], fh, indent=2)


if __name__ == '__main__':
    main()
    
