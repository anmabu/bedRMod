__version__ = "1.8.0-rc3"

from .convert2bedRMod import df2bedRMod, csv2bedRMod
from .helper import parse_row
from .read import read_header, read_data, read_bedRMod
from .write import write_header_from_dict, write_header_from_config,  write_data_from_df, write_bedRMod
