__version__ = "1.8.0-rc3"

from .convert2bedRMod import df2bedRMod, csv2bedRMod, parse_row
from .read import read_header, read_data, read_bedRMod
from .write import write_header, write_data, write_bedRMod
