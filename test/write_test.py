import pandas as pd
from pathlib import Path
from ruamel.yaml import YAML

yaml = YAML()
yaml.sort_base_mapping_type_on_output = False  # disable sorting of keys

from bedRMod import write_header_from_dict, write_header_from_config, write_data_from_df, write_single_data_row
from bedRMod import read_header, read_data

test_dir = Path(__file__).parent

def test_write_header_from_dict():
    header_dict = read_header(f"{test_dir}/test_df2bedrmod.bedrmod")
    write_header_from_dict(header_dict, f"{test_dir}/test_write2bedrmod.bedrmod")
    assert header_dict == read_header(f"{test_dir}/test_write2bedrmod.bedrmod")

def test_write_data_from_df():
    data_df = read_data(f"{test_dir}/test_df2bedrmod.bedrmod")
    write_data_from_df(data_df, f"{test_dir}/test_write2bedrmod.bedrmod")
    diff_df = data_df.compare(read_data(f"{test_dir}/test_write2bedrmod.bedrmod"))
    assert diff_df.empty

def test_write_header_from_config():
    write_header_from_config(f"{test_dir}/test_config.yaml", f"{test_dir}/test_output_config.bedrmod")
    assert read_header(f"{test_dir}/test_write2bedrmod.bedrmod") == read_header(f"{test_dir}/test_output_config.bedrmod")

def test_write_single_data_row():
    pass