import pandas as pd
from pathlib import Path
from ruamel.yaml import YAML

yaml = YAML()
yaml.sort_base_mapping_type_on_output = False  # disable sorting of keys

from bedRMod import write_header_from_dict, write_header_from_config, write_data_from_df, write_single_data_row, \
    read_single_data_row
from bedRMod import read_header, read_data_to_df

test_dir = Path(__file__).parent

def test_write_header_from_dict():
    header_dict = read_header(f"{test_dir}/test_static_df2bedrmod.bedrmod")
    write_header_from_dict(header_dict, f"{test_dir}/test_write2bedrmod.bedrmod")
    assert header_dict == read_header(f"{test_dir}/test_write2bedrmod.bedrmod")

def test_write_data_from_df():
    data_df = read_data_to_df(f"{test_dir}/test_static_df2bedrmod.bedrmod")
    write_data_from_df(data_df, f"{test_dir}/test_write2bedrmod.bedrmod")
    diff_df = data_df.compare(read_data_to_df(f"{test_dir}/test_write2bedrmod.bedrmod"))
    assert diff_df.empty

def test_write_header_from_config():
    write_header_from_config(f"{test_dir}/test_config.yaml", f"{test_dir}/test_output_config.bedrmod")
    assert read_header(f"{test_dir}/test_write2bedrmod.bedrmod") == read_header(f"{test_dir}/test_output_config.bedrmod")

def test_write_single_data_row():
    ref_seg = "1"
    start = 1000
    end = 1001
    name = "m1A"
    score = 30
    strand = "+"
    thick_start = 1000
    thick_end = 1001
    itemRgb = "0,128,128"
    coverage = 30
    frequency = 90

    # write a complete bedrmod file for testing
    write_header_from_config(f"{test_dir}/test_config.yaml", f"{test_dir}/test_output_single_data_row.bedrmod")
    write_single_data_row(f"{test_dir}/test_output_single_data_row.bedrmod", ref_seg, start, name, score, strand, coverage, frequency)
    # for the comparison the values have to be strings
    assert (read_single_data_row(f"{test_dir}/test_output_single_data_row.bedrmod", 0)
            == ['1', '1000', '1001', 'm1A', '30', '+', '1000', '1001', '0,128,128', '30', '90'])
