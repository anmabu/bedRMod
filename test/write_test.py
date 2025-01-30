import pandas as pd
from pathlib import Path

from bedRMod import write_header, write_data
from bedRMod import read_header, read_data

test_dir = Path(__file__).parent

def test_write_header():
    header_dict = read_header(f"{test_dir}/test_df2bedrmod.bedrmod")
    write_header(f"{test_dir}/test_write2bedrmod.bedrmod", header_dict)

    assert header_dict == read_header(f"{test_dir}/test_write2bedrmod.bedrmod")

def test_write_data():
    data_df = read_data(f"{test_dir}/test_df2bedrmod.bedrmod")
    write_data(f"{test_dir}/test_write2bedrmod.bedrmod", data_df)
    diff_df = data_df.compare(read_data(f"{test_dir}/test_write2bedrmod.bedrmod"))
    assert diff_df.empty
