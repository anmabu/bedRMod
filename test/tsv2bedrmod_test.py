import pytest
import pandas as pd
import filecmp

from tsv2bedRMod.tsv2bedRMod import csv2bedRMod, df2bedRMod, parse_row


def test_parse_row_working_example():
    columns = ['chrom', 'start_col', 'end', 'name', 'score_column', 'strandedness', 
           'thick_start', 'thick_end', 'item_rgb', 'coverage_col', 'frequency_col']
    data = [['chr1', 1000, 1001, 'm1A', 900, '+', 1000, 1001, '0,128,128', 30, 90]]
    df = pd.DataFrame(data, columns=columns)
    row = df.head(1).iloc[0]
    result = parse_row(row, columnnames=df.columns, ref_seg="chrom", start="start_col", start_function=None, modi="name", modi_column=True, score="score_column", score_function=None, strand="strandedness", coverage="coverage_col", coverage_function=None, frequency="frequency_col", frequency_function=None)
    assert result == ("1", 1000, 1001, "m1A", 900, "+", 1000, 1001, "0,128,128", 30, 90)


def test_parse_row_freq_none():
    def freq_func(param):
        return None
    columns = ['chrom', 'start_col', 'end', 'name', 'score_column', 'strandedness', 'thick_start', 'thick_end',
               'item_rgb', 'coverage_col', 'frequency_col']
    data = [['chr1', 1000, 1001, 'm1A', 900, '+', 1000, 1001, '0', 30, 90]]
    df = pd.DataFrame(data, columns=columns)
    row = df.head(1).iloc[0]
    result = parse_row(row, columnnames=df.columns, ref_seg="chrom", start="start_col", start_function=None,
                       modi="name", modi_column=True, score="score_column", score_function=None, strand="strandedness",
                       coverage="coverage_col", coverage_function=None, frequency="frequency_col",
                       frequency_function=freq_func)
    assert result is None


def test_df2bedrmod():
    columns = ['chrom', 'start_col', 'end', 'name', 'score_column', 'strandedness',
               'thick_start', 'thick_end', 'item_rgb', 'coverage_col', 'frequency_col']
    data = [
        ['chr1', 1000, 1001, 'm1A', 900, '+', 1000, 1001, '0', 30, 90],
        ['2', 1001, 1002, 'm5C', 850, '-', 1001, 1002, '0', 45, 83],
        ['chrX', 3001, 3002, 'm1A', 700, '+', 3001, 3002, '0', 6, 79],
        ['chrY', 4001, 4002, 'm3C', 920, '-', 4001, 4002, '0', 50, 37],
        ['chrM', 5001, 5002, 'm3C', 920, '-', 5001, 5002, '0', 23, 65]
    ]

    df = pd.DataFrame(data, columns=columns)
    df2bedRMod(df, "test_config.yaml", "test_df2bedrmod.bedrmod", ref_seg="chrom", start="start_col",
               start_function=None, modi="name", modi_column=True, score="score_column", score_function=None,
               strand="strandedness", coverage="coverage_col", coverage_function=None,
               frequency="frequency_col", frequency_function=None)
    assert filecmp.cmp("test_static_df2bedrmod.bedrmod", "test_df2bedrmod.bedrmod")
