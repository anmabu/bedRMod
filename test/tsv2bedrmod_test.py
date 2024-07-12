import pytest
import pandas as pd
import filecmp

from flat2euf.utils import check_keys_exist
from tsv2bedRMod.tsv2bedRMod import csv2bedRMod, df2bedRMod, parse_row
from tsv2bedRMod.helper import write_bioinformatics_keys


def test_parse_row_working_example():
    columns = ['chrom', 'start_col', 'end', 'name', 'score_column', 'strandedness', 
           'thick_start', 'thick_end', 'item_rgb', 'coverage_col', 'frequency_col']
    data = [['chr1', 1000, 1001, 'm3C', 900, '+', 1000, 1001, '0,128,128', 30, 90]]
    df = pd.DataFrame(data, columns=columns)
    row = df.head(1).iloc[0]
    result = parse_row(row, columnnames=df.columns, ref_seg="chrom", start="start_col", modi="m1A", modi_column=False,
                       score="score_column", strand="strandedness", coverage="coverage_col",
                       frequency="frequency_col")
    assert result == ("1", 1000, 1001, "m1A", 900, "+", 1000, 1001, "0,128,128", 30, 90)


def test_parse_row_freq_none():
    def freq_func(param):
        return None
    columns = ['chrom', 'start_col', 'end', 'name', 'score_column', 'strandedness', 'thick_start', 'thick_end',
               'item_rgb', 'coverage_col', 'frequency_col']
    data = [['chr1', 1000, 1001, 'm1A', 900, '+', 1000, 1001, '0', 30, 90]]
    df = pd.DataFrame(data, columns=columns)
    row = df.head(1).iloc[0]
    result = parse_row(row, columnnames=df.columns, ref_seg="chrom", start="start_col", modi="name", modi_column=True,
                       score="score_column", strand="strandedness", coverage="coverage_col", frequency="frequency_col",
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

    def score_func(params):
        score, cov = params
        return round(score / cov)

    df = pd.DataFrame(data, columns=columns)
    df2bedRMod(df, "test_config.yaml", "test_df2bedrmod.bedrmod", ref_seg="chrom", start="start_col",
               modi="name", modi_column=True, score=["score_column", "coverage_col"], score_function=score_func,
               strand="strandedness", coverage="coverage_col", frequency="frequency_col")
    assert filecmp.cmp("test_static_df2bedrmod.bedrmod", "test_df2bedrmod.bedrmod")


def test_csv2bedrmod():
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
    df.to_csv("test_csv2bedrmod.csv")

    def cov_func(param):
        return param * 2

    def start_func(param):
        return param - 1
    csv2bedRMod("test_csv2bedrmod.csv", "test_config.yaml", ref_seg="chrom", start="start_col",
                start_function=start_func, modi="name", modi_column=True, score="score_column",
                strand="strandedness", coverage="coverage_col", coverage_function=cov_func,
                frequency="frequency_col")
    assert filecmp.cmp("test_static_csv2bedrmod.bedrmod", "test_csv2bedrmod.bedrmod")


def test_check_bioinformatics_keys():
    score_func_str = "round(-log10(score))"
    coverage_func_str = "round(coverage)"
    frequency_func_str = "round(frequency)"
    write_bioinformatics_keys("test_config.yaml", score_function=score_func_str,
                              coverage_function=coverage_func_str, frequency_function=frequency_func_str)
