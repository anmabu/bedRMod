import pytest
import pandas as pd

from tsv2bedRMod.tsv2bedRMod import csv2bedRMod, df2bedRMod, parse_row

# def testdf2bedrmod():
def test_parse_row():
    columns = ['chrom', 'start_col', 'end', 'name', 'score_column', 'strandedness', 
           'thick_start', 'thick_end', 'item_rgb', 'coverage_col', 'frequency_col']
    data = [
        ['chr1', 1000, 1001, 'm1A', 900, '+', 1000, 1001, '0,128,128', 30, 90],
        #['chr2', 2000, 2500, 'gene2', 850, '-', 2000, 2500, '0,255,0', 45, 0.8],
        #['chr3', 3000, 3500, 'gene3', 700, '+', 3000, 3500, '0,0,255', 60, 0.9],
        #['chr4', 4000, 4500, 'gene4', 920, '-', 4000, 4500, '255,255,0', 50, 0.7]
            ]
    df = pd.DataFrame(data, columns=columns)
    row = df.head(1).iloc[0]
    # print(row)
    result = parse_row(row, columnnames=df.columns, ref_seg="chrom", start="start_col", start_function=None, modi="name", modi_column=True, score="score_column", score_function=None, strand="strandedness", coverage="coverage_col", coverage_function=None, frequency="frequency_col", frequency_function=None)
    assert result == ("1", 1000, 1001, "m1A", 900, "+", 1000, 1001, "0,128,128", 30, 90)