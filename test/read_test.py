import pandas as pd
from pathlib import Path

from bedRMod import read_header, read_data, read_bedRMod

test_dir = Path(__file__).parent


def test_read_header():
    test_file = f"{test_dir}/test_df2bedrmod.bedrmod"
    header = read_header(test_file)
    assert header == {'fileformat': 'bedRModv1.8', 'organism': '9606', 'modification_type': 'RNA', 'assembly': 'GRCh38', 'annotation_source': 'Ensemble', 'annotation_version': '93', 'sequencing_platform': 'Illumina NovaSeq 6000', 'basecalling': '', 'bioinformatics_workflow': 'workflow:https://github.com/anmabu/bedRMod;coverage_function:round(coverage);frequency_function:round(frequency);score_function:round(-log10(score))', 'experiment': 'https://doi.org/test', 'external_source': 'GEO;GSETEST', 'methods': 'TEST', 'references': 'pubmed_id:12345678', 'conversion_information': 'rgb value in df is ignored and calculated directly'}


def test_read_data():
    test_file = f"{test_dir}/test_df2bedrmod.bedrmod"
    columns = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "coverage", "frequency"]
    data = [
        ['1', 1000, 1001, 'm1A', 30, '+', 1000, 1001, '0,128,128', 30, 90],
        ['2', 1001, 1002, 'm5C', 19, '-', 1001, 1002, '0,139,139', 45, 83],
        ['X', 3001, 3002, 'm1A', 117, '+', 3001, 3002, '0,128,128', 6, 79],
        ['Y', 4001, 4002, 'm3C', 18, '-', 4001, 4002, '128,255,0', 50, 37],
        ['MT', 5001, 5002, 'm3C', 40, '-', 5001, 5002, '128,255,0', 23, 65]
    ]

    compare_df = pd.DataFrame(data, columns=columns)
    
    data = read_data(test_file)
    df_diff = data.compare(compare_df)
    assert df_diff.empty

