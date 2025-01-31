from pathlib import Path

from bedRMod.helper import write_bioinformatics_keys, read_bioinformatics_keys
from bedRMod.helper import funcify

test_dir = Path(__file__).parent

def test_bioinformatics_keys():
    workflow = "https://github.com/anmabu/bedRMod"
    score_func_str = "round(-log10(score))"
    coverage_func_str = "round(coverage)"
    frequency_func_str = "round(frequency)"
    write_bioinformatics_keys(f"{test_dir}/test_config.yaml", workflow=workflow, score_function=score_func_str,
                              coverage_function=coverage_func_str, frequency_function=frequency_func_str)
    workflow, cov_f, freq_f, score_f = read_bioinformatics_keys(f"{test_dir}/test_config.yaml")
    assert workflow == workflow
    assert score_f == score_func_str
    assert cov_f == coverage_func_str
    assert freq_f == frequency_func_str


def test_funcify():
    func = funcify("lambda x: round(-log10(x))")
    print(func(0.1))

