from tsv2bedRMod.helper import write_bioinformatics_keys, read_bioinformatics_keys
from tsv2bedRMod.helper import funcify


def test_bioinformatics_keys():
    score_func_str = "round(-log10(score))"
    coverage_func_str = "round(coverage)"
    frequency_func_str = "round(frequency)"
    write_bioinformatics_keys("test_config.yaml", score_function=score_func_str,
                              coverage_function=coverage_func_str, frequency_function=frequency_func_str)
    score_f, cov_f, freq_f = read_bioinformatics_keys("test_config.yaml")
    assert score_f == score_func_str
    assert cov_f == coverage_func_str
    assert freq_f == frequency_func_str


def test_funcify():
    func = funcify("lambda x: round(-log10(x))")
    print(func(0.1))

