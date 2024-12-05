import os
import pandas as pd

from ruamel.yaml import YAML

yaml = YAML()
yaml.sort_base_mapping_type_on_output = False  # disable sorting of keys

from bedRMod.helper import write_config_header, get_modification_color, check_value_range


def parse_row(row, columnnames=[], ref_seg="ref_seg", start="pos", start_function=None, modi="m1A", modi_column=False,
              score=None, score_function=None, strand="strand", coverage=None, coverage_function=None, frequency=None,
              frequency_function=None):
    """
    parses a dataframe/csv row and return the values needed for a row in the bedRMod format
    """
    chrom = row[ref_seg]
    has_alpha = any(c.isalpha() for c in chrom)
    has_digit = any(c.isdigit() for c in chrom)
    if chrom == "chrY" or (chrom == "Y"):
        chrom = "Y"
    elif chrom == "chrX" or (chrom == "X"):
        chrom = "X"
    elif chrom == "chrMT" or (chrom == "MT") or (chrom == "M") or (chrom == "chrM"):
        chrom = "MT"    
    elif has_alpha and has_digit:
        chrom = ''.join(c for c in chrom if c.isdigit())
    elif has_digit and not has_alpha:
        chrom = chrom
    else: 
        print(f"something is weird in chrom {chrom}") 
    if start_function is not None:
        if type(start) == list:
            params = [row[col] for col in start]
        elif isinstance(start, str):
            params = row[start]
        else:
            params = start
        start_col = start_function(params)
    else:
        start_col = int(row[start])
    if start_col is None:
        return None
    end = start_col + 1
    name = row[modi] if modi_column else modi
    if score_function is not None:
        if type(score) == list:
            params = [row[col] for col in score]
        elif isinstance(score, str):
            params = row[score]
        else:
            params = score
        score_column = score_function(params)
    else:
        if isinstance(score, str):
            score_column = round(row[score])
        else:
            score_column = score
    if strand == "+":
        strandedness = "+"
    elif strand == "-":
        strandedness = "-"
    else:
        strandedness = row[strand]
    if coverage_function is not None:
        if type(coverage) == list:
            params = [row[col] for col in coverage]
        elif isinstance(coverage, str):
            params = row[coverage]
        coverage_col = coverage_function(params)
    else:
        if coverage in columnnames:
            coverage_col = round(row[coverage])
        else:
            coverage_col = coverage
    if frequency_function is not None:
        if type(frequency) == list:
            params = [row[col] for col in frequency]
        elif isinstance(frequency, str):
            params = row[frequency]
        frequency_col = frequency_function(params)
    else:
        if frequency in columnnames:
            frequency_col = round(row[frequency])
        elif isinstance(frequency, (int, float)):
            frequency_col = round(frequency)
    thick_start = start_col
    thick_end = end
    item_rgb = get_modification_color(name)
    result = (chrom, start_col, end, name, score_column, strandedness, thick_start, thick_end, item_rgb, coverage_col,
            frequency_col)
    bedrmod_columns = ("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd",
                       "itemRgb", "coverage", "frequency")
    for index, item in enumerate(result):
        if item is None:
            print(f"The data has not been converted. \n"
                  f"Please check the input value/function for the {bedrmod_columns[index]} column.")
            result = None
    return result


def csv2bedRMod(input_file, config_yaml, output_file=None, delimiter=None, ref_seg="ref_seg", start="pos",
                start_function=None, modi="m1A", modi_column=False, score=None, score_function=None, strand="strand",
                coverage=None, coverage_function=None, frequency=None, frequency_function=None):

    """
    converts arbitrary csv files into bedRMod format.
    The parameters usually pass the column name of the csv which contains the respective information.
    The name of the output file is infered from the input file and put in the same directory as the input file.
    :param input_file:(path to) input csv file.
    :param config_yaml: (path to) config file containing the information on the metadata
    :param output_file: (path to) output bedrmod file. If "None" it is the path to input file
    :param delimiter: delimiter of the passed csv file. If "None" it is infered by pandas.
    :param ref_seg: column name of the column containing the reference sequence. i.e. the chromosome
    :param start: column name of the column that contains the positions of the modification
    :param start_function: fix value of column eg. off-by-one errors
    :param modi: contains the column name of the column containing the modification or the name of the
    modification for the whole file.
    :param modi_column: indicates whether the value passed to "modi" contains the column name containing the
    modification (True) or denominates the modifiation itself (False)
    :param score: can either be a fixed value e.g. 0 or 1000 if score is totally unknown, or can be a calculation of the
    score e.g. int(1000 - (row["FDR"] * 1000)) or indicate a column name containing the score.
    :param score_function:
    :param strand: indicates the column name of the column containing the strandedness. can be "+" or "-" to indicate
    same strandedness for whole file.
    :param coverage: column name of column containing the coverage at this position.
    :param coverage_function:
    :param frequency:
    :param frequency_function:
    :return:
    """
    file = pd.read_csv(input_file, delimiter=delimiter)
    if output_file is None:
        output_file = input_file

    path, ending = os.path.splitext(output_file)
    if not ending == ".bedrmod":
        output_file = path + ".bedrmod"
        print(f"output file: {output_file}")

    #config = yaml.safe_load(open(config_yaml, "r"))
    config = yaml.load(open(config_yaml, "r"))

    colnames = file.columns
    try:
        with open(output_file, 'w') as f:
            header_written = write_config_header(config, f)
            if not header_written:
                raise TypeError("Header could not be written.")
            f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                    "\tfrequency\n")

            for _, row in file.iterrows():
                result = parse_row(row, colnames, ref_seg, start, start_function, modi, modi_column, score,
                                   score_function,
                                   strand, coverage, coverage_function, frequency, frequency_function)
                if not any(item is None for item in result) or (result is not None):
                    chrom, start_col, end, name, score_column, strandedness, thick_start, thick_end, item_rgb, \
                        coverage_col, frequency_col = result
                    f.write(f'{chrom}\t{start_col}\t{end}\t{name}\t{score_column}\t{strandedness}\t{thick_start}'
                            f'\t{thick_end}\t{item_rgb}\t{coverage_col}\t{frequency_col}\n')
            print("Done!")
    except TypeError:
        os.remove(output_file)

            
def df2bedRMod(df, config_yaml, output_file, ref_seg="ref_seg", start="pos", start_function=None, modi="m1A",
               modi_column=False, score=None, score_function=None, strand="strand", coverage=None,
               coverage_function=None, frequency=None, frequency_function=None, from_gui=False):

    """
    converts arbitrary pandas_dataframes into bedRMod format.
    The parameters usually pass the column name of the csv which contains the respective information.
    :param df: input pandas dataframe.
    :param config_yaml: (path to) config file containing the information on the metadata
    :param output_file:
    :param ref_seg: column name of the column containing the reference sequence. i.e. the chromosome
    :param start: column name of the column that contains the positions of the modification
    :param start_function: fix value of column eg. off-by-one errors
    :param modi: contains the column name of the column containing the modification or the name of the
    modification for the whole df.
    :param modi_column: indicates whether the value passed to "modi" contains the column name containing the
    modification (True) or denominates the modifiation itself (False)
    :param score: can either be a fixed value e.g. 0 or 1000 if score is totally unknown, or can be a calculation of the
    score e.g. int(1000 - (row["FDR"] * 1000)) or indicate a column name containing the score.
    :param score_function:
    :param strand: indicates the column name of the column containing the strandedness. can be "+" or "-" to indicate
    same strandedness for whole df.
    :param coverage: column name of column containing the coverage at this position.
    :param coverage_function:
    :param frequency:
    :param frequency_function:
    :param from_gui: (bool) indicates whether it was called from the GUI.
    :return:
    """

    path, ending = os.path.splitext(output_file)
    if not ending == ".bedrmod":
        output_file = path + ".bedrmod"
        print(f"output file: {output_file}")

    # config = yaml.safe_load(open(config_yaml, "r"))
    config = yaml.load(open(config_yaml, "r"))

    colnames = df.columns
    try:
        with open(output_file, 'w') as f:
            header_written = write_config_header(config, f)
            if not header_written:
                raise TypeError("Header could not be written.")
            f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                    "\tfrequency\n")

            for _, row in df.iterrows():
                result = parse_row(row, colnames, ref_seg, start, start_function, modi, modi_column, score, score_function,
                                   strand, coverage, coverage_function, frequency, frequency_function)
                if not any(item is None for item in result) or (result is not None):
                    chrom, start_col, end, name, score_column, strandedness, thick_start, thick_end, item_rgb, \
                        coverage_col, frequency_col = result
                    check_value_range(result)
                    f.write(f'{chrom}\t{start_col}\t{end}\t{name}\t{score_column}\t{strandedness}\t{thick_start}'
                            f'\t{thick_end}\t{item_rgb}\t{coverage_col}\t{frequency_col}\n')
            print("Done!")
            if from_gui:
                return "Done"
    except TypeError:
        print("Something went wrong! The bedrmod file was not generated!")
        os.remove(output_file)
        if from_gui:
            return "Error"
