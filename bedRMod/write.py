import pandas as pd

def write_header(file, header_dict):
    """
    Write header dict in correct format to new bedRMod file
    :param file:
    :param header_dict:
    :return:
    """
    with open(file, 'w') as f:
        for key, value in header_dict.items():
            f.write('#' + key + '=' + value + '\n')
    return

def write_data(file, data_df):
    """
    append pandas data frame to bedRMod file
    :param file:
    :param data_df:
    :return:
    """
    column_headers = True
    with open(file, 'r') as f:
        last = f.readlines()[-1]
        if not last.startswith("#chrom"):
            column_headers = False

    with open(file, 'a') as f:
        if not column_headers:
            f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage"
                    "\tfrequency\n")
        data_df.to_csv(f, sep='\t', index=False, header=False, mode='a')

    return


def write_bedRMod(file, header_dict, data_df):
    """
    write header dict and pandas dataframe to new bedRMod file
    :param file:
    :param header_dict:
    :param data_df:
    :return:
    """
    write_header(file, header_dict)
    write_data(file, data_df)
    return