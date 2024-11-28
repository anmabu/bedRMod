# bedRMod

A python project about storing RNA modification data into the new epitranscriptomics unified data format, **bedRMod**. 
This repository contains the code for conversion of pandas dataframes and plain-text files (csv, tsv, etc.) into bedRMod format.

A bedRMod file consists of two parts: A header section and a data section. 
In the header section meta-information about the data is stored. 
This includes information about the RNA itself such as source organism, assembly and annotation information as well as more technical aspects such as sequencing plattform, basecalling model, bioinformatic workflow and information about the experimental protocol. 
Further header fields can be added freely. 

Each data row of a bedRMod file contains information of one RNA modification site. 
This includes the reference segment/chromosome, positions, modification name (MODOMICS short name), modification confidence, strand, coverage and modification frequency. 
As bedRMod follows the BED9+2 convention column for visualizing the file with e.g. IGV are also contained. 
The last header row before the data section contains the column names of the bedRMod file. 


## Specification
For the data specification, please refer to the [bedRModv1.8.pdf](bedRModv1.8.pdf).

# Installation
For ease of use a python package has been developed of this repository. 
In a python(>=3.9) environment use `pip install bedRMod` to download the package from pypi. 


# Usage Information
To convert RNA sequencing data into bedRMod a few requirements have to be met. 
Those differ for the input formats. 

## 1. Config file setup
Regardless of using the GUI or the API, a config file needs to be configured. 
An example of a config.yaml file can be found [here](examples/example_config.yaml).
Only the header fields with the `# required` comment need a value but is it **highly recommended** to use the other fields as well to keep track of the data.

## 2. Converting the data
### 2.1 Using the GUI
Run `bedRMod-gui` from the command line after installing with pip.

When starting the GUI, the user has to select the input file, config file and output file, individually. 
If a config file does not exist yet, a new one can be created from a template. 
It is highly recommended for the input file to have a header aka column names in the first row, as the first row is parsed to give selectable options for the required information.
As the columns cannot be processed further in the GUI, e.g. split a column if there are several values, all more sophisticated operations have to be done on the input file beforehand.
Minor changes/adaptations of the values in the columns can still be done in the GUI, though. 
This includes selecting whether the position is 0- or 1-indexed  (counting start from 0 like birthdays or 1 like enumeration).
If the input file does not contain information on the modification type or the strand these can be set for the whole file, in the GUI.
Also functions can be passed to adapt score, coverage and frequency e.g. rounding for converting a float to an integer or scaling of the values. 

Using the GUI is recommended for converting already existing, single files into bedRMod and for users getting to know the conversion toolkit. 

### 2.2 Using the API
Once installed, the 
Creating bedRMod files can be achived using the `df2bedrmod` function. 
This takes a pandas.DataFrame as well as a config.yaml file as input and creates a bedRMod file. 
For this, the data needs to be manipulated beforehand to comply with the specifications as defined in [bedRModv1.8.pdf](bedRModv1.8.pdf).
Then, the dataframe and the config.yaml can be given to the `df2bedrmod` function in a similar way as follows:
```angular2html
import pandas as pd
from bedRMod import df2bedRMod, csv2bedRMod

# define column names
columns = ['chrom', 'start_col', 'end', 'name', 'score_column', 'strandedness', 'thick_start', 'thick_end', 'item_rgb', 'coverage_col', 'frequency_col']
data = [
    ['chr1', 1000, 1001, 'm1A', 900, '+', 1000, 1001, '0', 30, 90],
    ['2', 1001, 1002, 'm5C', 850, '-', 1001, 1002, '0', 45, 83],
    ['chrX', 3001, 3002, 'm1A', 700, '+', 3001, 3002, '0', 6, 79],
    ['chrY', 4001, 4002, 'm3C', 920, '-', 4001, 4002, '0', 50, 37],
    ['chrM', 5001, 5002, 'm3C', 920, '-', 5001, 5002, '0', 23, 65]
]

def score_func(params):
    """ calculating the score using two columns """
    score, cov = params
    return round(score / cov)

# creating the dataframe from the test data with column names
df = pd.DataFrame(data, columns=columns)

# converting pandas.DataFrame to bedRMod by passing the columnnames as keyword arguments
df2bedRMod(df, "test_config.yaml", "output_df2bedrmod.bedrmod", ref_seg="chrom", start="start_col",
           modi="name", modi_column=True, score=["score_column", "coverage_col"], score_function=score_func,
           strand="strandedness", coverage="coverage_col", frequency="frequency_col")
```

If the data already exists in a plain-text format with a delimiter (e.g. csv, tsv) and only needs minimal manipulation to meet the bedRMod specifications, it can be converted using the `csv2bedrmod` function. 
The function still requires the column names to be keyword arguments. 
Also, the first line of the input file should contain the column names. 
Pandas inferes the delimiter of the file. To avoid ambiguity it can be passed as keyword-argument to the function, as well. 
If the keyword-argument `output_file` is unused, the .bedrmod file will be written to the same directory as the input file. 
```angular2html
# continuing the code above
df.to_csv("test_csv2bedrmod.csv", index=False)

def cov_func(param):
    """ adapt values in coverage column to reflect true coverage value """
    return param * 2

def start_func(param):
    """ correcting index as bedRMod has 0-based indices """
    return param - 1

# converting csv file to bedRMod
csv2bedRMod("test_csv2bedrmod.csv", "test_config.yaml", ref_seg="chrom", start="start_col",
            start_function=start_func, modi="name", modi_column=True, score="score_column",
            strand="strandedness", coverage="coverage_col", coverage_function=cov_func,
            frequency="frequency_col")
```

# Issues
Feedback is welcome! If you encounter any bugs, have a question, or would like to request a new feature, please open an [issue](https://github.com/anmabu/bedRMod/issues) in this repository.