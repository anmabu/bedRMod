# bedRMod

This project is about converting RNA modification data into the new epitranscriptomics unified data format(EUF), **bedRMod**. 
bedRMod contains the read data of RNA modifications per site (as opposed to SAM/BAM with Mm tags which contain the information per read-level).  

The available options are conversion from: 
- xlsx and its relatives (e.g. xls)
- odf and its relatives (e.g. ods)
- csv, tsv, basically any common seperator


## Specification
For the data specification, please refer to the bedRModv*.pdf.

# Installation
In a python environment use `pip install bedRMod` to download the package from pypi. 

# Usage Information
To convert RNA sequencing data into bedRMod a few requirements have to be met. 
Those differ for the input formats. 

## 1. Config file
Independent of the input format, a config file is needed for successful conversion. 
This file contains metadata of the input file.
Please have a look at the example config file in /examples/example_config.yaml. 
For further information, please refer to the specification

## 2. Converting the data
### 2.1 Using the GUI
When starting the GUI, the user has to select the input file, config file and output file, individually. 
If a config file does not exist yet, a new one can be created from a template. 
It is highly recommended for the input file to have a header aka column names in the first row, as the first row is parsed to give selectable options for the required information.
As the columns cannot be processed further in the GUI, e.g. split a column if there are several values, all more sophisticated operations have to be done on the input file beforehand.
Minor changes/adaptations of the values in the columnscan still be done in the GUI, though. 
This includes selecting whether the position is 0- or 1-indexed  (counting start from 0 like birthdays or 1 like enumeration).
If the input file does not contain information on the modification type or the strand these can be set for the whole file, in the GUI.
Also functions can be passed to adapt score, coverage and frequency e.g. rounding for converting a float to an integer or scaling of the values. 

Using the GUI is recommended for converting single files into bedRMod and users getting to know the conversion toolkit. 

### 2.2 Using the API
tbd