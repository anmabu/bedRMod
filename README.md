# bedRMod

This project is about converting RNA sequencing data into the new epitranscriptomics unified data format(EUF), **bedRMod**. 
bedRMod does contain the read data of modifications a per site-level (as opposed to SAM/BAM with Mm tags which contain the information per read-level).  

The available options are conversion from: 
- xlsx and its relatives (e.g. xls)
- odf and its relatives (e.g. ods)
- csv, tsv, basically any common seperator


# Specification
For the data specification, please refer to the bedRMod.pdf.

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
tbd

### 2.2 Using the command line
tbd