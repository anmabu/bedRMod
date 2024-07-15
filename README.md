# bedRMod

This project is about converting RNA sequencing data into the new epitranscriptomics unified data format(EUF), **bedRMod**. 
bedRMod does contain the read data of modifications a per site-level (as opposed to SAM/BAM with Mm tags which contain the information per read-level).  

The available options are conversion from: 
- structured formats e.g. xlsx, odf, tsv, csv
- pileup data (converted from fastq/fasta)

## Specification
For the data specification, please refer to the bedRMod.pdf.

## Instruction

To convert RNA sequencing data into bedRMod a few requirements have to be met. 
Those differ for the input formats. 

### Converting from structured formats

A config file is needed in which the metadata of the files are stored. 
Please have a look at the /test/test_config.yaml file to get a better impression. 

This description is a work in progress.

### Converting from pileup
A pileup file contains read results per site and can be directly converted from fasta/fastq files using [SAMtools](http://www.htslib.org/).
Once the data has been converted to pileup (using `fastx2pileup`), a proEUF file has to be constructed using `pileup2proEUF`.
The resulting file is the starting point to convert into EUF using `proEUF2bedRMod`. 

To not only store the per site information but also modification data, an additional file containing the modified sites is required. 
The path to this file needs to be specified in the config file and the modification data has to meet some format definitions. 
Please have a look at `mod_indices.csv` to see an example. 
For now, only information on reference segment/chromosome, position, strandedness and modification type are included. 


This function `proEUF2bedRMod` takes also 3 arguments: 
- path to proEUF (input) file
- path to config.yaml file, optional: with path to file containing modification indices. 
- path to output file 


#### Converting into pileup
Converting into pileup format can be most easily achieved by using [SAMtools](http://www.htslib.org/).
When BAM files are available, they (can be merged and) have to be sorted and indexed before their conversion into pileup format. 
It is recommended to set the options `-A -Q 0 -d 1000000 -x` when calling `samtools mpileup`. 

