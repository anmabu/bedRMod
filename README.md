# EUF

This project is about converting RNA sequencing data into the new epitranscriptomics unified data format (EUF). 
EUF does contain the read data on a per site-level (as opposed to SAM/BAM which contain the information per read-level) 
as well as information on modification of the site. 

The available options are conversion from: 
- flat (tsv) format
- pileup data (converted from fastq/fasta)

### Instruction

To convert RNA sequencing data into EUF a few requirements have to be met. 
Those differ for the input formats. 

#### Convert from tsv (flat) format
**!!It is currently not possible to store modification data in the output files with this conversion method!!**

A config file is needed in which the metadata of the files are stored. 
Please have a look at the config.yaml file to get a better impression. 
To convert the file, call the `tsv2euf` function with the following arguments: 
- path to tsv (input) file e.g. "/flat2euf/m6aSACseq/GSE198246/GSE198246_2ng_sites.tsv.gz"
- path to config.yaml file e.g. "/tsv2euf/config.yaml"
- path to output file e.g. "/flat2euf/m6aSACseq/euf/output_file.bed"

#### Convert from pileup
A pileup file contains read results per site and can be directly converted from fasta/fastq files using [SAMtools](http://www.htslib.org/).
Once the data has been converted to pileup (workflow will be provided), a proEUF file has to be constructed using `pileup2proEUF`.
The resulting file is the starting point to convert into EUF using `proEUF2euf`. 

To not only store the per site information but also modification data, an additional file containing the modified sites is required. 
The path to this file needs to be specified in the config file and the modification data has to meet some format definitions. 
Please have a look at `mod_indices.csv` to see an example. 
For now, only information on reference segment/chromosome, position, strandedness and modification type are included. 


This function `proEUF2euf` takes also 3 arguments: 
- path to proEUF (input) file
- path to config.yaml file, optional: with path to file containing modification indices. 
- path to output file 

## Known Issues
- It is not possible to open the file "as is" in IGV. This is due to the file having too many column for IGV. If only the first 11 columns are used, IGV can open them without problems
