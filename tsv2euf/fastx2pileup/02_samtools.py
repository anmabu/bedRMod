import os
from multiprocessing import Pool
import subprocess


# convert sam to bam
def sam2bam(file):
    if file.endswith(".sam"):
        outfile = file[:-3] + "bam"
        subprocess.call(f'samtools view -h -o {outfile} {file}', shell=True)


# merge both bam files
def merge_bams(file, second_file):
    if file.endswith(".bam") and ("R1" in file) and ("R1_R2" not in file):
        outfile = file.replace("R1", "both")
        subprocess.call(f"samtools merge -o {outfile} {file} {second_file}", shell=True)


# sort merged bam files
def sort_bam(file):
    if file.endswith(".bam") and "both" in file:
        outfile = file[:-3] + "sorted.bam"
        subprocess.call(f"samtools sort {file} -o {outfile}", shell=True)


# index sorted bam files
def index_bam(file):
    if file.endswith(".sorted.bam"):
        subprocess.call(f"samtools index {file}", shell=True)


def bam2pileup(bamfile, fastafile):
    """
    converts bamfile into a pileup file using fastafile as a reference.
    The working directory should be the directory where the bamfiles are stored.
    :param bamfile: file in bam format
    :param fastafile: file in fasta format (or fasta derivative) with the reference sequence
    """
    output_file = bamfile[:-3] + "pileup"
    subprocess.call(f"samtools mpileup -f {fastafile} -o {output_file} {bamfile}", shell=True)


def bam2pileup_MT(bamfile):
    """
    converts bam into pilupe without using a fastafile as reference.
    It was done like this in the Master's Thesis.
    :param bamfile:
    :return:
    """
    output_file = bamfile[:-3] + "pileup"
    # the -B option disables the calculation of the Base Alignment Quality which would be calulated
    # if a reference i.e. GCA or GCF was supplied
    subprocess.call(f"samtools mpileup -B -d 1000000 -Q 0 -o {output_file} {bamfile}", shell=True)
