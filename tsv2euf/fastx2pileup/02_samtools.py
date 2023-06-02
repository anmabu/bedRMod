import os
import multiprocessing
# from multiprocessing import Pool
import subprocess


# convert sam to bam
def sam2bam(file):
    if file.endswith(".sam"):
        outfile = file[:-3] + "bam"
        subprocess.call(f'samtools view -h -o {outfile} {file}', shell=True)
        print(f"converted {file} to {outfile}")


def sam2bam_dir(input_dir):
    os.chdir(input_dir)
    file_list = [file for file in os.listdir() if file.endswith(".sam")]
    num_processes = multiprocessing.cpu_count() - 2  # don't use all threads to not disable use of computer
    with multiprocessing.Pool(num_processes) as p:
        p.map(sam2bam, file_list)


# merge both bam files
def merge_bams(file, second_file):
    if os.path.exists(file) and os.path.exists(second_file):
        if file.endswith(".bam") and ("R1" in file) and ("R1_R2" not in file):
            outfile = file.replace("R1", "both")
            subprocess.call(f"samtools merge -f -o {outfile} {file} {second_file}", shell=True)
            print(f"merged {file} and {second_file} to {outfile}")


def merge_bams_dir(input_dir):
    os.chdir(input_dir)
    file_list = [file for file in os.listdir() if (file.endswith(".bam") and ("R1" in file) and ("both" not in file) and ("R1_R2" not in file))]
    second_file_list = [file.replace("R1", "R2") for file in file_list]
    second_file_list = [file.replace("nofw", "norc") for file in second_file_list]  # condition used with bowtie2
    paired_files = list(zip(file_list, second_file_list))
    num_processes = multiprocessing.cpu_count() - 2  # don't use all threads to not disable use of computer
    with multiprocessing.Pool(num_processes) as p:
        p.starmap(merge_bams, paired_files)


# sort merged bam files
def sort_bam(file):
    if file.endswith(".bam") and "both" in file and not file.endswith(".sorted.bam"):
        outfile = file[:-3] + "sorted.bam"
        subprocess.call(f"samtools sort {file} -o {outfile}", shell=True)
        print(f"sorted {file}")


def sort_bam_dir(input_dir):
    os.chdir(input_dir)
    file_list = [file for file in os.listdir() if (file.endswith(".bam") and ("both" in file))]
    num_processes = multiprocessing.cpu_count() - 2  # don't use all threads to not disable use of computer
    with multiprocessing.Pool(num_processes) as p:
        p.map(sort_bam, file_list)


# index sorted bam files
def index_bam(file):
    if file.endswith(".sorted.bam"):
        subprocess.call(f"samtools index {file}", shell=True)
        print(f"indexed {file}")


def index_bam_dir(input_dir):
    os.chdir(input_dir)
    file_list = [file for file in os.listdir() if file.endswith(".sorted.bam")]
    num_processes = multiprocessing.cpu_count() - 2  # don't use all threads to not disable use of computer
    with multiprocessing.Pool(num_processes) as p:
        p.map(index_bam, file_list)


def bam2pileup(bamfile, fastafile):
    """
    converts bamfile into a pileup file using fastafile as a reference.
    The working directory should be the directory where the bamfiles are stored.
    :param bamfile: file in bam format
    :param fastafile: file in fasta format (or fasta derivative) with the reference sequence
    """
    p1, sep, p2 = bamfile.partition("sorted")
    output_file = p1 + "pileup"
    # https://github.com/samtools/samtools/issues/910
    # enabling -s returns 7th column in pileup which contains the Mapping Quality in a Phred+33 Score
    # instead: using output-extra MAPQ returns value of MAPQ between 0 and 255. Maybe use as display color?
    # https://github.com/samtools/samtools/issues/1129#issuecomment-585172313
    subprocess.call(f"samtools mpileup -A -Q 0 -x -f {fastafile} -o {output_file} {bamfile}", shell=True)


def bam2pileup_dir(input_dir, ref_sequence, ref_sequence_path):
    os.chdir(input_dir)
    file_list = [file for file in os.listdir() if "both" in file and file.endswith(".sorted.bam") and (ref_sequence in file)]
    ref_file_list = [ref_sequence_path for file in file_list]
    paired_files = list(zip(file_list, ref_file_list))
    print(paired_files)
    num_processes = multiprocessing.cpu_count() - 2  # don't use all threads to not disable use of computer
    with multiprocessing.Pool(num_processes) as p:
        p.starmap(bam2pileup, paired_files)


def bam2pileup_MT(bamfile):
    """
    converts bam into pilupe without using a fastafile as reference.
    It was done like this in the Master's Thesis.
    :param bamfile:
    :return:
    """
    p1, sep, p2 = bamfile.partition("sorted")
    output_file = p1 + "pileup"
    # the -B option disables the calculation of the Base Alignment Quality which would be calulated
    # if a reference i.e. GCA or GCF was supplied
    subprocess.call(f"samtools mpileup -B -d 1000000 -Q 0 -o {output_file} {bamfile}", shell=True)
    print(f"converted {bamfile} to {output_file}")


def bam2pileup_MT_dir(input_dir):
    os.chdir(input_dir)
    file_list = [file for file in os.listdir() if "both" in file and file.endswith(".sorted.bam")]
    num_processes = multiprocessing.cpu_count() - 2  # don't use all threads to not disable use of computer
    with multiprocessing.Pool(num_processes) as p:
        p.map(bam2pileup_MT, file_list)


if __name__ == "__main__":
    in_dir = "/home/annebusch/Documents/PyCharmProjects/EUF/tsv2euf/test_files/"
    # sam2bam_dir(in_dir)
    # merge_bams_dir(in_dir)
    # sort_bam_dir(in_dir)
    # index_bam_dir(in_dir)
    # bam2pileup_MT_dir(in_dir)
    # bam2pileup_dir(in_dir, "GCF_new_ref", "/home/annebusch/anne02/m1A_Marco/S_cerv_R64/ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna")
    os.chdir(in_dir)
    bam2pileup("MH1601_both_GCF_ref_localN1L10nofwD20R3k1.sorted.bam", "/home/annebusch/anne02/m1A_Marco/S_cerv_R64/ncbi_dataset/data"
                                              "/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna")
