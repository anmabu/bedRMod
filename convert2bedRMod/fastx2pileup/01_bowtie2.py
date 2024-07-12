import multiprocessing
from multiprocessing import Pool

import os
import subprocess


def call_bowtie2(ref_sequence, fastq_sequence_path, sam_file, bowtie_args=None, threads=3, paired=False,
                 bowtie2_dir="BT2_HOME"):
    """
    calls bowtie2 for read alignement of sequence to reference.
    step into directory containing the ref_sequence or bowtie2 will not work!!
    :param ref_sequence: name of the reference sequence file.
    :param fastq_sequence_path: path to fastq/fasta file that should be aligned to ref_sequence. The "working sequence"
    :param sam_file: path or name of resulting sam file. Should end with ".sam"
    :param bowtie_args: list or None. Arguments passed to bowtie2 which are the parameters for the alignment
    :param threads: number of threads used in alignment. Important to control if multiple files are aligned simultaneously.
    :param paired: boolean that indicates whether the read was paired (forward and reverse read exist) or whether it should be handled as a single read
    :param bowtie2_dir: path to BT2_HOME environment variable as specified in the bowtie2 documentation
    This function aligns the fastq/fasta file with the reference sequence and outputs a SAM file of the alignment as well as a plain text file containing metadata on the process.
    When setting the arguments, keep in mind that the speed does not scale linearly with the number of threads, but utilizing more threads still speeds up the process considerably.
    Having a forward and reverse read also helps to get a better alignemnt. When setting "paired" to True this function expects both reads in the same directory
    and only to differ in the read number "R1" or "R2", respectively.
    """

    bowtie_directory = os.environ.get(f"{bowtie2_dir}")
    if not bowtie_directory.endswith("/"):
        bowtie_directory += "/"

    if bowtie_args is None:
        bowtie_args = ""
    elif len(bowtie_args) > 0:
        bowtie_args = " ".join(bowtie_args)

    try:
        if os.stat(fastq_sequence_path).st_size > 0:  # non-empty file exists
            print(os.path.basename(fastq_sequence_path))
        else:
            print(f"The file at {fastq_sequence_path} is empty!")
    except FileNotFoundError:
        print(f"The file {fastq_sequence_path} does not exists!")

    if paired:
        paired_fastq_sequence = fastq_sequence_path.replace("R1", "R2")
        output = subprocess.run(f"{bowtie_directory}bowtie2 -p {threads} {bowtie_args} -x {ref_sequence} -1 "
                                f"{fastq_sequence_path} -2 {paired_fastq_sequence} -S {sam_file}",
                                shell=True, capture_output=True, text=True)
    else:
        output = subprocess.run(f"{bowtie_directory}bowtie2 -p {threads} {bowtie_args} -x {ref_sequence} -U "
                                f"{fastq_sequence_path} -S {sam_file}", shell=True, capture_output=True, text=True)
    with open(sam_file[:-3] + "txt", "w") as meta_file:
        meta_file.write(str(output))

    print(f"Converted {fastq_sequence_path} to {sam_file}!")


def mp_bowtie2(ref_sequence_dir, ref_sequence, fastq_sequence_path, output_dir, bowtie_args=None, threads=3,
               paired=False, bowtie2_dir="BT2_HOME"):
    """multiprocessing bowtie2.
    :param ref_sequence_dir: path to directory of the indexed reference sequences
    :param ref_sequence: indexed reference sequence
    :param fastq_sequence_path: path to fastq sequences' input directory.
    :param output_dir: path to output directory
    :param bowtie_args: (list of strings) contains the arguments used by bowtie. Each argument is a passed in its own
    string.
    :param threads: (int) number of threads
    :param paired: (bool) paired fastq reads
    :param bowtie2_dir: directory where bowtie2 has been installed. If installing of bowtie2 was done according to its
    user manual and bowtie2 has been added to the PATH, this is "BT2_HOME". Otherwise, the path to the directory
    containing the bowtie2 executable has to be passed as an argument.
    This function changes the working directory to the ref_sequence_directory. So make sure that the other (file)paths
    are absolut or an error occurs.
    Adjustments to (output) filenames and bowtie2 parameters also take place in this function.
    When everything is set up properly the function executing bowtie2 gets called.
    """
    os.chdir(ref_sequence_dir)
    fastq_filename = os.path.basename(fastq_sequence_path)

    run_number = "undefined"
    if paired and ("R1" in fastq_filename):
        run_number = "R1_R2"
    else:
        if "R1" in fastq_filename:
            run_number = "R1"
        if "R2" in fastq_filename:
            run_number = "R2"
            if "--nofw" in bowtie_args:
                replace_index = bowtie_args.index("--nofw")
                bowtie_args[replace_index] = "--norc"

    # print(bowtie_args)
    if bowtie_args is None:
        bowtie_args_filename = "no_args"
    else:
        bowtie_args_filename = " ".join(bowtie_args)
        bowtie_args_filename = bowtie_args_filename.replace(" ", "")
        bowtie_args_filename = bowtie_args_filename.replace("-", "")

    # this condition is redundant when only files to be converted are in the input directory
    if fastq_filename.endswith(".fastq"):
        sam_file = output_dir + "/" + fastq_filename[:6] \
                   + f"_{run_number}_{ref_sequence}_{bowtie_args_filename}" + ".sam"

        if paired:
            call_bowtie2(ref_sequence, fastq_sequence_path, sam_file, bowtie_args, threads=threads, paired=True,
                         bowtie2_dir=bowtie2_dir)
        else:
            call_bowtie2(ref_sequence, fastq_sequence_path, sam_file, bowtie_args, threads=threads, paired=False,
                         bowtie2_dir=bowtie2_dir)


def run_mp_bowtie2_dir(ref_dir, ref_sequence, fastq_dir, output_dir, bowtie_args=None, threads=2, paired=False,
                       bowtie2_dir="BT2_HOME"):
    """
    run bowtie 2 of the fastq files in the fastq_dir on indexed genome in the ref dir
    :param ref_dir: path to directory of the indexed reference sequences
    :param ref_sequence: path to indexed reference sequence
    :param fastq_dir: directory where the fastq input files are located
    :param output_dir: directory where the .sam files of the output are stored
    :param bowtie_args: (list of strings) containing the arguments passed to bowtie2
    :param threads: (int) number of threads used in each subprocess
    :param paired: (bool) indicating whether the passed sequences are paired or unpaired
    :param bowtie2_dir: directory where bowtie2 has been installed. If installing of bowtie2 was done according to its
    user manual and bowtie2 has been added to the PATH, this is "BT2_HOME". Otherwise, the path to the directory
    containing the bowtie2 executable has to be passed as an argument.
    """
    potential_files = [file for file in os.listdir(fastq_dir) if file.endswith(".fastq")]
    if len(potential_files) == 0:
        raise FileNotFoundError(f"No .fastq files in {fastq_dir}. Please check if fastq_dir is correct.")

    if bowtie_args is None:
        bowtie_args = []
    argument_list = [[ref_dir,
                      ref_sequence,
                      f"{fastq_dir}/{x}",
                      output_dir,
                      bowtie_args,
                      threads,
                      paired,
                      bowtie2_dir] for x in os.listdir(fastq_dir)]
    pool_size = int(multiprocessing.cpu_count() / threads)
    with Pool(pool_size) as p:
        p.starmap(mp_bowtie2, argument_list)


if __name__ == "__main__":
    fastq_input_dir = os.path.abspath("../other_input")
    test_output_dir = os.path.abspath("../test_files")
    ref_sequence_path = os.path.abspath("../../../../../anne02/m1A_Marco/S_cerv_R64")
    bowtie_dir = "BT2_HOME"
    run_mp_bowtie2_dir(ref_sequence_path, "GCF_ref", fastq_input_dir, test_output_dir,
                       bowtie_args=["-q", "--local", "-N 1", "-L 10", "--nofw", "-D 20", "-R 3", "-k 1"],
                       threads=7, paired=False, bowtie2_dir=bowtie_dir)
