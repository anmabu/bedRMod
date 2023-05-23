# this is modified compared to modA_ML!!

import sys
import pandas as pd


def parse_pileline(aligned_nts: str, ref_nt: str, count_dict_current: dict):
    """
    Counts the number of occurences for each base given the pileup data. Also, counts the number of reference skips in
    the pile, since '<' and '>' are included in the coverage and have to be removed.
    :param aligned_nts: The 5th column of the pileup-line. Bases at that position from aligned reads
    :param ref_nt: The reference nucleotide of the line. This is done to calculate the arrest rate of the
    RT in respect to the nt.
    :param count_dict_current: Dictionary to count instances of bases at that position in the current line
    :return:
    """
    skip_counter = 0
    quality_caret = False
    number_indels = 0
    for index, elem in enumerate(aligned_nts):
        if elem == "$":  # marks end of read segment and contains no other information
            continue
        if elem == "^":  # when this is encountered, the following character is a quality score and not an alignment
            quality_caret = True
            continue
        if quality_caret:
            quality_caret = False
            continue
        if number_indels != 0:  # skip iterations when indel happens to avoid counting bases multiple times
            number_indels -= 1
            continue
        if elem in 'aAcCtTgG.,*':  # i.e. is [aAcCtTgG.,*]
            count_dict_current[elem] += 1  # Counts all aligned bases within the line
        if elem in '+-<>':
            if elem in '+-':  # If indel, count jumps in jump_dict
                number_indels = int(aligned_nts[index + 1])
                for j in range(index + 2, len(aligned_nts)):
                    if aligned_nts[j] in "1234567890":
                        number_indels += int(aligned_nts[j])
                    else:
                        break
                number_indels += 1  # add one to skip the right amount of positions when indel happened
            # $-case requires no action; $ marks the end of a read
            elif elem in ['<', '>']:
                skip_counter += 1

    count_dict_current[ref_nt] = count_dict_current["."] + count_dict_current[","]
    return skip_counter, count_dict_current


def calculate_features(splitted_line: list):
    """

    :param splitted_line:
    :return:
    """
    char_string = 'agtcnAGTCNXRKS.,*_'
    # Initialize count-dict with zeroes
    count_dict_current = {char_x: 0 for char_x in char_string}

    aligned_bases = splitted_line[4]
    ref_nucleotide = splitted_line[2]
    num_aligned_reads = float(splitted_line[3])

    skips, count_dict_current = parse_pileline(aligned_bases, ref_nucleotide, count_dict_current)
    coverage = num_aligned_reads - float(skips)  # aka sequence coverage or depth
    if coverage != num_aligned_reads:
        print(num_aligned_reads)
        print(coverage)
    # coverage -= count_dict_current["*"]  # * is a placeholder for a deleted base.
    # this does not make sense! the coverage is already calulated properly!!!

    # select strandedness depending on coverage at this position
    strand = "+" if count_dict_current["."] >= count_dict_current[","] else "-"

    return splitted_line[0], splitted_line[1], ref_nucleotide, coverage, strand


def pileup2proEUF(input_file, output_file):
    """
    converts input file in pileup format into a proEUF.
    The output file contains the following columns:
    reference_segment/chromosome  position  reference_base  coverage  strand
    Features needed to convert into euf:
    all lines starting with # will not be worked on in this function
    chrom: chromosome number
    chromStart: position of base
    # chromEnd: position of base +1
    # name: modification name (".") if no modification
    # score: modification confidence at this position
    strand: strand orientation at that position
    # thickStart: = chromStart
    # thickEnd: = chromEnd
    # (itemRgb): color value for displaying modification, leave blank here
    coverage: number of reads at that position
    # frequency: percentage of modified reads at that position, read from modification file
    refBase: reference base at this position.
    :param input_file: (path to) input file
    :param output_file: (path to) output file
    """
    features_list = []
    with open(input_file, "r") as infile:
        for index, line in enumerate(infile):
            prev_line = line.split("\t")
            features = calculate_features(prev_line)
            features_list.append(features)
    df = pd.DataFrame(features_list, columns=["ref_seg", "pos", "ref_base", "cov", "strand"])
    df.to_csv(output_file, index=False, sep="\t")
    print(f"{input_file} converted to {output_file}!")


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    pileup2proEUF(infile, outfile)
