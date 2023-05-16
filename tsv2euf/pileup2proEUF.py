# this is modified compared to modA_ML!!

import sys


def check_if_intron(pileup_line: list):
    """
    Function checks whether a given line from the pileup file only contains sequence skips (e.g. '<' and '>')
    :param pileup_line: Given line from the pileup file
    :return: Boolean value 'is_intron'. True if pileup_line only consists of skips, else false.
    """
    is_intron = True
    if len(pileup_line) <= 4:  # Check if pileup is empty
        return False
    else:
        for symbol in pileup_line[4]:
            if symbol not in ['<', '>']:
                return False
    return is_intron


def parse_pileline(aligned_nts: str, ref_nt: str, curr_line: bool, count_dict: dict, count_dict_2: dict):
    """
    Counts the number of occurences for each base given the pileup data. Also, counts the number of reference skips in
    the pile, since '<' and '>' are included in the coverage and have to be removed.
    :param aligned_nts: The 5th column of the pileup-line. Bases at that position from aligned reads
    :param ref_nt: The reference nucleotide of the line. This is done to calculate the arrest rate of the
    RT in respect to the nt.
    :param curr_line: Boolean, indicating whether the current_line (True) or the following line (False) is processed
    :param count_dict: Dictionary to count instances of bases at that position in the current line
    :param count_dict_2: Dictionary to count instances of bases at that position in the next line
    :return:
    """
    skip_counter = 0
    quality_caret = False
    number_indels = 0
    for index, elem in enumerate(aligned_nts):
        if elem == "$":  # marks end of read segment an contains no other information
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
            if curr_line:
                count_dict[elem] += 1  # Counts all aligned bases within the line
            else:
                count_dict_2[elem] += 1
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

    # add matches for the reference nucleotide
    if curr_line:
        count_dict[ref_nt] = count_dict["."] + count_dict[","]
    else:
        count_dict_2[ref_nt] = count_dict_2["."] + count_dict_2[","]
    return skip_counter, count_dict, count_dict_2


def write_output(outfile, splitted_line, coverage, strand):
    """
    Totally not sure this is the best way to do it...
    Writes the given values to the output file which has been opened e.g in a context manager.
    Output format: ref_seg pos ref_base cov
    :param outfile: name of opened output file
    :param splitted_line: list containing the pileup line, seperated at \t
    :param coverage: number of aligned nts at this position
    """
    outfile.write(f"{splitted_line[0]}\t{splitted_line[1]}\t{splitted_line[2]}\t{coverage}\t{strand}\n")


def calculate_features(outfile, splitted_line_prev: list):
    """
    :param outfile:
    :param splitted_line_prev:
    :return:
    """
    char_string = 'agtcnAGTCNXRKS.,*_'
    # Initialize count-dicts with zeroes
    count_dict = {char_x: 0 for char_x in char_string}
    count_dict_2 = {char_x: 0 for char_x in char_string}

    aligned_bases_prev = splitted_line_prev[4]
    ref_nucleotide_prev = splitted_line_prev[2]
    num_aligned_reads_prev = float(splitted_line_prev[3])

    skips, count_dict, count_dict_2 = parse_pileline(aligned_bases_prev, ref_nucleotide_prev, True,
                                                     count_dict, count_dict_2)
    coverage = num_aligned_reads_prev - float(skips)  # aka sequence coverage or depth
    # but num_alinged_reads_prev is already defined as the depth at this position. Aren't skips already deducted?

    strand = "+" if sum(count_dict.values()) >= sum(count_dict_2.values()) else "-"

    coverage -= count_dict["*"]  # * is a placeholder for a deleted base.
    # Why adjust coverage after calculating mismatch and jump rates?

    # does it make sense to use splitted_line_prev but current ref nt?
    # maybe it's cDNA -> inverted order of sequences
    write_output(outfile, splitted_line_prev, coverage, strand)

    return


def write_header(file):
    """
    Writes header to a file in an open context manager.

    :param file:
    :return:
    """
    file.write('ref_seg\tpos\tref_base\tcov\tstrand\n')


def pileup2proEUF(input_file, output_file):
    """
    converts input file in pileup format into a profile. This is used for conversion into EUF!
    Features needed to convert into euf:
    all lines starting with # will not be worked on in this function
    chrom: chromosome number
    chromStart: position of base
    # chromEnd: position of base +1
    # name: modification name (".") if no modification
    # score: (use quality score of read later?) 0 for now
    # maybe use mismatched score at the position or something like that?
    # thickStart: = chromStart
    # thickEnd: = chromEnd
    # (itemRgb): color value for displaying modification, leave blank here
    coverage: number of reads at that position
    # frequency: percentage of modified reads at that position, read from modification file
    refBase: reference base at this position.
    :param input_file: (path to) input file
    :param output_file: (path to) output file
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        write_header(outfile)  # initialze output file
        prev_line = []
        current_line = []
        for index, line in enumerate(infile):
            if index == 0:  # initialize first line.
                prev_line = line.split("\t")
                continue
            current_line = line.split("\t")
            if check_if_intron(current_line):  # if current line is part of an intron
                # while part of an intron, don't calculate features and check in next iteration
                continue
            # current line is not part of an intron
            # check for same segment of both lines and for continuous sequence positions
            calculate_features(outfile, prev_line)
            prev_line = current_line
    print(f"{input_file} converted to {output_file}!")


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    pileup2proEUF(infile, outfile)
