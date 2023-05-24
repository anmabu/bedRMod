# this is modified compared to modA_ML!!

import sys


def parse_pileup_line(aligned_nts: str, ref_nt: str, count_dict: dict):
    """
    Counts the number of occurences for each base given the pileup data. Also, counts the number of reference skips in
    the pile, since '<' and '>' are included in the coverage and have to be removed.
    :param aligned_nts: The 5th column of the pileup-line. Bases at that position from aligned reads
    :param ref_nt: The reference nucleotide of the line. This is done to calculate the arrest rate of the
    RT in respect to the nt.
    :param count_dict: Dictionary to count instances of bases at that position in the current line
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
            count_dict[elem] += 1  # Counts all aligned bases within the line
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

    count_dict[ref_nt] = count_dict["."] + count_dict[","]
    return skip_counter, count_dict


def get_proEUF_features(line):
    """
    this function gets the features from a pileup line for the proEUF format. 
    :param line: line from a pileup file.
    :return: ref_segment, position, ref_base, coverage, strand
    """
    splitted_line = line.split("\t")
    # don't split directly into vars because pileup contains up to 6 columns
    ref_segment = splitted_line[0]
    position = splitted_line[1]
    ref_nucleotide = splitted_line[2]
    num_aligned_reads = float(splitted_line[3])
    aligned_bases = splitted_line[4]

    char_string = 'agtcnAGTCNXRKS.,*_'
    # Initialize count-dict with zeroes
    count_dict = {char_x: 0 for char_x in char_string}

    skips, count_dict = parse_pileup_line(aligned_bases, ref_nucleotide, count_dict)
    coverage = num_aligned_reads - float(skips)  # aka sequence coverage or depth
    if coverage != num_aligned_reads:
        print(num_aligned_reads)
        print(coverage)

    # select strandedness depending on coverage at this position
    strand = "+" if count_dict["."] >= count_dict[","] else "-"

    return ref_segment, position, ref_nucleotide, coverage, strand


def pileup2proEUF(input_file, output_file):
    """
    converts input file in pileup format into a proEUF.
    The output file contains the following columns:
    reference_segment/chromosome  position  reference_base  coverage  strand
    :param input_file: (path to) input file
    :param output_file: (path to) output file
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        outfile.write("ref_seg\tpos\tref_base\tcov\tstrand\n")
        for index, line in enumerate(infile):
            ref_seg, pos, ref_base, cov, strand = get_proEUF_features(line)
            outfile.write(f"{ref_seg}\t{pos}\t{ref_base}\t{cov}\t{strand}\n")
    print(f"{input_file} converted to {output_file}!")


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    pileup2proEUF(infile, outfile)
