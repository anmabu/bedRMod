import sys


def check_if_intron(pileup_line: list):
    """
    Function checks whether a given line from the pileup file only contains sequence skips (e.g. '<' and '>')
    :param pileup_line: Given line from the pileup file
    :return: Boolean value 'is_intron'. True if pileup_line only consists of skips, else false.
    """
    is_intron = True
    # pileup_line = pileup_line.split()
    if len(pileup_line) <= 4:  # Check if pileup is empty
        return False
    else:
        for symbol in pileup_line[4]:
            if symbol not in ['<', '>']:
                return False
    return is_intron


def parse_pileline(aligned_nts: str, ref_nt: str, curr_line: bool, count_dict: dict, count_dict_2: dict,
                   jump_dict: dict):
    """
    Counts the number of occurences for each base given the pileup data. Also, counts the number of reference skips in
    the pile, since '<' and '>' are included in the coverage and have to be removed.
    :param aligned_nts: The 5th column of the pileup-line. Bases at that position from aligned reads
    :param ref_nt: The reference nucleotide of the line. This is done to calculate the arrest rate of the
    RT in respect to the nt.
    :param curr_line: Boolean, indicating whether the current_line (True) or the following line (False) is processed
    :param count_dict: Dictionary to count instances of bases at that position in the current line
    :param count_dict_2: Dictionary to count instances of bases at that position in the next line
    :param jump_dict: Dictionary counting the jumps caused by deletion
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
                if elem == '-' and number_indels <= 3 and curr_line:  # why are the indels only counted when less than 3? -> because the jump dict only goes to 3???
                    (jump_dict['0'])[number_indels] += 1  # Count jumps caused by del.
                number_indels += 1  # add one to skip the right amount of positions when indel happened
            # $-case requires no action; $ marks the end of a read
            elif elem in ['<', '>']:
                skip_counter += 1

    # add matches for the reference nucleotide
    if curr_line:
        count_dict[ref_nt] = count_dict["."] + count_dict[","]
    else:
        count_dict_2[ref_nt] = count_dict_2["."] + count_dict_2[","]
    return skip_counter, count_dict, count_dict_2, jump_dict


def write_output(outfile, splitted_line, coverage, pre_base, rel_mism, single_jump_rate_direct,
                 single_jump_rate_delayed,
                 double_jump_rate, arrest, count_dict):
    """
    Totally not sure this is the best way to do it...
    Writes the given values to the output file which has been opened e.g in a context manager.
    Output format: ref_seg pos ref_base cov pre_base rel_mism  A  G  T  C  N  a  g  t  c  n  arrest_rate
    :param outfile: name of opened output file
    :param splitted_line: list containing the pileup line, seperated at \t
    :param coverage: number of aligned nts at this position  # TODO: This does not make sense whatsoever! Look up definition of coverage and convert!
    :param pre_base: reference_nt of the following pileup line.
    :param rel_mism: relative mismatch of nts at this position
    :param single_jump_rate_direct: jumps compared to previous postion
    :param single_jump_rate_delayed: jumps compared to previous previous position
    :param double_jump_rate: jumps compared to previous previous positon where jump already occured TODO: WTF?????
    :param arrest: arrest rate TODO: is this calculated properly?
    :param count_dict: Dictionary, containing the counted occurences of each nt
    """

    outfile.write(f"{splitted_line[0]}\t{splitted_line[1]}\t{splitted_line[2]}\t{coverage}\t{pre_base}\t{rel_mism}"
                  f"\t{count_dict['A']}\t{count_dict['G']}\t{count_dict['T']}\t{count_dict['C']}\t{count_dict['N']}"
                  f"\t{count_dict['a']}\t{count_dict['g']}\t{count_dict['t']}\t{count_dict['c']}\t{count_dict['n']}"
                  f"\t{single_jump_rate_direct}\t{single_jump_rate_delayed}\t{double_jump_rate}\t{arrest}\n")


def calculate_features(outfile, jump_dict: dict, splitted_line_prev: list, splitted_line_current: list,
                       consecutive: bool):
    """
    :param outfile:
    :param jump_dict:
    :param splitted_line_prev:
    :param splitted_line_current:
    :param consecutive: True if lines are consecutive. Arrest rate can be calculated if true, because count/coverage
    of current and next position are known.
    :return jump_dict:
    """
    char_string = 'agtcnAGTCNXRKS.,*_'
    # Initialize count-dicts with zeroes
    count_dict = {char_x: 0 for char_x in char_string}
    count_dict_2 = {char_x: 0 for char_x in char_string}

    aligned_bases_prev = splitted_line_prev[4]
    ref_nucleotide_prev = splitted_line_prev[2]
    num_aligned_reads_prev = float(splitted_line_prev[3])
    aligned_bases_current = splitted_line_current[4]
    ref_nucleotide_current = splitted_line_current[2]
    num_aligned_reads_current = float(splitted_line_current[3])

    skips, count_dict, count_dict_2, jump_dict = parse_pileline(aligned_bases_prev, ref_nucleotide_prev, True,
                                                                count_dict, count_dict_2, jump_dict)
    coverage = num_aligned_reads_prev - float(skips)  # aka sequence coverage or depth
    # but num_alinged_reads_prev is already defined as the depth at this position. Aren't skips already deducted?

    rel_mism = 0.0
    single_jump_rate_direct = 0.0
    single_jump_rate_delayed = 0.0
    double_jump_rate = 0.0
    arrest = 0.0

    if coverage > 0.0:
        # "* (asterisk) is a placeholder for a deleted base in a multiple basepair deletion
        # that was mentioned in a previous line"
        if (coverage - count_dict["*"]) != 0:
            rel_mism = 1 - (count_dict[ref_nucleotide_prev] / (coverage - count_dict["*"]))
        else:
            rel_mism = 0
        # rel_mism = (coverage - count_dict["."] - count_dict[","] - count_dict["*"]) / (coverage - count_dict["*"])
        single_jump_rate_direct = (jump_dict['-1'])[1] / coverage
        single_jump_rate_delayed = (jump_dict['-2'])[1] / coverage
        double_jump_rate = (jump_dict['-2'])[2] / coverage

    coverage -= count_dict["*"]  # * is a placeholder for a deleted base.
    # Why adjust coverage after calculating mismatch and jump rates?

    if consecutive:
        skips_next, _, _, jump_dict = parse_pileline(aligned_bases_current, ref_nucleotide_current, False,
                                                     count_dict, count_dict_2, jump_dict)
        coverage_next = num_aligned_reads_current - float(skips_next)
        if num_aligned_reads_current > 0 and coverage_next > 0:  # if the given lines are successive
            # -> calculate arrest rate.
            arrest = (aligned_bases_current.count('^') / coverage_next)

    # does it make sense to use splitted_line_prev but current ref nt?
    # maybe it's cDNA -> inverted order of sequences
    write_output(outfile, splitted_line_prev, coverage, ref_nucleotide_current, rel_mism, single_jump_rate_direct,
                 single_jump_rate_delayed, double_jump_rate, arrest, count_dict)

    return jump_dict


def write_header(file):
    """
    Writes header to a file in an open context manager.
    :param file:
    :return:
    """
    file.write('ref_seg\tpos\tref_base\tcov\tpre_base\tmism_rate\tA\tG\tT\tC\tN\ta\tg\tt\tc\tn'
               '\tsingle_jump_rate_direct\tsingle_jumprate_delayed\tdouble_jump_rate\tarrest_rate\n')


def pileup2profile(input_file, output_file):
    """
    :param input_file: (path to) input file
    :param output_file: (path to) output file
    """
    # Dictionary containing information on the jumps over multiple positions.
    # TODO: check whether this is makes sense at all!
    jump_dict = {'0': {1: 0, 2: 0, 3: 0}, '-1': {1: 0, 2: 0, 3: 0}, '-2': {1: 0, 2: 0, 3: 0}, '-3': {1: 0, 2: 0, 3: 0}}

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
            if prev_line[0] == current_line[0] and (int(prev_line[1]) == (int(current_line[1]) - 1)):
                jump_dict = calculate_features(outfile, jump_dict, prev_line, current_line, True)
            else:  # if intron has been inbetween, set consecutive bool to False
                # Does this make sense?
                jump_dict = calculate_features(outfile, jump_dict, prev_line, current_line, False)
            # Move the jump-dicts by one position for the next line to be processed
            for i in range(-3, 0):
                jump_dict[str(i)] = jump_dict[str(i + 1)]
            jump_dict['0'] = {1: 0, 2: 0, 3: 0}

            prev_line = current_line
    print(f"{input_file} converted to {output_file}!")


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    pileup2profile(infile, outfile)
