from copy import deepcopy
import pandas as pd


def get_SAMtags():
    # the ambiguity codes are not from og samtags but from Modomics
    SAMtags_df = pd.DataFrame([["C", "m", "5mC", "5-Methylcytosine", "27551"],
                               ["C", "h", "5hmC", "5-Hydroxymethylcytosine", "76792"],
                               ["C", "f", "5fC", "5-Formylcytosine", "76794"],
                               ["C", "c", "5caC", "5-Carboxylcytosine", "76793"],
                               ["C", "C", "xC", "Ambiguity code; any C mod", ""],
                               ["T", "g", "5hmU", "5-Hydroxymethyluracil", "16964"],
                               ["T", "e", "5fU", "5-Formyluracil", "80961"],
                               ["T", "b", "5caU", "5-Carboxyluracil", "17477"],
                               ["T", "T", "xT", "Ambiguity code; any T mod", ""],
                               ["U", "U", "xU", "Ambiguity code; any U mod", ""],
                               ["A", "a", "6mA", "6-Methyladenine", "28871"],
                               ["A", "A", "xA", "Ambiguity code; any A mod", ""],
                               ["G", "o", "8oxoG", "8-Oxoguanine", "44605"],
                               ["G", "G", "xG", "Ambiguity code; any G mod", ""],
                               ["N", "n", "Xao", "Xanthosine", "18107"],
                               ["N", "N", "xX", "Ambiguity code; any mod", ""],
                               ], columns=["Unmodified base", "Code", "Abbreviation", "Name", "ChEBI"])
    return SAMtags_df


def scale_score_ML_tag(prob):
    """
    Convert the probability value given in the ML-Tag into the range of 0 - 1000 to fit the score of bedRMod
    "Thus the probability range corresponding to integer value N is N/256 to (N + 1)/256." - according to SAMtags
    :param prob:
    :return:
    """
    return round(prob/256 * 1000)


def read_flag(flag):
    """
    Reads flag of the line and returns whether the line should be skipped or not and signalises whether the seq is
    reverse complemented.
    Skipping the line is set to true if mapping mistakes were recorded, multiple mappings are in the same file or
    (non-)optical duplicates are recorded.
    Especially the (non-)duplicate part has to be checked, as this ignores all reads marked as duplicates even though
    at least one should be read!
    :param flag:
    :return: Bool(skip), Bool(reverse complement)
    """
    skip = False
    rc = False
    # if flag == 0, do everything as before
    if flag > 0:
        bit_flag = bin(flag)
        i = 1
        while i < len(bit_flag) - 1:
            if i == 1 and bit_flag[-i] == "1":  # relevance
                print("template having multiple sequences in sequencing")
            if i == 2 and bit_flag[-i] == "1":  # okay? I guess?
                print("each segment properly aligned according to the aligner")
            if i == 3 and bit_flag[-i] == "1":
                print("segment unmapped")
                skip = True
            if i == 4 and bit_flag[-i] == "1":  # doesn't matter
                print("next segment in the sequence is unmapped")
            if i == 5 and bit_flag[-i] == "1":
                print("Seq being reverse complemented, qual being reversed")
                rc = True
            if i == 6 and bit_flag[-i] == "1":  # doesn't matter
                print("Seq of next segment in template being reverse complemented")
            if i == 7 and bit_flag[-i] == "1":  # doesn't matter
                print("first segment in the template")
            if i == 8 and bit_flag[-i] == "1":  # doesn't matter
                print("last segment in the template")
            if i == 9 and bit_flag[-i] == "1":  # line not to be used when multiple alignments are in same file
                print("secondary alignment")
                skip = True
            if i == 10 and bit_flag[-i] == "1":  # bad!
                print("not passing filters such as platform/vendor quality controls")
                skip = True
            if i == 11 and bit_flag[-i] == "1":  # not okay!
                print("PCR or optical duplicate")
                skip = True
            if i == 12 and bit_flag[-i] == "1":  # part of chimeric alignment. not bad
                print("supplementary alignment")
            i += 1

    return skip, rc


def single_base_consecutive_to_absolute_sequence_index(mod_type, mod_index, seq):
    """
    This function converts the "consecutive positions" ((C+m, 0, 2, 1), meaning the first, fourth and sixth positions
    are modified), into the absolute 0-based positions of the modifications on the read sequence.
    :param mod_type: list containing the modification types e.g. ['C+m', 'N+n']
    :param mod_index: list of lists containing the consecutive positions of the modifications
    :param seq: read sequence where the modifications were read on
    :return: list of list containing the absolute 0-based positions of the modifications on the sequence
    """
    # print(mod_index)
    mod_index = [[int(y) for y in x.split(",")] for x in mod_index]
    new_mod_index = deepcopy(mod_index)
    for idx, indices in enumerate(mod_index):
        # the first char in mod_type is the reference base in seq
        base_type = mod_type[idx][0]
        for ind, single_index in enumerate(indices):
            curr_counter = 0
            A_counter = 0
            C_counter = 0
            G_counter = 0
            T_counter = 0
            # do I also have to add DNA/RNA ambiguity codes to the counters?
            # only canonical bases are counted when setting the index for mutated canonical bases.
            for i, base in enumerate(seq):
                if base.upper() == "A":
                    A_counter += 1
                    curr_counter = A_counter
                elif base.upper() == "C":
                    C_counter += 1
                    curr_counter = C_counter
                elif base.upper() == "G":
                    G_counter += 1
                    curr_counter = G_counter
                elif base.upper() == "T":
                    T_counter += 1
                    curr_counter = T_counter
                else:
                    curr_counter = i

                if base_type.upper() == "N":
                    curr_counter = i
                    if ind == 0:
                        if curr_counter == single_index:
                            new_mod_index[idx][ind] = i
                            # print(i)
                            # print(seq)
                            # print(i * " " + "^")
                    if ind > 0:
                        if curr_counter == single_index + indices[ind-1] + 1:
                            new_mod_index[idx][ind] = i
                            indices[ind] = single_index + indices[ind-1] + 1
                            # print(i)
                            # print(seq)
                            # print(i * " " + "^")
                elif base.upper() == base_type.upper():
                    if ind == 0:
                        if curr_counter == single_index + 1:
                            new_mod_index[idx][ind] = i
                            # print(i)
                            # print(seq)
                            # print(i * " " + "^")
                    if ind > 0:
                        if curr_counter == single_index + indices[ind-1] + 2:
                            # print(i)
                            # print(seq)
                            # print(i*" " + "^")
                            new_mod_index[idx][ind] = i
                            indices[ind] = single_index + indices[ind - 1] + 1
    # print(new_mod_index)
    return new_mod_index


def consecutive_to_positional_single_base(mod_index):
    """
    This function adjust these "consecutive positions" (C+m, 0, 2, 1) into "counting positions" (C+m, 0, 3, 5) on the
    read sequence.
    :param mod_index:
    :return:
    """
    counting_positions = []
    for occur in mod_index:
        consecutive_positions = occur.split(",")
        counting_positions_local = []
        i = 0
        while i < len(consecutive_positions):
            if i == 0:
                counting_positions_local.append(int(consecutive_positions[i]) + 1)
                i += 1
            else:
                new_pos = counting_positions_local[i - 1] + int(consecutive_positions[i]) + 1
                counting_positions_local.append(new_pos)
                i += 1
        counting_positions.append(counting_positions_local)
    return counting_positions


def single_base_to_sequence_indices(mod_df, seq):
    """
    The function converts the modification position in accordance to the complete read sequence.
    The current mod_indices reflect the number of e.g. Cs until the modified C is reached (0-based counting).
    This function changes the indices counting only one base (e.g. C) into the positional index along the sequence,
    counting all bases along the sequence.
    I'd like to further extend this function to return the positions on the reference segment of the read, but the now
    available test files don't contain neither reference segments nor read starting positions.
    :param mod_df: Dataframe with "wrong" (consecutive) mod_index information of one SAM line.
    :param seq: Read sequence of this SAM line.
    :return: Dataframe with correct (positional) mod_index information of the SAM line.
    """
    # count length of sequences as well as occurrences of bases and assign correct 0-based positions to modifications
    for j in range(len(mod_df)):
        ref_base, strand, mod_type, abs_pos, score = mod_df.iloc[j]
        print(ref_base, strand, abs_pos)
        curr_counter = 0
        A_counter = 0
        C_counter = 0
        G_counter = 0
        T_counter = 0
        for i, base in enumerate(seq):
            if base == "A" or base == "a":
                A_counter += 1
                curr_counter = A_counter
            elif base == "C" or base == "c":
                C_counter += 1
                curr_counter = C_counter
            elif base == "G" or base == "g":
                G_counter += 1
                curr_counter = G_counter
            elif base == "T" or base == "t":
                T_counter += 1
                curr_counter = T_counter
            else:
                curr_counter = i

            if base == ref_base and curr_counter == abs_pos:
                mod_df.loc[j, "mod_index"] = i
            elif ref_base == "N" and i == abs_pos:  # nothing changes?
                continue
                # pos = i
    return mod_df


if __name__ == "__main__":
    # print(get_SAMtags())
    # print(scale_score_ML_tag(229))
    # ['C+m', 'C+h', 'N+n'] mod_type
    print(single_base_consecutive_to_absolute_sequence_index(['C+m', 'C+h', 'N+n'],
                                                             ['2,2,1,4,1', '6,7', '15,2'], "AGCTCTCCAGAGTCGNACGCCATYCGCGCGCCACCA"))
