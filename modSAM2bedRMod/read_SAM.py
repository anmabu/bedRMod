import pandas as pd
from modSAM2bedRMod import samtags_helper

sam_file = pd.read_csv("example_files/MM-multi.sam", comment="@", delimiter="\t", header=None)
# don't use comment=@, because this cuts off quality strings, containing @ just skip lines starting with @
# the base modifications might also differ wrt the sequences stored in the SEQ field of the SAM file if
# the sequence has been reverse complemented, in this case, the 0x10 flag is set
# (in the second column of the alignment line)
# but that means that the modified positions are recorded in their original position and count the original bases

# incorporate all test cases/files from here to enable conversion
# https://github.com/samtools/htslib/tree/61b037bb881e85259f8df30c78d99ad3a357ed52/test/base_mods

# THE COORDINATES IN SAM ARE 1-BASED!


def sam2proEUF():
    sam_file = pd.read_csv("example_files/MM-multi.sam", comment="@", delimiter="\t", header=None)
    for index, row in sam_file.iterrows():
        # first 11 columns in SAM file are fixed
        qname, flag, rname, abs_pos, mapq, cigar, rnext, pnext, tlen, seq, qual = row[:11]

        # find mm and ml tags dynamically
        mm_column = []
        ml_column = []
        for column in row:
            if type(column) is str and column.startswith("Mm"):
                mm_column = column
            elif type(column) is str and column.startswith("Ml"):
                ml_column = column

        if not mm_column:
            print(f"no Mm tag found in row {index}")
        if not ml_column:
            print(f"no Ml tag found in row {index}")

        ml, probs = ml_column.split(",", 1)
        probabilites = probs.split(",")
        # convert probabilites of sam into the score
        score = [samtags_helper.scale_score_ML_tag(int(x)) for x in probabilites]

        mm, dt, mods = mm_column.split(":")
        mod_list = mods.split(";")

        # which modifications occur
        mod_type = [x.split(",", 1)[0] for x in mod_list if len(x) != 0]
        # at which relative indices these modifications happen
        # The hard and soft clips are missing!!!
        mod_index = [x.split(",", 1)[1] for x in mod_list if len(x) != 0]
        print(mod_index)

        # there are multiple modifications possible at each index
        # mod_number indicates how many consecutive probabilites belong to each modified type
        # determine how many symbols are after the + or - (indicating the strand)
        mod_number = []
        strand = []
        for mod in mod_type:
            if "+" in mod:
                strand.append("+")
                mod_number.append(len(mod.split("+")[1]))
            elif "-" in mod:
                strand.append("-")
                mod_number.append(len(mod.split("-")[1]))

        # check FLAGs of line and act accordingly
        bit_flag = bin(flag)

        # if flag == 0, do everything as before
        if flag > 0:
            i = 1
            while i < len(bit_flag) - 1:
                if i == 1 and bit_flag[-i] == "1":  # relevance
                    print("template having multiple sequences in sequencing")
                if i == 2 and bit_flag[-i] == "1":  # okay? i guess?
                    print("each segment properly aligned according to the aligner")
                if i == 3 and bit_flag[-i] == "1":
                    print("segment unmapped")
                    # bad
                if i == 4 and bit_flag[-i] == "1":  # doesn't matter
                    print("next segment in the sequence is unmapped")
                if i == 5 and bit_flag[-i] == "1":
                    print("Seq being reverse complemented")
                    # take action here
                if i == 6 and bit_flag[-i] == "1":  # doesn't matter
                    print("Seq of next segment in template being reverse complemented")
                if i == 7 and bit_flag[-i] == "1":  # doesn't matter
                    print("first segment in the template")
                if i == 8 and bit_flag[-i] == "1":  # doesn't matter
                    print("last segment in the template")
                if i == 9 and bit_flag[-i] == "1":  # ?
                    print("secondary alignment")
                if i == 10 and bit_flag[-i] == "1":  # bad!
                    print("not passing filters such as platform/vendor quality controls")
                if i == 11 and bit_flag[-i] == "1":  # okay?
                    print("PCR or optical duplicate")
                if i == 12 and bit_flag[-i] == "1":  # not a bad thing
                    print("supplementary alignment")
                i += 1

        # change counting relative index to be consecutive
        abs_positions = []
        for occur in mod_index:
            rel_positions = occur.split(",")
            abs_positions_local = []
            i = 0
            while i < len(rel_positions):
                if i == 0:
                    abs_positions_local.append(int(rel_positions[i]) + 1)
                    i += 1
                else:
                    new_pos = abs_positions_local[i - 1] + int(rel_positions[i]) + 1
                    abs_positions_local.append(new_pos)
                    i += 1
            abs_positions.append(abs_positions_local)
        print(abs_positions)

        mod_position_df = pd.DataFrame(columns=["ref_base", "strand", "mod_type", "mod_index", "score"])
        i = 0
        # iterate over length of score as positions are doubly nested lists
        while i < len(score):
            for idx, abs_pos in enumerate(abs_positions):
                for j in abs_pos:
                    if mod_number[idx] > 1:
                        for k in range(mod_number[idx]):
                            ref_base, modification = mod_type[idx][0], mod_type[idx][2 + k]
                            new_row = [ref_base, strand[idx], modification, j, score[i]]
                            mod_position_df.loc[i] = new_row
                            i += 1
                    else:
                        ref_base, modification = mod_type[idx][0], mod_type[idx][2:]
                        new_row = [ref_base, strand[idx], modification, j, score[i]]
                        mod_position_df.loc[i] = new_row
                        i += 1

        mod_position_df = mod_position_df.sort_values(by="mod_index")
        print(mod_position_df)

        # the current mod_indices reflect e.g. the number of Cs until this modified C is reached
        # the following code adjusts the modification position in accordance to the complete read
        # count length of sequences as well as occurrences of bases and assign correct 0-based positions to modifications
        for j in range(len(mod_position_df)):
            ref_base, strand, mod_type, abs_pos, score = mod_position_df.iloc[j]
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
                    mod_position_df.loc[j, "mod_index"] = i
                elif ref_base == "N" and i == abs_pos:  # nothing changes?
                    continue
                    # pos = i

        print(mod_position_df)

        # translate code to Abbreviation of modification
        SAMtags = samtags_helper.get_SAMtags()

        for ix, r in mod_position_df.iterrows():
            base, modification = r["ref_base"], r["mod_type"]
            mod_abbrev = SAMtags.loc[(SAMtags['Unmodified base'] == base)
                                     & (SAMtags['Code'] == modification), 'Abbreviation'].iloc[0]
            mod_position_df.at[index, "mod_type"] = mod_abbrev

        print(mod_position_df)

        # the mod index in relativ in relation to the read segment
        # to make it absolute, add 1-based starting index minus one to mod_index and add column of ref_segment to df
        print(abs_pos)
        print(rname)
        # the example files don't have either...


if __name__ == "__main__":
    # sam2proEUF()
    sam_file = pd.read_csv("example_files/MM-orient.sam", comment="@", delimiter="\t", header=None)
    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = sam_file.iloc[1, 0:11]
    print(flag)
    bit_flag = bin(24)
    i = 1
