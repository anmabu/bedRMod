import pandas as pd
from Bio.Seq import Seq

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

        # check FLAGs of line and act accordingly
        skip, reverse_complemented = samtags_helper.read_flag(flag)

        if skip:  # sequences is unmapped and/or not usable in this context
            continue  # jump to next line in sam file
        if reverse_complemented:
            seq = Seq(seq).reverse_complement()
            qual = qual[::-1]  # qual is not used as of now
            # change reference letter in mod_type
            mod_type = list(map(lambda x: str(Seq(x[0]).reverse_complement() + x[1:]), mod_type))

        # there are multiple modifications possible at each index
        # mod_number indicates how many consecutive probabilites belong to each modified type
        # determine how many symbols are after the + or - (indicating the strand)
        # Check if after +/- are numbers (ChEBI)
        mod_number = []
        strand = []
        for mod in mod_type:
            if "+" in mod:
                if reverse_complemented:
                    strand.append("-")
                else:
                    strand.append("+")
                mod_code = mod.split("+")[1]
                if any(char.isdigit() for char in mod_code):
                    # if there is any digit, the assumtion is that the whole code is ChEBI
                    # I don't think something like "C+m76794" or multiple ChEBI codes concatenated are allowed
                    mod_number.append(1)
                else:
                    mod_number.append(len(mod.split("+")[1]))
            elif "-" in mod:
                if reverse_complemented:
                    strand.append("+")
                else:
                    strand.append("-")
                mod_code = mod.split("-")[1]
                if any(char.isdigit() for char in mod_code):
                    mod_number.append(1)
                else:
                    mod_number.append(len(mod.split("-")[1]))

        # change counting index, starting from 0 after finding the next modification to starting for all
        # modifications from the first base (e.g. C)
        abs_positions = samtags_helper.consecutive_to_positional_single_base(mod_index)
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
        # the following function adjusts the modification position in accordance to the complete sequence
        mod_position_df = samtags_helper.single_base_to_sequence_indices(mod_position_df, seq)
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
    sam2proEUF()
    sam_file = pd.read_csv("example_files/MM-multi.sam", comment="@", delimiter="\t", header=None)
    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = sam_file.iloc[1, 0:11]
