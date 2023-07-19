import pandas as pd
from modSAM2bedRMod import samtags_helper

sam_file = pd.read_csv("example_files/MM-multi.sam", comment="@", delimiter="\t", header=None)
# don't use comment=@, because this cuts off quality strings, containing @ just skip lines starting with @

print(sam_file.iloc[:, 11:])

# iterate this through all lines!
# find column that starts with Ml and Mm
mm_tag = sam_file.iloc[0, 11]
ml_tag = sam_file.iloc[0, 12]

nt_seq = sam_file.iloc[0, 9]

if not mm_tag.startswith("Mm:"):
    print("No valid Mm tag")

if not ml_tag.startswith("Ml:"):
    print("No valid Ml tag")

ml, probs = ml_tag.split(",", 1)
probabilites = probs.split(",")
probabilites = [samtags_helper.scale_probability_ML_tag(int(x)) for x in probabilites]
print(probabilites)

mm, dt, mods = mm_tag.split(":")
mod_list = mods.split(";")
print(mod_list)

mod_type = [x.split(",", 1)[0] for x in mod_list if len(x) != 0]
mod_occur = [x.split(",", 1)[1] for x in mod_list if len(x) != 0]

print(mod_type)
print(mod_occur)


# change weird index to be consecutive
abs_positions = []
for occur in mod_occur:
    rel_positions = occur.split(",")
    abs_positions_local = []
    i = 0
    while i < len(rel_positions):
        if i == 0:
            abs_positions_local.append(int(rel_positions[i]))
            i += 1
        else:
            new_pos = abs_positions_local[i-1] + int(rel_positions[i]) + 1
            abs_positions_local.append(new_pos)
            i += 1
    abs_positions.append(abs_positions_local)

print(abs_positions)


# zip probabilities with respective position
mod_position_df = pd.DataFrame(columns=["ref_base", "strand", "mod_type", "mod_index", "probability"])
i = 0
while i < len(probabilites):
    for index, pos in enumerate(abs_positions):
        tmp_occur = []
        for j in pos:
            ref_base, strand, modification = mod_type[index][0], mod_type[index][1], mod_type[index][2:]
            new_row = [ref_base, strand, modification, j, probabilites[i]]
            mod_position_df.loc[len(mod_position_df)] = new_row
            i += 1


# count length of sequences as well as occurrences of bases and assign correct 0-based positions to modifications
mod_position_df = mod_position_df.sort_values(by="mod_index")
print(mod_position_df)

i = 0
A_counter = 0
C_counter = 0
G_counter = 0
T_counter = 0
# N_counter = 0  == i
j = 0
n_rows = len(mod_position_df)
print(nt_seq)
while i < len(nt_seq):
    curr_nt = nt_seq[i]
    curr_counter = ""
    if curr_nt == "A":
        A_counter += 1
        curr_counter = A_counter
    elif curr_nt == "C":
        C_counter += 1
        curr_counter = C_counter
    elif curr_nt == "G":
        G_counter += 1
        curr_counter = G_counter
    elif curr_nt == "T":
        T_counter += 1
        curr_counter = T_counter
    else:
        curr_counter = i

    ref_seg, strand, mod_type, pos, prob = mod_position_df.iloc[j]
    if curr_nt == ref_seg and curr_counter == pos + 1:
        pos = i
        mod_position_df.iloc[j] = [ref_seg, strand, mod_type, pos, prob]
        j += 1
    i += 1


print(mod_position_df)

# translate code to Abbreviation of modification
SAMtags = samtags_helper.get_SAMtags()

for index, row in mod_position_df.iterrows():
    base, modification = row["ref_base"], row["mod_type"]
    mod_abbrev = SAMtags.loc[(SAMtags['Unmodified base'] == base)
                             & (SAMtags['Code'] == modification), 'Abbreviation'].iloc[0]
    mod_position_df.at[index, "mod_type"] = mod_abbrev

print(mod_position_df)

# the mod index in relativ in relation to the read segment
# to make it absolute, add 1-based starting index minus one to mod_index and add column of ref_segment to df

starting_position_1_indexed = sam_file.iloc[0, 3]
ref_segment = sam_file.iloc[0, 2]

print(starting_position_1_indexed)
print(ref_segment)
# the example files don't have either...
