# create mod indices file.

import pandas as pd
from Bio import SeqIO


# reads in the file containing the modomics abbrevations and one letter codes
def get_modomics_shortvalues_df(input_file):
    """
    processes the information of Modomics' modification abbreviations and one letter codes with their "translation" as input and return dataframe to easily parse the information.
    """
    # convert list of shortvalues into Dataframe for easier conversion

    modomics_shorts = open(input_file, "r").read()

    modomics_shorts = modomics_shorts.split("\n")
    modomics_shorts = [x.split(" ", 2) for x in modomics_shorts]
    modomics_shorts = pd.DataFrame(modomics_shorts, columns=["one_letter_code", "symbol", "name"])
    modomics_shorts = modomics_shorts.iloc[:-1]  # last row contains only None

    modomics_shorts = modomics_shorts
    new_column = [str(x) for y in modomics_shorts["symbol"] for x in y if str(x).isupper()]  # convert symbols into "original" base
    # for "-, _, ., P, ], Z" are no conversions done. -> add last 3 manually (pseudouridine -> U)
    new_column.extend(['U', 'U', 'U'])

    # delete first 8 rows of mod_shortvalues. Contains U, C, A, G, T, empty, underline and unknown
    modomics_shorts = modomics_shorts[8:]
    modomics_shorts = modomics_shorts.reset_index(drop=True)

    # delete first 5 rows of new_column to get correct values
    new_column = new_column[5:]
    modomics_shorts["orig_base"] = new_column

    # adjust mods dataframe. I -> A and D -> U
    modomics_shorts["orig_base"].replace("I", "A", inplace=True)
    modomics_shorts["orig_base"].replace("D", "U", inplace=True)
    return modomics_shorts


# read in the genes from tRNA
def get_all_tRNA_genes(fasta_file):
    # read in alltRNA sequences from fasta file
    records = SeqIO.parse(fasta_file, "fasta")
    data = []
    for record in records:
        gene_id, chrom, gene_start, gene_stop, trans_start, strand, exon_start, exon_stop = record.id.split("|")
        data.append([gene_id, chrom, str(record.seq), int(gene_start), int(gene_stop), int(trans_start), int(strand), exon_start, exon_stop])
    df = pd.DataFrame(data, columns=["Gene stable ID", "Chromosome/scaffold name", "cDNA Sequence", "Gene start (bp)", "Gene end (bp)", "Transcription start", "Strand", "Exon region start (bp)", "Exon region end (bp)"])
    # substitute "Mito" with "mitochondrion" to match expression in all_tRNA
    tRNA_df = df.replace("Mito", "mitochondrion")
    bk_dict = {"I": "BK006935.2", "II": "BK006936.2", "III": "BK006937.2", "IV": "BK006938.2", "V": "BK006939.2",
              "VI": "BK006940.2", "VII": "BK006941.2", "VIII": "BK006934.2", "IX": "BK006942.2", "X": "BK006943.2",
              "XI": "BK006944.2", "XII": "BK006945.2", "XIII": "BK006946.2", "XIV": "BK006947.3", "XV": "BK006948.2",
              "XVI": "BK006949.2"}
    # convert roman numerals to BK-number
    # does not make sense to do this. Better use roman numerals throughout. To make it more compatible with new data
    # tRNA_df["Chromosome/scaffold name"] = [bk_dict[chrom] if chrom in bk_dict.keys() else chrom for chrom in tRNA_df["Chromosome/scaffold name"]]
    # sort df by exon region start for the right order
    tRNA_df = tRNA_df.sort_values(by=["Chromosome/scaffold name", "Exon region start (bp)"]).reset_index(drop=True)
    return tRNA_df


def select_modification(all_tRNA_df, mod_symbol):
    """
    helper function for clean_up_modification_df.
    select for m1A modification (") in all tRNA reference sequences
    Selects the sequences that contain the modification which is indicated by the one letter mod_symbol.
    Returns a dataframe with information on the tRNAs which have this modification.
    """
    select_mod_df = all_tRNA_df[all_tRNA_df["seq"].str.contains(f'{mod_symbol}')]
     # clean up select_mod_df dataframe
    select_mod_df = select_mod_df[select_mod_df["subtype"] != "None"]
    del select_mod_df["genebank"]
    del select_mod_df["mintbase"]
    # select_mod_df = select_mod_df.drop_duplicates(subset=["anticodon"])  # don't drop multiple occurences of same anticodon! Different genes on different chromosomes!
    select_mod_df = select_mod_df.reset_index(drop=True)
    print(select_mod_df.columns)
    return select_mod_df


def clean_up_modification_df(modification_tRNA_df, modomics_shortvalues_df, mod_symbol):
    """
    tidy up the dataframe which contains the sequences with contain the (selected) modification.
    """
    modification_tRNA_df = select_modification(modification_tRNA_df, mod_symbol)
    # substitue other modifications in mod_df to only get one modification
    for index, modi in enumerate(modomics_shortvalues_df["one_letter_code"]):
        for number, sequence in enumerate(modification_tRNA_df["seq"]):
            if modi in sequence and (modi != f'{mod_symbol}'):
                modification_tRNA_df.at[number, "seq"] = sequence.replace(modi, modomics_shortvalues_df.at[index, "orig_base"])
    mod_name = modomics_shortvalues_df.loc[modomics_shortvalues_df["one_letter_code"] == mod_symbol, "symbol"].values[0]
    modification_tRNA_df[f"{mod_name}_pos"] = [sequence.index(f'{mod_symbol}') for sequence in modification_tRNA_df["seq"]]

    # substitue one letter information in anticodon to only get one modification
    for index, modi in enumerate(modomics_shortvalues_df["one_letter_code"]):
        for number, sequence in enumerate(modification_tRNA_df["anticodon"]):
            if modi in sequence:
                modification_tRNA_df.at[number, "anticodon"] = sequence.replace(modi, modomics_shortvalues_df.at[index, "orig_base"])
    return modification_tRNA_df


def get_modified_tRNAs(modifications_df, all_tRNA_genes_df):
    data = []
    # compare start of sequences if anticodon matches to find the correct positions in the genome of modifications
    for ind, values in modifications_df.iterrows():
        values["seq"] = values["seq"].replace("U", "T")
        for k, rows in all_tRNA_genes_df.iterrows():
            if values["anticodon"] in rows["Gene stable ID"]:
                data.append(rows)
    data_f = pd.DataFrame(data).reset_index(drop=True)
    return data_f


def get_modification_indices(m1A_modification_df, tRNA_start_stop_df, mod_name="m1A"):
    """
    Creates a Dataframe containing the positions of the specified modifications on the whole genome. i.e. at which position on the chromosome the modification is.
    :param m1A_modification_df: pd.DataFrame containing the anticodons and the positions of the specified modification within the tRNA sequence
    :param tRNA_start_stop_df: pd.DataFrame containing the ID, Chromosome, Strandedness, Exon start and end
    :return: Dataframe containing the positions of modifications in the tRNA
    """
    data = []
    for i, values in m1A_modification_df.iterrows():
        for j, rows in tRNA_start_stop_df.iterrows():
            values["seq"] = values["seq"].replace("U", "T")
            if values["anticodon"] in rows["Gene stable ID"]:
                # find the position of modification in the exon
                split_exon_start = list(map(int, rows["Exon region start (bp)"].split(";")))
                if len(split_exon_start) > 1:
                    split_exon_end = list(map(int, rows["Exon region end (bp)"].split(";")))
                    length_one = min(split_exon_end) - min(split_exon_start) + 1
                    length_two = max(split_exon_end) - max(split_exon_start) + 1  # adjust for the start and end nt being in the sequence
                    if rows["Strand"] > 0:  # aka positive strand
                        new_pos = values[f"{mod_name}_pos"] - length_one
                        mod_pos = max(split_exon_start) + new_pos
                    else:  # the negative strand
                        new_pos = values[f"{mod_name}_pos"] - length_two
                        mod_pos = min(split_exon_end) - new_pos
                else:
                    if rows["Strand"] > 0:
                        mod_pos = int(rows["Exon region start (bp)"]) + values[f"{mod_name}_pos"]
                    else:
                        mod_pos = int(rows["Exon region end (bp)"]) - values[f"{mod_name}_pos"]
                data.append([values["anticodon"], rows["Gene stable ID"], rows["Chromosome/scaffold name"], rows["cDNA Sequence"], rows["Exon region start (bp)"], rows["Exon region end (bp)"], rows["Strand"], mod_pos, values[f"{mod_name}_pos"]])

    return pd.DataFrame(data, columns=["anticodon", "gene_id", "ref_seg", "cDNA Sequence", "gene_start", "gene_stop", "strand", "mod_index", f"{mod_name}_original_pos"])


# this is the "pipeline" using the helper functions above
def get_reference_indices_sequences(modomics_tRNA, modomics_mods, tRNA_genes, mod_symbol, mod_name):
    """
    :param modomics_tRNA: df containing the unprocessed tRNA sequences as obtained from modomics
    :param modomics_mods: df containing the abbreviations and single-letter abbreviations of modifications used in modomics sequences
    :param tRNA_genes: df containing infomation about tRNA genes such as start and end of gene and exon(s), strandedness, ids...
    :param mod_symbol: one-letter code of selected modification
    :return: df with modification position on gene assigned to the correct gene, chromosome and anticodon
    """
    modifications_df = clean_up_modification_df(modomics_tRNA, modomics_mods, mod_symbol)
    tRNA_genes_df = get_modified_tRNAs(modifications_df, tRNA_genes)
    mod_pos_df = get_modification_indices(modifications_df, tRNA_genes_df, mod_name)
    mod_pos_df = mod_pos_df.drop_duplicates(ignore_index=True)
    return mod_pos_df


def create_mod_indices_file(output_file="../tsv2euf/example_files/mod_indices.csv", mod_symbol='"', mod_name="m1A"):
    # both from /m1A_Marco/tRNA_S_cerv_modimics. Write somewhere that these files are.
    mods = get_modomics_shortvalues_df("../tsv2euf/example_files/modomics_shortvalues.txt")
    # read in ref sequences obtained from: https://www.genesilico.pl/modomics/api/sequences?RNAtype=tRNA&organism=saccharomyces+cerevisiae&format=json
    tRNA_s_cerv = pd.read_json("../tsv2euf/example_files/tRNA_sequences.json")
    tRNA_s_cerv = tRNA_s_cerv.transpose().reset_index(drop=True)
    tRNA_s_cerv.sort_values("anticodon")
    # ensembl search of S_cerv_tRNA sequences
    tRNA_df = get_all_tRNA_genes(
        "../tsv2euf/example_files/alltRNA_sequences_chrom_start_stop_transstart_strand_exonstart_exonend.txt")
    reference_indices = get_reference_indices_sequences(tRNA_s_cerv, mods, tRNA_df, mod_symbol, mod_name)
    mod_indices = pd.DataFrame(reference_indices[["ref_seg", 'mod_index', "strand"]])
    strand_dict = {1: "+", -1: "-"}
    mod_indices["strand"] = mod_indices["strand"].replace(strand_dict)
    mod_indices["mod_type"] = [mod_name for i in mod_indices.index]
    mod_indices.to_csv(output_file, index=False)


if __name__ == "__main__":
    create_mod_indices_file()
