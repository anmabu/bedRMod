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


def read_bitflag(flag):
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


if __name__ == "__main__":
    print(get_SAMtags())
    print(scale_score_ML_tag(229))