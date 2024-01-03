import pandas as pd


def get_SAMtags():
    # the ambiguity codes are not from og samtags but stolen from Modomics
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


def scale_probability_ML_tag(prob):
    """
    Convert the probability value given in the ML-Tag into the range of 0 - 1
    "Thus the probability range corresponding to integer value N is N/256 to (N + 1)/256." - according to SAMtags
    :param prob:
    :return:
    """
    return prob/256


if __name__ == "__main__":
    print(get_SAMtags())
    print(scale_probability_ML_tag(229))