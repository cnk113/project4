# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    input = [gg_seq, mm_seq, br_seq, tt_seq]
    species = ['Chicken', 'Mouse', 'Stork', "Dolphin"]
    needle = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    align = []
    for i in range(len(species)):
        x,y,z = needle.align(hs_seq, input[i])
        align.append((x,species[i]))
    align.sort(key = lambda x: x[0], reverse=True)
    for x,y in align:
        print(str(x) + " " + y)


if __name__ == "__main__":
    main()
