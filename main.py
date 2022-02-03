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

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    # create lists for species sequences and headers, as well as an empty list for scores
    other_species_BRD2 = [gg_seq, mm_seq, br_seq, tt_seq]
    other_species_headers = [gg_header, mm_header, br_header, tt_header]
    scores = []
    
    # align each species to human, record scores
    Alignment=NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    for species_num in range(4):
        compare_to_human = Alignment.align(hs_seq, other_species_BRD2[species_num])
        scores.append(compare_to_human[0]) # add score to list
    
    # sort scores from largest to smallest
    arranged_scores = sorted(scores, reverse=True)
    
    # print species for each score from largest to smallest
    for score in arranged_scores:
        score_orig_index = scores.index(score)
        print(other_species_headers[score_orig_index])

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print(arranged_scores)

if __name__ == "__main__":
    main()
