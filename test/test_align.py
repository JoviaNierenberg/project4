# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    Alignment=NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    Alignment.align(seq1, seq2)

    align_matrix_by_hand = np.array([[  0, -np.inf, -np.inf, -np.inf],
                                     [-np.inf,   5, -11, -13,],
                                     [-np.inf, -12,   4,  -8],
                                     [-np.inf, -12,  -1,   5],
                                     [-np.inf, -14,  -6,   4]])

    gapA_matrix_by_hand = np.array([[-10, -11, -12, -13],
                                     [-np.inf, -22,  -6,  -7],
                                     [-np.inf, -23, -17,  -7],
                                     [-np.inf, -24, -18, -12],
                                     [-np.inf, -25, -19, -17]])

    gapB_matrix_by_hand = np.array([[-10, -np.inf, -np.inf, -np.inf],
                                     [-11, -22, -23, -24],
                                     [-12,  -6, -17, -18],
                                     [-13,  -7,  -7, -18],
                                     [-14,  -8,  -8,  -6]])

    assert np.array_equal(Alignment._align_matrix, align_matrix_by_hand) 
    assert np.array_equal(Alignment._gapA_matrix, gapA_matrix_by_hand)
    assert np.array_equal(Alignment._gapB_matrix, gapB_matrix_by_hand)


def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    Alignment=NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    seq3_seq4_alignemnt = Alignment.align(seq3, seq4)

    assert seq3_seq4_alignemnt == (17, 'MAVHQLIRRP', 'M---QLIRHP')




