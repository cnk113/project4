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
    needle = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    x,y,z = needle.align(seq1, seq2) 
    assert needle._align_matrix[0, 0] == 0, "Initial match isn't 0"
    assert np.all(needle._back[0, :] == np.inf), "Not inf backtraces for initial row"
    assert np.all(needle._back[:, 0] == np.inf), "Not inf backtraces for initial column"


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
    needle = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    x,y,z = needle.align(seq3, seq4)
    # Tests last position of backtrace
    mx_list = [needle._align_matrix[-1, -1], needle._gapA_matrix[-1, -1], needle._gapB_matrix[-1, -1]]
    idx = np.argmax(mx_list)
    if idx == 0:
        assert y[-1] != "-" and z[-1] != "-", "Last position is incorrect"
    elif idx == 1:
        assert y[-1] == "-", "Last position is incorrect"
    else:
        assert z[-1] == "-", "Last position is incorrect"
    assert x == 18, "Score is incorrect"


