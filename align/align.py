# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        # TODO: Fill in the Needleman-Wunsch Algorithm below
        to perform global sequence alignment of seqA and seqB
        and return a tuple with the following format
        (alignment score, seqA alignment, seqB alignment)
        Also, write up a docstring for this function using the
        _read_sub_matrix as an example.
        Don't forget to comment your code!
        """
        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # TODO Implement the global sequence alignment here
        for a in range(len(seqA)+1): # num rows
            for b in range(len(seqB)+1): # num cols
                if a==0:
                    self._gapA_matrix[a,b] = self.gap_open + b*self.gap_extend # first row/col
                    if b==0:
                        self._align_matrix[a,b] = 0 # set upper left corner to zero
                elif a>0:
                    if b>0:
                        # determine the match score
                        char_a = seqA[a-1] 
                        char_b = seqB[b-1] 
                        match = self.sub_dict[(char_a, char_b)]

                        # calculate maximums for M, A, and B matrices, storing options
                        m_options = [self._align_matrix[a-1,b-1], self._gapA_matrix[a-1,b-1], self._gapB_matrix[a-1,b-1]]
                        m_max = max(m_options)
                        a_options = [self.gap_open + self.gap_extend + self._align_matrix[a,b-1],
                                    self.gap_extend + self._gapA_matrix[a,b-1],
                                    self.gap_open + self.gap_extend + self._gapB_matrix[a,b-1]]
                        a_max = max(a_options)
                        b_options = [self.gap_open + self.gap_extend + self._align_matrix[a-1,b],
                                    self.gap_open + self.gap_extend + self._gapA_matrix[a-1,b],
                                    self.gap_extend + self._gapB_matrix[a-1,b]]
                        b_max = max(b_options)

                        # update pointers
                        self._back[a,b] = m_options.index(m_max)
                        self._back_A[a,b] = a_options.index(a_max)
                        self._back_B[a,b] = b_options.index(b_max)

                        # calculate values for each align/gap matrix
                        self._align_matrix[a,b] = match + m_max
                        self._gapA_matrix[a,b] = a_max
                        self._gapB_matrix[a,b] = b_max
                if b==0:
                    self._gapB_matrix[a,b] = self.gap_open + a*self.gap_extend # first row/col

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        # TODO Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.
        """
        # Implement this method based upon the heuristic chosen in the align method above.
        
        # start with bottom right square
        curr_a = len(self._seqA)
        curr_b = len(self._seqB)
        
        # trace path from botom right to top left
        while curr_a>0 and curr_b>0:
            # determine maximum and index of maximum in cell
            curr_options = [self._align_matrix[curr_a, curr_b],
                             self._gapA_matrix[curr_a, curr_b],
                             self._gapB_matrix[curr_a, curr_b]]
            curr_max = max(curr_options)
            curr_max_index = curr_options.index(curr_max)
            
            # assign alignment score from bottom right
            if curr_a==len(self._seqA) and curr_b==len(self._seqB):
                align_score = curr_max
            
            # match/mismatch - add to sequences, set new cell
            if curr_max_index==0:
                self.seqA_align = self.seqA_align + self._seqA[curr_a-1]
                self.seqB_align = self.seqB_align + self._seqB[curr_b-1]
                curr_a -= 1
                curr_b -= 1
            
            # skip in sequence A
            elif curr_max_index==1:
                self.seqA_align = self.seqA_align + "-"
                self.seqB_align = self.seqB_align + self._seqB[curr_b-1]
                curr_b -= 1
            
            # skip in sequence B
            elif curr_max_index==2:
                self.seqA_align = self.seqA_align + self._seqA[curr_a-1]
                self.seqB_align = self.seqB_align + "-"
                curr_a -= 1
            
        # reverse strings
        self.seqA_align = self.seqA_align[::-1]
        self.seqB_align = self.seqB_align[::-1]
        
        return (align_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
