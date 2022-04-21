"""
Currently meant to be a future replacement of Align class
"""
from Bio import SeqIO
import numpy as np
from collections import Counter
import itertools
import copy
import blosum as bl

DELETION, INSERTION, MATCH = range(3)


class MonkeyAlign():

    def __init__(self):

        # Numpy array of character arrays (with trailing gaps)
        self.seqs = []

    # this calls read_fasta
    def _read_fasta(self, files):
        """ Read in fasta file(s)

            Args:
                files (str, list[str]): path to FASTA file or list of file paths

            Return:
                self.seqs (list[ndarray]): list of amino acid sequence ndarrays
        """

        # recursively convert amino acid sequence to numpy character array for each sequence in fasta file
        if type(files) == list:
            fastas = [self._read_fasta(file) for file in files]
            self.seqs = [seq for fasta in fastas for seq in fasta]
            return self.seqs

        # parse fasta file
        fasta_sequences = SeqIO.parse(open(files), 'fasta')

        # get list of individual string characters for each seq in file and convert to list
        self.seqs = [list(str(fasta.seq)) for fasta in fasta_sequences]

        # convert sequence strings to numpy character arrays
        self.seqs = list(map(np.array, self.seqs))

        # return seqs
        return self.seqs

    def _add_trailing(self):
        """ Add trailing '*' characters to end of sequence arrays """

        # get max length of sequences
        max_seq = max([len(seq) for seq in self.seqs])

        # add trailing characters to arrays
        self.seqs = [np.concatenate(
            (seq, np.full(fill_value='*', shape=max_seq - len(seq)))) for seq in self.seqs]

    def read_fasta(self, files):
        """ User function to read in and format fasta files of amino acid sequences """

        # read in fasta file
        self._read_fasta(files)

        # add trailing sequences
        self._add_trailing()

        # convert list of arrays to 2D ndarray
        self.seqs = np.array(self.seqs)

    def load_str(self, *str_seq):
        """ Allows user to input sequences as string

        Args:
            str_seq (str): *args as many strings as the user wants

        Returns:

        """
        # set self.seqs to a list of np.arrays
        self.seqs = [np.array([char for char in current_str],
                                       dtype='<U1') for current_str in str_seq]

        # add trailing sequences so they match lengths
        self._add_trailing()

        # convert list of arrays to 2D ndarray
        self.seqs = np.array(self.seqs, dtype='<U1')

        return self.seqs


    def get_seqs(self):
        """ Retrieve array of sequences

            Return:
                self.seqs (ndarray): array of amino acid sequence arrays
        """

        return self.seqs

    # FITNESS CRITERIA for entire matrix
    def sum_pairs_score(self, matrix=bl.BLOSUM(62), gap_cost=1):
        """ Calculates score across all columns for matches, mismatches, and gaps

            Args:
                matrix (bl.BLOSUM matrix): BLOSUM matrix as eval system

            Return:
                score (int): sum of pairs score for the fasta array
        """

        # Initialize score
        score = 0

        # Iterates through each position in the alignment
        for pos in self.seqs.T[::-1]:
            # use substition matrix to score matches and mismatches
            score += sum([matrix[x + y] for i, x in enumerate(pos)
                          for j, y in enumerate(pos) if i > j])
        return score * -1

    # Modification Agent, aligns two random seqs, returns entire alignment
    def smith_waterman(self, insertion_penalty=-1, deletion_penalty=-1,
                       mismatch_penalty=-1, match_score=2):
        """
        Find the optimum local sequence alignment for the sequences `seq1`
        and `seq2` using the Smith-Waterman algorithm. Optional keyword
        arguments give the gap-scoring scheme:

        `insertion_penalty` penalty for an insertion (default: -1)
        `deletion_penalty`  penalty for a deletion (default: -1)
        `mismatch_penalty`  penalty for a mismatch (default: -1)
        `match_score`       score for a match (default: 2)

        See <http://en.wikipedia.org/wiki/Smith-Waterman_algorithm>.

        #>>> for s in smith_waterman('AGCAGACT', 'ACACACTA'): print s
        ...
        AGCAGACT-
        A-CACACTA
        """
        # get target variables from _two_rand_seqs() method
        # static_lst is a list of np arrays of all the other alignments
        seq1, seq2, static_lst = self._two_rand_seqs()

        # get the lengths
        m, n = len(seq1), len(seq2)

        # Construct the similarity matrix in p[i][j], and remember how
        # we constructed it -- insertion, deletion or (mis)match -- in
        # q[i][j].
        p = np.zeros((m + 1, n + 1))
        q = np.zeros((m + 1, n + 1))
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                deletion = (p[i - 1][j] + deletion_penalty, DELETION)
                insertion = (p[i][j - 1] + insertion_penalty, INSERTION)
                if seq1[i - 1] == seq2[j - 1]:
                    match = (p[i - 1][j - 1] + match_score, MATCH)
                else:
                    match = (p[i - 1][j - 1] + mismatch_penalty, MATCH)
                p[i][j], q[i][j] = max((0, 0), deletion, insertion, match)

        # Yield the aligned sequences one character at a time in reverse
        # order.
        def backtrack():

            i, j = m, n
            # CHANGED HERE FROM 'OR' TO 'AND'
            while i > 0 and j > 0:
                # IS THIS REDUNDANT?
                assert i >= 0 and j >= 0
                if q[i][j] == MATCH:
                    i -= 1
                    j -= 1
                    yield seq1[i], seq2[j]
                elif q[i][j] == INSERTION:
                    j -= 1
                    yield '*', seq2[j]
                elif q[i][j] == DELETION:
                    i -= 1
                    yield seq1[i], '*'
                else:
                    assert (False)

        # TRY TO SHORTEN CODE HERE
        seq1_aligned, seq2_aligned = [''.join(reversed(s)) for s in zip(*backtrack())]
        # from string to character array
        seq1_aligned = np.array([char for char in seq1_aligned])
        seq2_aligned = np.array([char for char in seq2_aligned])

        # call _combine_again method to convert to full alignment
        self._combine_again(seq1_aligned, seq2_aligned, static_lst=static_lst)
        return self


    def _two_rand_seqs(self):
        """ finds two random sequences to align, saves the other sequences as a list of np.arrays

        Returns:
            seq1
            seq2
            static_lst

        """
        # Shuffle the current alignment
        np.random.shuffle(self.seqs)

        # Select 2 random sequences (first 2 because sequences were shuffled)
        return self.seqs[0], self.seqs[1], [static for static in self.seqs[2:]]

    def _combine_again(self, seq1, seq2, static_lst):
        """ self.seqs is now the 2 new alignments and the static remainders, adjusted for new length
        """
        # first we place them all in the same list
        static_lst.append(seq1)
        static_lst.append(seq2)

        # set list to self.seqs
        self.seqs = static_lst

        # add trailing sequences
        self._add_trailing()

        # convert list of arrays to 2D ndarray
        self.seqs = np.array(self.seqs, dtype='<U1')

    def __repr__(self):
        return str(self.seqs)









def main():

    ma = MonkeyAlign()
    ma.load_str('ABC', 'ABD', 'ABCD', 'BD')
    print(ma)
    #print(ma._two_rand_seqs())
    print(ma.sum_pairs_score())
    print(ma.smith_waterman())
    print(ma.sum_pairs_score())


if __name__ == '__main__':
    main()