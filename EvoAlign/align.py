from Bio import SeqIO
import numpy as np
from collections import Counter
import itertools
import copy
import blosum as bl


DELETION, INSERTION, MATCH = range(3)

class Align():

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

    def get_seqs(self):
        """ Retrieve array of sequences

            Return:
                self.seqs (ndarray): array of amino acid sequence arrays
        """

        return self.seqs

    #FITNESS CRITERIA
    def sum_pairs_score(self, matrix=bl.BLOSUM(62)):
        """ Calculates the sum of pairs for matches, mismatches, and gaps

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
        return score

    #Not for evo framework
    def count_gaps(self):
        """ Calculates the sum of pairs for matches, mismatches, and gaps

            Return:
                score (int): sum of pairs score for the fasta array
        """

        # Initialize score
        gap_score = 0

        # Iterates through each position in the alignment
        for pos in self.seqs.T:

            # Subtracts 1 at each position if a gap is present
            gap_score = sum([-1 if '*' in (x, y) else 0 for i, x in enumerate(pos)
                             for j, y in enumerate(pos) if i > j])

        # return score
        return gap_score

    #Not for evo framework
    def count_matches(self):
        """ Calculates the sum of pairs for matches, mismatches, and gaps

            Return:
                score (int): sum of pairs score for the fasta array
        """

        # Initialize score
        match_score = 0

        # Iterates through each position in the alignment
        for pos in self.seqs.T:

            # 1 for a match and then -1 if mismatch
            match_score = sum([1 if x == y else -1 for i, x in enumerate(pos)
                               for j, y in enumerate(pos) if i > j])

        # return score
        return match_score

    # https://tiefenauer.github.io/blog/smith-waterman/

    @staticmethod
    def _seq_comparison(a, b, matrix=bl.BLOSUM(62)):
        """ Calculates the comparison score of the inputted sequences

            Args:
                a (np.array): first sequence
                b (np.array): second sequence

            Return:
                score (int): comparison score for inputted sequences
        """

        seq_array = np.array([a, b])

        # Initialize score
        score = 0

        # Iterates through each position in the alignment
        for pos in seq_array.T[::-1]:
            # 1 for a match and then -1 if mismatch
            score += sum([matrix[x + y] for i, x in enumerate(pos)
                               for j, y in enumerate(pos) if i > j])
        return score

    def get_star_matrix(self):
        """Gets star alignment matrix (used to get root sequence)

            Returns: self.star_matrix (np.array): star alignment matrix
        """

        # Check to make sure all sequences are the same length
        # https://stackoverflow.com/questions/35791051/better-way-to-check-if-all-lists-in-a-list-are-the-same-length
        it = iter(self.seqs)
        the_len = len(next(it))
        if not all(len(l) == the_len for l in it):
            raise ValueError('not all sequences have same length!')

        # Initialize the matrix (with a non-integer value that _seq_comparison won't return)
        # np.empty((len(self.seqs), len(self.seqs))) * np.nan
        star_matrix = np.ones((len(self.seqs), len(self.seqs))) * 0.5

        # Iterate through the rows and columns of the star matrix
        for idx_row, row in enumerate(self.seqs):
            for idx_col, col in enumerate(self.seqs):

                # For non-diagonal comparisons
                if idx_row != idx_col:

                    # If score hasn't been calculated (default value is 0.5)
                    if star_matrix[idx_row, idx_col] == 0.5:

                        # Calculate score and assign to position and reverse position
                        star_matrix[idx_row, idx_col] = self._seq_comparison(
                            row, col)
                        star_matrix[idx_col,
                                    idx_row] = star_matrix[idx_row, idx_col]

                # Reassign diagonal score from 0.5 to 0
                else:
                    star_matrix[idx_row, idx_col] = 0

        return star_matrix

    def get_root_seq(self):
        """ Creates star alignment matrix to select root sequence
            Uses sum_pair_scores to compare sequence similarity

            Returns:
                (np.array): Sequence with the highest similarity score
                (most similar to all the sequences (root sequence))
        """
        mat = self.get_star_matrix()

        root_idx = np.sum(mat, axis=1).argmax()

        if self.get_star_matrix()[root_idx].argmax() == root_idx:
            max = np.sort(mat[root_idx])[-2]
            max_idx = list(self.get_star_matrix()[root_idx]).index(max)

        return root_idx, max_idx

    # Agent
    def smith_waterman(self, seq1, seq2, insertion_penalty = -1, deletion_penalty = -1,
                    mismatch_penalty = -1, match_score = 2):
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
            while i > 0 or j > 0:
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
                    assert(False)

        # seq1_aligned, seq2_aligned = [''.join(reversed(s)) for s in zip(*backtrack())]

        # return np.array(list(seq1_aligned)), np.array(list(seq2_aligned))

        return [''.join(reversed(s)) for s in zip(*backtrack())]

    def sw(self, seq1, seq2, insertion_penalty = -1, deletion_penalty = -1, mismatch_penalty = -1, match_score = 2):
        seq1_aligned, seq2_aligned = self.smith_waterman(seq1, seq2, insertion_penalty, deletion_penalty,
                    mismatch_penalty, match_score)

        return np.array(list(seq1_aligned)), np.array(list(seq2_aligned))

