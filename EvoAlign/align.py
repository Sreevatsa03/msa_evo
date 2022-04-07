from Bio import SeqIO
import numpy as np
from collections import Counter
import itertools


class Align():

    def __init__(self):

        # Numpy array of character arrays (with trailing gaps)
        self.seqs = []

        # Star alignment matrix for selecting root sequence
        self.star_matrix = []

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
        """ Add trailing '-' characters to end of sequence arrays """

        # get max length of sequences
        max_seq = max([len(seq) for seq in self.seqs])

        # add trailing characters to arrays
        self.seqs = [np.concatenate(
            (seq, np.full(fill_value='-', shape=max_seq - len(seq)))) for seq in self.seqs]

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

    def sum_pairs_score(self):
        """ Calculates the sum of pairs for matches, mismatches, and gaps

            Return:
                score (int): sum of pairs score for the fasta array
        """

        # Initialize score
        score = 0

        # Iterates through each position in the alignment
        for pos in self.seqs.T[::-1]:

            # 1 for a match and then -1 if mismatch
            match_score = sum([1 if x == y else -1 for i, x in enumerate(pos)
                               for j, y in enumerate(pos) if i > j])

            # Subtracts 1 at each position if a gap is present
            gap_score = sum([-1 if '-' in (x, y) else 0 for i, x in enumerate(pos)
                             for j, y in enumerate(pos) if i > j])

            score += match_score + gap_score
        return score

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
            gap_score = sum([-1 if '-' in (x, y) else 0 for i, x in enumerate(pos)
                             for j, y in enumerate(pos) if i > j])

        # return score
        return gap_score

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
    def _seq_comparison(a, b):
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
            match_score = sum([1 if x == y else -1 for i, x in enumerate(pos)
                               for j, y in enumerate(pos) if i > j])

            # Subtracts 1 at each position if a gap is present
            gap_score = sum([-1 if '-' in (x, y) else 0 for i, x in enumerate(pos)
                             for j, y in enumerate(pos) if i > j])

            score += match_score + gap_score
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
        self.star_matrix = np.ones((len(self.seqs), len(self.seqs))) * 0.5#np.empty((len(self.seqs), len(self.seqs))) * np.nan

        # Iterate through the rows and columns of the star matrix
        for idx_row, row in enumerate(self.seqs):
            for idx_col, col in enumerate(self.seqs):

                # For non-diagonal comparisons
                if idx_row != idx_col:

                    # If score hasn't been calculated (default value is 0.5)
                    if self.star_matrix[idx_row, idx_col] == 0.5:

                        # Calculate score and assign to position and reverse position
                        self.star_matrix[idx_row, idx_col] = Align._seq_comparison(row, col)
                        self.star_matrix[idx_col, idx_row] = self.star_matrix[idx_row, idx_col]

                # Reassign diagonal score from 0.5 to 0
                else:
                    self.star_matrix[idx_row, idx_col] = 0

        return self.star_matrix

    def get_root_seq(self):
        """ Creates star alignment matrix to select root sequence
            Uses sum_pair_scores to compare sequence similarity

            Returns:
                (np.array): Sequence with the highest similarity score
                (most similar to all the sequences (root sequence))
        """
        return self.seqs[np.sum(self.get_star_matrix(), axis=1).argmax()]

    @staticmethod
    def _matrix(a, b, match_score=3, gap_cost=2):
        """ gets matrix of alignment for smith waterman algorithm

            Args:

            Returns:

        """

        # H is the np array of scores
        H = np.zeros((len(a) + 1, len(b) + 1), int)

        # itertools makes nested for loop
        for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):

            # if a and b align at this spot, add the match score to the score from left diagonal
            # if they don't, then subtract the match score from left diagonal
            match = H[i - 1, j - 1] + \
                (match_score if a[i - 1] == b[j - 1] else - match_score)

            # delete score, subtract from score immediately left
            delete = H[i - 1, j] - gap_cost

            # insert score, subtract from score above
            insert = H[i, j - 1] - gap_cost

            # set the matrix equal to the highest val (min is 0)
            H[i, j] = max(match, delete, insert, 0)

        return H

    @staticmethod
    def _traceback(H, b, b_='', old_i=0):
        # flip H to get index of **last** occurrence of H.max() with np.argmax()
        # 2 np flips (0 then 1) https://numpy.org/doc/stable/reference/generated/numpy.flip.html
        # maintains individual row order, but row positions are reversed (top row as bottom row, etc.)
        H_flip = np.flip(np.flip(H, 0), 1)

        # iterates through matrix to find the the location (i_, j_) of the max val
        i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)

        # convert i_, j_ from H_flip back into i,j from H
        # (i, j) are **last** indexes of H.max()
        i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))

        # if the max val is zero, return empty string and j
        if H[i, j] == 0:
            return b_, j

        #
        b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
        return Align._traceback(H[0:i, 0:j], b, b_, i)

    @staticmethod
    def _smith_waterman(a, b, match_score, gap_cost):
        a, b = a.upper(), b.upper()
        H = Align._matrix(a, b, match_score, gap_cost)
        b_, pos = Align._traceback(H, b)
        return pos, pos + len(b_)

    def align_smith_waterman(self, match_score=3, gap_cost=2):
        """

        """
        # create string version

        for a, b in zip(self.str_seqs, self.str_seqs[1:]):
            start, end = Align._smith_waterman(
                a, b, match_score=match_score, gap_cost=gap_cost)
            print(a[start:end])
            return a[start:end]
