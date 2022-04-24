"""
Currently meant to be a future replacement of Align class
"""
from Bio import SeqIO
import numpy as np
from collections import Counter
import itertools
import copy
import blosum as bl
import random as rnd

DELETION, INSERTION, MATCH = range(3)
BLOSUM_MATRICES = {45: bl.BLOSUM(45), 50: bl.BLOSUM(50), 62: bl.BLOSUM(62), 80: bl.BLOSUM(80), 90: bl.BLOSUM(90)}

class MonkeyAlign():

    def __init__(self):

        # Numpy array of character arrays (with trailing gaps)
        self.seqs = []
        self.num_seqs = 0

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

        # update number of seqs
        self.num_seqs = len(self.seqs)

        # return seqs
        return self.seqs

    def load_str(self, *str_seq):
        """ Allows user to input sequences as string

        Args:
            str_seq (str): *args as many strings as the user wants

        Returns:

        """
        # set self.seqs to a list of np.arrays
        self.seqs = [np.array([char for char in current_str],
                                       dtype='<U1') for current_str in str_seq]

        # update number of seqs
        self.num_seqs = len(self.seqs)

        # add trailing sequences so they match lengths
        self._add_trailing()

        # convert list of arrays to 2D ndarray
        self.seqs = np.array(self.seqs, dtype='<U1')

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


    # FITNESS CRITERIA
    def match_count(self):
        """ Calculates the number of matches

            Return:
                score (int): sum of pairs score for the fasta array
        """

        # Initialize score
        match_score = 0

        # Iterates through each position in the alignment
        for pos in self.seqs.T:

            # Subtracts 1 at each position if a gap is present
            match_score += sum([1 if (x == y and x != '*') else 0
                             for i, x in enumerate(pos)
                             for j, y in enumerate(pos) if i > j])

        # return score
        return match_score

    """# assists a fitness criteria
    def _split_align(self):
        ''' splits the alignment into first and second half to run sum pairs score
        we want to get the seq

        Returns:

        '''
        first_half = self.seqs[:int(0.5*len(self.seqs))]
        second_half = self.seqs[int(0.5*len(self.seqs)):]
        
        return first_half, second_half"""

    def get_seqs(self):
        """ Retrieve array of sequences

            Return:
                self.seqs (ndarray): array of amino acid sequence arrays
        """

        return self.seqs

    # FITNESS CRITERIA for entire matrix
    def sum_pairs_score(self, matrix=bl.BLOSUM(62)):
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
        return score

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

        # from string to character array - MAKE FUNCTIONAL
        seq1_aligned = np.array([char for char in seq1_aligned])
        seq2_aligned = np.array([char for char in seq2_aligned])

        # call _combine_again method to convert to full alignment
        self._combine_again(seq1_aligned, seq2_aligned, static_lst=static_lst)

        # return self
        return self

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

        # modification agent

    @staticmethod
    def _get_begin_trail(seq):
        """ gets you the index of the last character before trailing seqs begin"""

        for i in range(len(seq)-1, -1, -1):
            if seq[i] != '*':
                return i
        # just in case that the if statement is never true
        return len(seq)-1

    def move_end_gap(self):
        """ Take a trailing sequence and randomly insert it for one sequence """

        # get the index for a random sequence
        rand_seq_idx = rnd.randint(0, self.num_seqs - 1)

        # if there is at least one gap at the end
        if self.seqs[rand_seq_idx][-1] == '*':
            # pick a random index for the position, not including the other trailing gaps
            rand_pos_idx = rnd.randint(0, MonkeyAlign._get_begin_trail(self.seqs[rand_seq_idx]))

            # must use temp_storage because we are changing the length
            temp = self.seqs[rand_seq_idx][:-1]

            temp = np.insert(temp, rand_pos_idx, '*')

            # reset the seq to the new one
            self.seqs[rand_seq_idx] = temp

            # must return a MonkeyAlign object
            return self

        # if no trailing gaps
        else:
            return self

    def __repr__(self):
        return str(self.seqs)

