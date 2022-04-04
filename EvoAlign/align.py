from Bio import SeqIO
import numpy as np
from collections import Counter

class Align():

    def __init__(self):
        self.seqs = []

    def _read_fasta(self, files):
        """ Read in fasta file(s)

            Args:
                files (str, list[str]): path to FASTA file or list of file paths

            Return:
                self.seqs (list[ndarray]): list of amino acid sequence ndarrays
        """

        # recursively convert amino acid sequence to numpy character array for each sequence in fasta file
        if type(files) == list:
            fastas = [self.read_fasta(file) for file in files]
            self.seqs = [seq for fasta in fastas for seq in fasta]
            return self.seqs

        # parse fasta file
        fasta_sequences = SeqIO.parse(open(files),'fasta')

        # get string for each seq in file and convert to list
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
        self.seqs = [np.concatenate((seq, np.full(fill_value='-', shape=max_seq - len(seq)))) for seq in self.seqs]

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

        # If amino acids are the same, add 1; if different, subtract 1 from score
        for pos in self.seqs.T:
            score += sum([1 if x == y else -1 for i,x in enumerate(pos) for j,y in enumerate(pos) if i > j])

        # return score
        return score