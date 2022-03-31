from Bio import SeqIO
import numpy as np
from collections import Counter

class Align():

    def __init__(self):
        self.seqs = []

    # go back and add trailing sequences --> so they are same length
    def read_fasta(self, files):
        if type(files) == list:
            fastas = [self.read_fasta(file) for file in files]
            self.seqs = [seq for fasta in fastas for seq in fasta]
            return self.seqs

        fasta_sequences = SeqIO.parse(open(files),'fasta')
        self.seqs = [list(str(fasta.seq)) for fasta in fasta_sequences]
        self.seqs = list(map(np.array, self.seqs))
        return self.seqs

    def sum_pairs_score(self):
        """ Calculates the sum of pairs for matches, mismatches, and gaps

            Args:
                fasta_array (np.array): character array with each amino acid
                as an element

            Return:
                 score (int): sum of pairs score for the fasta array
        """

        # Initialize score
        score = 0

        # Construct column-wise comparison
        for i in range(min([len(seq) for seq in self.seqs])):
            column = [self.seqs[j][i] for j in range(len(self.seqs))]

            # If amino acids are the same, add 1; if different, subtract 1 from score
            for i in range(len(column)):
                for j in range(i+1, len(column)):
                    if column[i] == column[j]:
                        score += 1
                    else:
                        score -= 1
        return score


