from Bio import SeqIO
import numpy as np

class Align():

    def __init__(self):
        self.seqs = []

    def read_fasta(self, files):
        if type(files) == list:
            fastas = [self.read_fasta(file) for file in files]
            self.seqs = [seq for fasta in fastas for seq in fasta]
            return self.seqs

        fasta_sequences = SeqIO.parse(open(files),'fasta')
        self.seqs = [list(str(fasta.seq)) for fasta in fasta_sequences]
        self.seqs = list(map(np.array, self.seqs))
        return self.seqs