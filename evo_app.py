from re import A # is this a miscopy?
from EvoAlign import align, evo
import itertools
import numpy as np
import blosum as bl
import random as rnd


dog = 'data_sources/P53_test_data/canis_lupus_familiaris.fasta'
human = 'data_sources/P53_test_data/homo_sapiens.fasta'
mouse = 'data_sources/P53_test_data/mus_musculus.fasta'
all = 'data_sources/P53_test_data/all.fasta'
fastas = [dog, human, all]

alignment = align.Align()
alignment.read_fasta(all)

#print(alignment.get_seqs())

#print(alignment.sum_pairs_score())
#print(alignment.get_root_seq())

# alignment.clustal()


# prints correct scoring matrix from Wikipedia example
# print(alignment._matrix('GGTTGACTA', 'TGTTACGG'))

print(alignment.smith_waterman())


# a, b = 'GGTTGACTA', 'TGTTACGG'
# start, end = align._smith_waterman(a, b)
# print(a[start:end])  # GTTGAC
