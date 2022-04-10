from EvoAlign import Align
import itertools
import numpy as np

dog = 'data_sources/P53_test_data/canis_lupus_familiaris.fasta'
human = 'data_sources/P53_test_data/homo_sapiens.fasta'
mouse = 'data_sources/P53_test_data/mus_musculus.fasta'
all = 'data_sources/P53_test_data/all.fasta'
fastas = [dog, human, all]

align = Align()
align.read_fasta(all)

#print(align.get_seqs())

# print(align.sum_pairs_score())
# print(align.get_root_seq())

# align.clustal()


# prints correct scoring matrix from Wikipedia example
# print(align._matrix('GGTTGACTA', 'TGTTACGG'))

# a, b = 'GGTTGACTA', 'TGTTACGG'
# H = align._matrix(a, b)
# print(align._traceback(H, b))  # ('gtt-ac', 1)

a, b = 'GGTTGACTA', 'TGTTACGG'
start, end = align._smith_waterman(a, b)
print(a[start:end])  # GTTGAC
