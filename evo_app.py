from re import A
from EvoAlign import Align
import itertools
import numpy as np
import blosum as bl

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

# print(bl.BLOSUM(62))

# prints correct scoring matrix from Wikipedia example
# print(align._matrix('GGTTGACTA', 'TGTTACGG'))

# a, b = 'GGTTGACTA', 'TGTTACGG'
a, b = align.get_seqs()[0], align.get_seqs()[1]
# H = align._matrix(a, b)
# print(H)
# print(align._traceback(H, b)[0])

a_align, b_align = align.align(a, b)

print(align._seq_comparison(a, b))
print(align._seq_comparison(a_align, b_align))

