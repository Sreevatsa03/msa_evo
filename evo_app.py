 
from re import A # is this a miscopy?
from EvoAlign import Align, MonkeyAlign
import itertools
import numpy as np
import blosum as bl
import random as rnd
 
 
dog = 'data_sources/P53_test_data/canis_lupus_familiaris.fasta'
human = 'data_sources/P53_test_data/homo_sapiens.fasta'
mouse = 'data_sources/P53_test_data/mus_musculus.fasta'
all = 'data_sources/P53_test_data/all.fasta'
fastas = [dog, human, all]
 
align = Align()
align.read_fasta(all)

seq1, seq2 = align.sw(align.get_seqs()[0], align.get_seqs()[1])

print(align._seq_comparison(seq1, seq2))
 
# Test MonkeyAlign
# ma = MonkeyAlign()
# ma.load_str('ABC', 'ABD', 'ABCD', 'BD')
# print(ma)
# #print(ma._two_rand_seqs())
# print(ma.sum_pairs_score())
# print(ma.smith_waterman())
# print(ma.sum_pairs_score())

