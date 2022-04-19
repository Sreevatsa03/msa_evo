from re import A
from EvoAlign import align, evo
import itertools
import numpy as np


def blosum_62_score(current_alignment):
    """ Evaluates the alignment score for the current_alignment as a whole
        using the BLOSUM 62 Matrix

        Return:
            score (int): BLOSUM 62 Matrix score
    """

    # Initialize score
    score = 0

    # Iterates through each position in the alignment
    for pos in current_alignment.T[::-1]:
        # 1 for a match and then -1 if mismatch
        match_score = sum([1 if x == y else -1 for i, x in enumerate(pos)
                           for j, y in enumerate(pos) if i > j])

        # Subtracts 1 at each position if a gap is present
        gap_score = sum([-1 if '-' in (x, y) else 0 for i, x in enumerate(pos)
                         for j, y in enumerate(pos) if i > j])

        score += match_score + gap_score

    return score

def main():
    # path to all.fasta data
    all = 'data_sources/P53_test_data/all.fasta'

    # create instance of Align
    current_alignment = align.Align()

    # read fasta into Align
    current_alignment.read_fasta(all)

    # getter function to see seqs
    a = current_alignment.get_seqs()

    print(blosum_62_score(a))

    # create Evo environment
    #E = evo.Evo()

    # register fitness criteria
    #E.add_fitness_criteria('blosum_62_score', blosum_62_score)

    # add modification agent
    #E.add_agent('smith-waterman', align.smith_waterman)

    # add initial solution



if __name__ == '__main__':
    main()












"""
dog = 'data_sources/P53_test_data/canis_lupus_familiaris.fasta'
human = 'data_sources/P53_test_data/homo_sapiens.fasta'
mouse = 'data_sources/P53_test_data/mus_musculus.fasta'
all = 'data_sources/P53_test_data/all.fasta'
fastas = [dog, human, all]

align = Align()
align.read_fasta(all)

print(align.get_seqs())

print(align.sum_pairs_score())
print(align.get_root_seq())

# align.clustal()


# prints correct scoring matrix from Wikipedia example
# print(align._matrix('GGTTGACTA', 'TGTTACGG'))

a, b = 'GGTTGACTA', 'TGTTACGG'
# after traceback

H = align._matrix(a, b)
print(align._traceback(H, b))

# a, b = 'GGTTGACTA', 'TGTTACGG'
# start, end = align._smith_waterman(a, b)
# print(a[start:end])  # GTTGAC
"""