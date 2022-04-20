from EvoAlign import align, evo
import itertools
import numpy as np
import blosum as bl
import random as rnd

DELETION, INSERTION, MATCH = range(3)

def blosum_62_score(current_alignment):
    """ Evaluates the alignment score for the current_alignment as a whole
        using the BLOSUM 62 Matrix

        Return:
            score (int): BLOSUM 62 Matrix score
    """

    # Initialize score
    score = 0

    # Get BLOSUM 62 Matrix
    matrix = bl.BLOSUM(62)

    # Iterates through each position in the alignment
    for pos in current_alignment.T[::-1]:
        # use substition matrix to score matches and mismatches
        score += sum([matrix[x + y] for i, x in enumerate(pos)
                      for j, y in enumerate(pos) if i > j])
    return score


def smith_waterman(current_alignment, insertion_penalty = -1, deletion_penalty = -1,
                mismatch_penalty = -1, match_score = 2):

    current_alignment = current_alignment[0]

    # Shuffle the current alignment
    np.random.shuffle(current_alignment)

    # Select 2 random sequences (first 2 because sequences were shuffled)
    seq1, seq2 = current_alignment[0], current_alignment[1]

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
        while i > 0 or j > 0:
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
                assert(False)

    seq1_aligned, seq2_aligned = [''.join(reversed(s)) for s in zip(*backtrack())]

    return np.array(list(seq1_aligned)), np.array(list(seq2_aligned))



def main():
    # path to all.fasta data
    all = 'data_sources/P53_test_data/all.fasta'

    # create instance of Align
    initial_object = align.Align()

    # read fasta into Align
    initial_object.read_fasta(all)

    # getter function to see seqs
    a = initial_object.get_seqs()

    # create Evo environment
    E = evo.Evo()

    # register fitness criteria
    E.add_fitness_criteria('blosum_62_score', blosum_62_score)

    # add modification agent
    E.add_agent('smith_waterman', smith_waterman)

    # add initial solution
    E.add_solution(a)

    E.evolve(10, dom=0, status=1)



if __name__ == '__main__':
    main()

