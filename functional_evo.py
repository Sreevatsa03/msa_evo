"""
getter functions for agents in Evo
"""

from EvoAlign import Align, Evo, MonkeyAlign
import itertools
import numpy as np
import blosum as bl
import random as rnd

def smith_waterman_align(curr_align):
    """ takes in MonkeyAlignment object

    Args:
        curr_align:

    Returns:

    """
    # assert?
    curr_align = curr_align[0]
    return curr_align.smith_waterman()

# fitness agent
def blosum_62_score(curr_align):
    return curr_align.sum_pairs_score(matrix=bl.BLOSUM(62), gap_cost=1)

def blosum_50_score(curr_align):
    return curr_align.sum_pairs_score(matrix=bl.BLOSUM(45), gap_cost=1)


def main():
    # path to all.fasta data
    all = 'data_sources/P53_test_data/all.fasta'

    # create instance of Align
    initial_object = MonkeyAlign()

    # read fasta into Align
    initial_object.read_fasta(all)

    #ma = monkey_align.MonkeyAlign()
    #ma.load_str('ABCF', 'ABDEH', 'ABCD', 'BD', 'AB')

    # create Evo environment
    E1 = Evo()

    # register fitness criteria
    E1.add_fitness_criteria('blosum_62_score', blosum_62_score)
    E1.add_fitness_criteria('blosum_50_score', blosum_50_score)

    # add modification agent
    E1.add_agent('smith_waterman', smith_waterman_align)

    # add initial solution
    E1.add_solution(initial_object)
    #E1.add_solution(ma)

    E1.evolve(10, dom=100, status=3)

    print(E1.pop)
    print(E1.fitness)
    print(E1.agents)






if __name__ == '__main__':
    main()

