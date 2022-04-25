"""
getter functions for agents in Evo
"""

from colorama import init
from EvoAlign import Align, Evo, MonkeyAlign, AminoAcid
import itertools
import numpy as np
import blosum as bl
import random as rnd
import os

BLOSUM = {45: bl.BLOSUM(45), 50: bl.BLOSUM(50), 62: bl.BLOSUM(62), 80: bl.BLOSUM(80), 90: bl.BLOSUM(90), 'HYDRO': AminoAcid().hydropthy_dict, 'VOL': AminoAcid().volume_dict}

# agent


def smith_waterman_half_align(curr_align):
    """ takes in MonkeyAlign object

    Args:
        curr_align:

    Returns:

    """
    # ASSERT HERE?
    curr_align = curr_align[0]
    return curr_align.smith_waterman(0.5)


def smith_waterman_quarter_align(curr_align):
    """ takes in MonkeyAlign object

    Args:
        curr_align:

    Returns:

    """
    # ASSERT HERE?
    curr_align = curr_align[0]
    return curr_align.smith_waterman(0.25)


def smith_waterman_eighth_align(curr_align):
    """ takes in MonkeyAlign object

    Args:
        curr_align:

    Returns:

    """
    # ASSERT HERE?
    curr_align = curr_align[0]
    return curr_align.smith_waterman(0.125)

# dumb mod agent


def move_end_gap(curr_align):
    """ removes a trailing sequence of one randomly selected seqeunce and randomly adds in a gap (maintains length)
    takes in MonkeyAlign object and calls .move_end_gap() method
    """
    curr_align = curr_align[0]
    return curr_align.move_end_gap()


# fitness criteria
def blosum_62_score(curr_align):
    return curr_align.sum_pairs_score(matrix=BLOSUM[62])


def blosum_45_score(curr_align):
    return curr_align.sum_pairs_score(matrix=BLOSUM[45])


def blosum_90_score(curr_align):
    return curr_align.sum_pairs_score(matrix=BLOSUM[90])


def hydropathy_score(curr_align):
    return curr_align.sum_pairs_score(matrix=BLOSUM['HYDRO'])


def volume_score(curr_align):
    return curr_align.sum_pairs_score(matrix=BLOSUM['VOL'])


# def match_count(curr_align):
#     return curr_align.match_count()


def main():
    if os.path.exists("solutions.dat"):
        os.remove("solutions.dat")

    # path to all.fasta data
    # dash = 'data_sources/P53_test_data/dash.fasta'
    all = 'data_sources/P53_test_data/all.fasta'
    # dog = 'data_sources/P53_test_data/canis_lupus_familiaris.fasta'
    # human = 'data_sources/P53_test_data/homo_sapiens.fasta'
    # mouse = 'data_sources/P53_test_data/mus_musculus.fasta'

    # create instance of Align
    initial_object = MonkeyAlign()

    # read fasta into Align
    initial_object.read_fasta(all)

    # initial_object._split_seqs(initial_object.seqs[0], initial_object.seqs[1], 0.25)
    initial_object.smith_waterman(0.5)

    # create Evo environment
    E1 = Evo()

    # register fitness criteria
    E1.add_fitness_criteria('blosum_62_score', blosum_62_score)
    E1.add_fitness_criteria('hydropathy_score', hydropathy_score)
    E1.add_fitness_criteria('volume_score', volume_score)
    # E1.add_fitness_criteria('blosum_45_score', blosum_45_score)
    # E1.add_fitness_criteria('blosum_90_score', blosum_90_score)
    # E1.add_fitness_criteria('match_count', match_count)

    # add modification agent
    E1.add_agent('smith_waterman', smith_waterman_half_align)
    E1.add_agent('smith_waterman', smith_waterman_quarter_align)
    E1.add_agent('smith_waterman', smith_waterman_eighth_align)
    # E1.add_agent('move_end_gap', move_end_gap)

    # add initial solution
    E1.add_solution(initial_object)

    E1.evolve(1000, dom=100, status=100)

    E1.save_solutions()
    E1.visualize()


if __name__ == '__main__':
    main()
