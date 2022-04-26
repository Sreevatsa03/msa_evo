from colorama import init
from EvoAlign import Align, Evo, AminoAcid
import blosum as bl
from argparse import ArgumentParser

MATRICES = {45: bl.BLOSUM(45), 50: bl.BLOSUM(50), 62: bl.BLOSUM(62), 80: bl.BLOSUM(80), 90: bl.BLOSUM(90), 'HYDRO': AminoAcid().hydropthy_dict, 'VOL': AminoAcid().volume_dict}

class EvoAlign():

    def __init__(self):
        self.Evo = Evo()
        self.Align = Align()

    @staticmethod
    def smith_waterman_half_align(curr_align):
        """ 

        Args:
            curr_align:

        Returns:

        """

        curr_align = curr_align[0]
        return curr_align.smith_waterman(0.5)

    @staticmethod
    def smith_waterman_quarter_align(curr_align):
        """ 

        Args:
            curr_align:

        Returns:

        """

        curr_align = curr_align[0]
        return curr_align.smith_waterman(0.25)

    @staticmethod
    def smith_waterman_eighth_align(curr_align):
        """ 

        Args:
            curr_align:

        Returns:

        """

        curr_align = curr_align[0]
        return curr_align.smith_waterman(0.125)

    @staticmethod
    def blosum_62_score(curr_align):
        return curr_align.sum_pairs_score(matrix=MATRICES[62])

    @staticmethod
    def blosum_45_score(curr_align):
        return curr_align.sum_pairs_score(matrix=MATRICES[45])

    @staticmethod
    def blosum_90_score(curr_align):
        return curr_align.sum_pairs_score(matrix=MATRICES[90])

    @staticmethod
    def hydropathy_score(curr_align):
        return curr_align.sum_pairs_score(matrix=MATRICES['HYDRO'])

    @staticmethod
    def volume_score(curr_align):
        return curr_align.sum_pairs_score(matrix=MATRICES['VOL'])

    def read_fasta(self, files):
        """ Read in fasta file(s)

            Args:
                files (str, list[str]): path to FASTA file or list of file paths

            Return:
                self.seqs (list[ndarray]): list of amino acid sequence ndarrays
        """

        self.Align.read_fasta(files)

    def _run(self, gens=1000, dom=100, status=100):
        # register fitness criteria
        self.Evo.add_fitness_criteria('blosum_62_score', EvoAlign.blosum_62_score)
        self.Evo.add_fitness_criteria('hydropathy_score', EvoAlign.hydropathy_score)
        self.Evo.add_fitness_criteria('volume_score', EvoAlign.volume_score)

        # add modification agent
        self.Evo.add_agent('smith_waterman', EvoAlign.smith_waterman_half_align)
        self.Evo.add_agent('smith_waterman', EvoAlign.smith_waterman_quarter_align)
        self.Evo.add_agent('smith_waterman', EvoAlign.smith_waterman_eighth_align)

        # add initial solution
        self.Evo.add_solution(self.Align)

        # evolve population
        self.Evo.evolve(gens, dom, status)

    def fitness_criteria(self):
        """ Show the fitness criteria used in evolution """
        print(list(self.Evo.fitness.keys()))

    def save_alignment(self, rankings=[]):
        """ Save a best alignment solution based on ranking of fitness criteria

            Args:
                rankings = (list(str)): list of fitness criteria in order of importance to alignment (most important to least important)
        """
        
        self.Evo.save_to_fasta(rankings)

    def align(self, gens=1000, dom=100, status=100, show=True):
        """ Aligns given amino acid sequences

            Args:
                gens (int): number of generations to evolve population
                dom (int): frequency at which solutions are dominated
                status (int): frequency at which the population is outputted
                show (bool): should the visualizations of tradeoffs be made
                save (bool): should the solution be saved to a fasta file
        """

        self._run(gens, dom, status)

        if show:
            self.Evo.visualize()
