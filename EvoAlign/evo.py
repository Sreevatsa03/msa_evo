"""
Created on Sun Nov  1 19:48:48 2020
@author: John Rachlin
@file: evo_v4.py: An evolutionary computing framework (version 4)
Assumes no Solutions class.
"""
import pickle
import random as rnd
import copy
from functools import reduce
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class Evo:

    def __init__(self):
        """ Population constructor """
        self.pop = {}  # The solution population eval -> solution
        self.fitness = {}  # Registered fitness functions: name -> objective function
        # Registered agents:  name -> (operator, num_solutions_input)
        self.agents = {}

    def size(self):
        """ The size of the current population """
        return len(self.pop)

    def add_fitness_criteria(self, name, f):
        """ Register a fitness criterion (objective) with the
        environment. Any solution added to the environment is scored
        according to this objective """
        self.fitness[name] = f

    def add_agent(self, name, op, k=1):
        """ Register a named agent with the population.
        The operator (op) function defines what the agent does.
        k defines the number of solutions the agent operates on. """
        self.agents[name] = (op, k)

    def add_solution(self, sol):
        """ Add a solution to the population """
        # eval = ((obj1, score1), (obj2, score2).....)
        eval = tuple((name, f(sol)) for name, f in self.fitness.items())
        self.pop[eval] = sol

    def run_agent(self, name):
        """ Invoke an agent against the population """
        op, k = self.agents[name]
        picks = self.get_random_solutions(k)
        new_solution = op(picks)
        self.add_solution(new_solution)

    def evolve(self, n=1, dom=100, status=100, sync=1000):
        """ Run n random agents (default=1)
        dom defines how often we remove dominated (unfit) solutions
        status defines how often we display the current population """

        agent_names = list(self.agents.keys())
        for i in range(n):
            pick = rnd.choice(agent_names)
            self.run_agent(pick)

            if i % sync == 0:
                # merge the saved solutions into my populations
                try:
                    with open('solutions.dat', 'rb') as file:
                        loaded = pickle.load(file)
                        for eval, sol in loaded.items():
                            self.pop[eval] = sol
                except Exception as e:
                    print(e)

                # remove dominated
                self.remove_dominated()

                # resave the population back to disk
                with open('solutions.dat', 'wb') as file:
                    pickle.dump(self.pop, file)

            if i % dom == 0:
                self.remove_dominated()

            if i % status == 0:
                self.remove_dominated()
                print("Iteration:", i)
                print("Population size:", self.size())
                print(self)

        # Clean up the population
        self.remove_dominated()

    def get_random_solutions(self, k=1):
        """ Pick k random solutions from the population """
        if self.size() == 0:
            return []
        else:
            solutions = tuple(self.pop.values())
            return [copy.deepcopy(rnd.choice(solutions)) for _ in range(k)]

    @staticmethod
    def _dominates(p, q):
        """ p = evaluation of solution: ((obj1, score1), (obj2, score2), ... )"""
        pscores = [score for _, score in p]
        qscores = [score for _, score in q]
        score_diffs = list(map(lambda x, y: y - x, pscores, qscores))
        min_diff = min(score_diffs)
        max_diff = max(score_diffs)
        return min_diff <= 0.0 and max_diff < 0.0

    @staticmethod
    def _reduce_nds(S, p):
        return S - {q for q in S if Evo._dominates(p, q)}

    def remove_dominated(self):
        """ Remove dominated solutions """
        nds = reduce(Evo._reduce_nds, self.pop.keys(), self.pop.keys())
        self.pop = {k: self.pop[k] for k in nds}

    def __str__(self):
        """ Output the solutions in the population """
        rslt = ""
        for eval, sol in self.pop.items():
            rslt += str(dict(eval)) + ":\t" + str(sol) + "\n"
        return rslt

    def _data_to_dict(self):
        """ Convert fitness criteria and solutions to dict """

        # Gets a list of solutions and fitness criteria
        solutions = list(self.pop.keys())
        fitness = list(self.fitness.keys())

        # Creates a dictionary where key is criteria and value is score list
        return {key: [[fit[1] for fit in val if fit[0] == key][0]
                      for val in solutions]
                for key in fitness}

    def save_solutions(self):
        """ Save the solutions in the population in a csv file """

        output_dict = self._data_to_dict()

        df = pd.DataFrame(output_dict)
        df.to_csv('solutions.csv')

    def visualize(self, axes=(0, 1, 2)):
        """ Create two visualizations to show the tradeoffs between agents: 3D scatterplot and pairplot """

        # plot 3D scatterplot and pairplot
        sns.set(style="darkgrid")

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Sets the axes based on user preference or default and makes dict
        fit = list(self.fitness.keys())
        axes = [fit.index(num) if num in fit else num for num in axes]
        dims_dict = self._data_to_dict()

        # Sets the labels for the scatter plot
        ax.set_xlabel(f'{fit[axes[0]]}')
        ax.set_ylabel(f'{fit[axes[1]]}')
        ax.set_zlabel(f'{fit[axes[2]]}')

        # Sets the dimensions for the scatter plot
        x = dims_dict[fit[axes[0]]]
        y = dims_dict[fit[axes[1]]]
        z = dims_dict[fit[axes[2]]]

        # Sets the title for the scatter plot
        plt.title(f'{fit[axes[0]]} vs {fit[axes[1]]} vs {fit[axes[2]]}')

        # Creates the scatter plot and shows it
        ax.scatter(x, y, z)
        plt.savefig("3D_scatter.png")

        # plot pairplot
        df = pd.read_csv('solutions.csv').drop(columns=['Unnamed: 0'])
        sns.pairplot(data=df)
        plt.savefig('pairplot.png')
