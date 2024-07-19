from abc import ABC, abstractmethod
from ..SettingsRetrievable import SettingsRetrievable
import evodesign.Population as Population
from ..Prediction.Predictor import Predictor
from ..GA.Selection.Selection import Selection
from ..GA.Recombination.Recombination import Recombination
from ..GA.Mutation.Mutation import Mutation
import pandas as pd
from typing import List
import math
from ..Metrics.Metric import Metric
from ..Context import Context





class Algorithm(SettingsRetrievable, ABC):

    def _params(self) -> dict:
        return {
            'max_generations': self._max_generations,
            'population_size': self._population_size,
            'predictor': self._predictor.settings(),
            'selection': self._selection.settings(),
            'recombination': self._recombination.settings(),
            'mutation': self._mutation.settings(),
            'metrics': [ metric.settings() for metric in self._metrics ]
        }



    def __init__(self,
                 max_generations: int,
                 population_size: int,
                 predictor: Predictor,
                 selection: Selection,
                 recombination: Recombination,
                 mutation: Mutation,
                 metrics: List[Metric],
                 sort_columns: List[str],
                 sort_ascending: List[bool]
                 ) -> None:
        """
        The generic workflow of a genetic algorithm for protein design.

        Parameters
        ----------
        max_generations : int
            The maximum number of generations for which the algorithm will be 
            executed (unless another termination condition is met first).
        population_size : int
            The number of sequences in the population in each generation.
        predictor : Predictor
            The protein structure prediction algorithm that will be used.
        selection : Selection
            The parents selection operation that will be used.
        recombination : Recombination
            The parents recombination operation that will be used.
        mutation : Mutation
            The sequence mutation operation that will be used.
        metrics : List[Metric]
            The metrics to be computed for each individual in the population, and
            to be used to compute the fitness function.
        sort_columns : List[str]
            The names of the columns according to which the population table will
            be sorted. 
        sort_ascending : List[bool]
            Indicate which of the columns in `sort_columns` should be sorted by
            ascending order and which by descending order.
        
        Raise
        -----
        RuntimeError
            If `sort_columns` and `sort_ascending` are if different length.
        """
        self._max_generations = max_generations
        self._population_size = population_size
        self._predictor = predictor
        self._selection = selection
        self._recombination = recombination
        self._mutation = mutation
        self._metrics = metrics
        if len(sort_columns) != len(sort_ascending):
            raise RuntimeError
        self._sort_columns = sort_columns
        self._sort_ascending = sort_ascending
        self._context = None



    def setup(self,
              context: Context,
              workspace_root: str,
              ) -> pd.DataFrame:
        """
        Initializes the workspace and the RNG, as well as the reference backbone
        and the population before running the evolutionary algorithm.

        Parameters
        ----------
        context : Context
            The context data to be used by the algorithm during its execution.
        workspace_root : str
            Path to the folder where all the output files generated by the 
            evolutionary algorithm will be stored.

        Returns
        -------
        pandas.DataFrame
            The population from which the algorithm should begin its execution.
        """
        # set the context and initialize the workspace and the RNG
        self._context = context.duplicate()
        self._context.init_workspace(workspace_root)
        self._context.workspace.save_commit_hash()
        self._context.init_rng()

        # store the sorting columns
        self._context.sort_columns = self._sort_columns
        self._context.sort_ascending = self._sort_ascending

        # save the algorithm's settings in the workspace
        self._context.workspace.save_settings(self.settings())

        # save the target PDB in the workspace
        self._context.workspace.save_target_pdb()

        # save the target FASTA in the workspace
        self._context.workspace.save_target_fasta()

        # save the sequence restrictions in the workspace
        self._context.workspace.save_sequence_restrictions(self._context.sequence_allowed_letters)

        # check if we are restoring from a previous population
        return self._context.workspace.load_population()



    def __call__(self, population: pd.DataFrame) -> None:
        """
        Starts the execution of the evolutionary algorithm.

        Parameters
        ----------
        population : pandas.DataFrame
            The population from which the evolutionary algorithm will begin its 
            execution.
        """
        # set the limit for how many generations execute
        num_generations = math.inf
        if self._context.num_generations > 0:
            num_generations = self._context.num_generations

        # set the generations counter
        t = 0

        # load the statistics file
        stats = self._context.workspace.load_statistics()

        # check if we are starting from an empty population or we are continuing
        # from the last population of an earlier execution
        if population.empty:
            population = self.initial_population()
            self._context.workspace.save_population(population)
            self._context.workspace.save_rng_state(self._context.rng.bit_generator.state)
            t = 1

        # if we are starting fresh, compute the fitness of all the initial 
        # individuals; if we are resuming from a previous population, continue
        # calculating the fitness of the individuals with missing fitness
        # (if there are any)
        population = self.compute_initial_fitness(population)

        # compute the statistics if we have not done it yet
        curr_gen_id = population.iloc[0]['generation_id']
        if stats.empty or stats.iloc[-1]['generation_id'] != curr_gen_id:
            new_stats = self.compute_statistics(population)
            stats = pd.concat([stats, new_stats.to_frame().T], ignore_index=True)

        # save progress
        self._context.workspace.save_population(population)
        self._context.workspace.save_statistics(stats)
        self.save_statistics_graph(stats)  # TODO mover la graficacion a Statistics

        # main loop
        while True:
            # check the termination conditions
            if t == num_generations:
                break
            if population.iloc[0]['generation_id'] == self._max_generations:
                break
            if stats.iloc[-1]['sequence_identity'] >= 0.95 * self._context.sequence_length:
                break
            if self.termination(population):
                break

            # check if we are resuming from an earlier execution and the children were
            # already generated; if not, apply the evolutionary operations
            children = self._context.workspace.load_population(temporary=True)
            if children.empty:
                children = self.next_population(population)
                self._context.workspace.save_rng_state(self._context.rng.bit_generator.state)
                self._context.workspace.save_population(children, temporary=True)

            # compute the children's fitness and select the surviving children
            children = self.compute_population_fitness(children, temporary=True)
            population = self.replacement(population, children)
            new_stats = self.compute_statistics(population)
            stats = pd.concat([stats, new_stats.to_frame().T], ignore_index=True)

            # save progress
            self._context.workspace.save_population(population)
            self._context.workspace.save_statistics(stats)
            self.save_statistics_graph(stats)
            self._context.workspace.delete_temporary_population()
            t += 1



    def initial_population(self) -> pd.DataFrame:
        """
        Creates a collection of random amino acid sequences to be used as the 
        population for the very first generation of the evolutionary algorithm.

        Returns
        -------
        pandas.DataFrame
            The generated population data.
        """
        population = Population.create_random(self._population_size,
                                              self._context.rng,
                                              self._context.sequence_length,
                                              1,
                                              self._context.sequence_allowed_letters)
        population['survivor'] = True
        return population



    @abstractmethod
    def compute_initial_fitness(self, population: pd.DataFrame) -> pd.DataFrame:
        raise NotImplementedError



    @abstractmethod
    def compute_statistics(self, population: pd.DataFrame) -> pd.Series:
        raise NotImplementedError



    @abstractmethod
    def save_statistics_graph(self, stats: pd.DataFrame) -> None:
        raise NotImplementedError



    def termination(self, population: pd.DataFrame) -> bool:
        """
        Checks if any additional termination conditions have been met.

        Parameters
        ----------
        population : pandas.DataFrame
            The population of the current generation.

        Returns
        -------
        bool
            A flag that indicates if any termination condition has been met or not.
        """
        return False



    def next_population(self, population: pd.DataFrame) -> pd.DataFrame:
        """
        Apply the selection, recombination and mutation operators to the given 
        population to create a new population for the next generation of the
        evolutionary algorithm.

        Parameters
        ----------
        population : pandas.DataFrame
            The population from which the new population will be created.

        Returns
        -------
        pandas.DataFrame
            The generated population.
        """
        next_gen_id = population.iloc[0]['generation_id'] + 1
        parents = self._selection(population, self._context)
        children = self._recombination(self._context.rng, parents, next_gen_id)
        children = self._mutation(self._context.rng, 
                                  children, 
                                  self._context.sequence_allowed_letters)
        return children



    @abstractmethod
    def compute_population_fitness(self,
                                   population: pd.DataFrame,
                                   temporary: bool = False
                                   ) -> pd.DataFrame:
        raise NotImplementedError



    @abstractmethod
    def replacement(self,
                    population: pd.DataFrame,
                    children: pd.DataFrame
                    ) -> pd.DataFrame:
        raise NotImplementedError
