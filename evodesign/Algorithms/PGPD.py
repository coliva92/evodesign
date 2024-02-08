from abc import ABC
from .Algorithm import Algorithm
from typing import List
from ..Prediction.Predictor import Predictor
from ..Fitness.FitnessFunction import FitnessFunction
from ..GA.Selection.Selection import Selection
from ..GA.Recombination.Recombination import Recombination
from ..GA.Mutation.Mutation import Mutation
from ..Population import Population
from ..Chain import Chain
from ..Random import Random
from ..Statistics import Statistics
from ..Sequence import Sequence
from Bio.PDB.Atom import Atom
import pandas as pd
import matplotlib.pyplot as plt
import os
import math





# TODO actualizar el README
class PGPD(Algorithm, ABC):

  @classmethod
  def _class_name(cls) -> str:
    return 'Algorithms.PGPD'
  


  def _params(self) -> dict:
    return {
      'numGenerations': self._max_generations,
      'popSize': self._pop_size,
      'elitismSize': self._elitism_size,
      'betterOffspringBias': self._better_bias,
      'predictor': self._predictor.settings(),
      'fitnessFn': self._fitness_fn.settings(),
      'selection': self._selection.settings(),
      'recombination': self._recombination.settings(),
      'mutation': self._mutation.settings()
    }
  


  def __init__(self,
               numGenerations: int,
               popSize: int,
               elitismSize: int,
               predictor: Predictor,
               fitnessFn: FitnessFunction,
               selection: Selection,
               recombination: Recombination,
               mutation: Mutation,
               betterOffspringBias: float = 1.0,
               ) -> None:
    """
    A genetic algorithm which optimizes some fitness function obtained from
    comparing the target backbone with the model backbone, which is predicted
    using a protein structure prediction algorithm like ESMFold or AlphaFold.

    Parameters
    ----------
    numGenerations : int
        The maximum number of generations for which the algorithm will be 
        executed.
    popSize : int
        The number of sequences in the population in each generation.
    elitismSize : int
        The number of sequences with the highest fitness value in the population
        that will be treated as the 'elite' population.
    betterOffspringBias : float
        The probability for selecting the offspring with better fitness for 
        survival to the next generation. All the recombination operations that
        can work with this algorithm are assumed to always produce two children,
        for each pair of parent sequences, but only one will be considered.
        The default value is 1.0.
    predictor : Predictor
        The protein structure prediction algorithm that will be used.
    fitnessFn : FitnessFunction
        The fitness funtion that will be optimized.
    selection : Selection
        The parents selection operation that will be used.
    recombination : Recombination
        The parents recombination operation that will be used.
    mutation : Mutation
        The sequence mutation operation that will be used.
    """
    super().__init__()
    self._max_generations = numGenerations
    self._pop_size = popSize
    self._elitism_size = elitismSize
    self._better_bias = betterOffspringBias
    self._predictor = predictor
    self._fitness_fn = fitnessFn
    self._selection = selection
    self._recombination = recombination
    self._mutation = mutation
    self._offspring_selection = _OffspringSelection(fitnessFn.column_name(), 
                                                    betterOffspringBias)

  
  
  def __call__(self, 
               reference: List[Atom], 
               population: pd.DataFrame, 
               **kwargs
               ) -> None:
    """
    Starts the execution of the evolutionary algorithm.

    Parameters
    ----------
    reference : List[Bio.PDB.Atom.Atom]
        The backbone of the target protein.
    population : pandas.DataFrame, optional
        The population from which the evolutionary algorithm will begin its 
        execution.
    numGenerations : int, optional
        If a value greater than zero is provided, the algorithm will be run 
        for the given amount of generations, unless any other termination
        condition is met. If a zero or lower value is provided, then the
        algorithm will until the configured maximum number of generations
        is met. The default is 0 (run until the maximum is met). 
    
    Raises
    ------
    Exceptions.HttpBadRequest
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with HTTP error code 
        400.
    Exceptions.HttpForbidden
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with HTTP error code 
        403.
    Exceptions.HttpInternalServerError
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with HTTP error code
        500.
    Exceptions.HttpGatewayTimeout
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with HTTP error code
        504.
    Exceptions.HttpUnknownError
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with any other HTTP
        error code not listed here.
    requests.exceptions.ConnectTimeout
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when communication attempt with the remote
        server times out.
    KeyboardInterrupt
        Raised any time the user presses Ctrl + C from the standard input.
    """
    num_generations = kwargs['numGenerations'] \
                      if kwargs['numGenerations'] > 0 \
                      else math.inf
    sequence_length = len(reference) // len(Chain.BACKBONE_ATOMS)
    i = 0

    # store the target PDB in the workspace
    self.workspace.save_target_pdb()

    # load the statistics file
    stats = self.workspace.load_statistics()

    # check if we are starting from an empty population or we are continuing
    # from the last population of an earlier execution
    if population.empty:
      population = self.initial_population(sequence_length)
      self.workspace.save_population(population)
      self.workspace.save_rng_state(Random.generator().bit_generator.state)
      i = 1
    
    c = self._fitness_fn.column_name()
    if c not in population.columns or population[c].isna().sum() > 0:
      # if we are starting fresh, compute the fitness of all the initial 
      # individuals; if we are resuming from a previous population, continue
      # calculating the fitness of the individuals with missing fitness
      # (if there are any)
      population = self.compute_fitness(population, reference)
      
      # sort the population by fitness in descending order
      population.sort_values(by=self._fitness_fn.column_name(), 
                            ascending=False, 
                            inplace=True)
    
    # compute the statistics if we have not done it yet
    curr_gen_id = population.iloc[0]['generation_id']
    if stats.empty or stats.iloc[-1]['generation_id'] != curr_gen_id:
      new_stats = self.compute_statistics(population)
      stats = pd.concat([ stats, new_stats.to_frame().T ], ignore_index=True)
    
    # save progress
    self.workspace.save_population(population)
    self.workspace.save_statistics(stats)
    self.save_statistics_graph(stats)

    # main loop
    f = self._fitness_fn
    while True:
      # check the termination conditions
      if i >= num_generations:
        break
      if population.iloc[0]['generation_id'] == self._max_generations:
        break
      if population.iloc[0][f.column_name()] >= f.upper_bound():
        break
      if stats.iloc[-1]['pop_sequence_identity'] >= 0.9 * sequence_length:
        break

      # check if we are resuming from an earlier execution and the children were
      # already generated; if not, apply the evolutionary operations
      children = self.workspace.load_population(temporary=True)
      if children.empty:
        children = self.next_population(population)
        self.workspace.save_rng_state(Random.generator().bit_generator.state)
        self.workspace.save_population(children, temporary=True)
      children = self.compute_fitness(children, 
                                      reference,
                                      temporary=True)
      children.sort_values(by=self._fitness_fn.column_name(), 
                           ascending=False, 
                           inplace=True)
      population = self.replacement(population, children)
      new_stats = self.compute_statistics(population)
      stats = pd.concat([ stats, new_stats.to_frame().T ], ignore_index=True)
      self.workspace.save_population(population)
      self.workspace.save_statistics(stats)
      self.save_statistics_graph(stats)
      self.workspace.delete_temporary_population()
      i += 1
  


  def initial_population(self, sequenceLength: int) -> pd.DataFrame:
    """
    Creates a collection of random amino acid sequences to be used as the 
    population for the very first generation of the evolutionary algorithm.

    Parameters
    ----------
    sequenceLength : int
        The number of amino acid residues that each sequence in the generated 
        population should have.

    Returns
    -------
    pandas.DataFrame
        The generated population data.
    """
    population = Population.create_random(self._pop_size, 
                                          sequenceLength, 
                                          generationId=1)
    population['survivor'] = True
    return population
  


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
    parents = self._selection(population)
    children = self._recombination(parents, next_gen_id)
    children = self._mutation(children)
    return children



  def compute_fitness(self, 
                      population: pd.DataFrame,
                      reference: List[Atom],
                      **kwargs
                      ) -> pd.DataFrame:
    """
    Computes the fitness for each sequence in the given population.

    Parameters
    ----------
    population : pandas.DataFrame
        The population for which sequences the fitness will be computed.
    reference : List[Bio.PDB.Atom.Atom]
        The target backbone against which each sequence in the population will 
        be compared in order to compute its fitness.
    temporary : bool, optional
        Indicates if the data of the current population data should be stored in 
        the temporary storage location or the permanent one in case an exception
        is raised during the fitness computation process.

    Returns
    -------
    pandas.DataFrame
        The original population with added columns to include the computed 
        fitness value and other associated results.

    Raises
    ------
    Exceptions.HttpBadRequest
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with HTTP error code 
        400.
    Exceptions.HttpForbidden
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with HTTP error code 
        403.
    Exceptions.HttpInternalServerError
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with HTTP error code
        500.
    Exceptions.HttpGatewayTimeout
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with HTTP error code
        504.
    Exceptions.HttpUnknownError
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when the server responds with any other HTTP
        error code not listed here.
    requests.exceptions.ConnectTimeout
        Raised only if the predictor used is 'Prediction.ESMFold.RemoteApi'.
        This exception is raised when communication attempt with the remote
        server times out.
    KeyboardInterrupt
        Raised any time the user presses Ctrl + C from the standard input.
    """
    subset = population
    if self._fitness_fn.column_name() in population.columns:
      indices = population[self._fitness_fn.column_name()].isna()
      subset = population[indices]
    
    # this loop makes me sad :( we could use DataFrame.apply, but then we
    # wouldn't be able to catch exceptions and save the progress mid-loop
    try:
      for idx, row in subset.iterrows():
        filename = f'{self.workspace.pdbs_dir}/{row["sequence_id"]}.pdb'
        model, row['plddt'] = self._predictor(row['sequence'], filename)
        results = self._fitness_fn(model=model, 
                                   reference=reference, 
                                   sequence=row['sequence'])
        row = row.combine_first(results)
        for column, value in row.items():
          population.at[idx, column] = value
    except BaseException as e:
      self.workspace.save_population(population, **kwargs)
      raise e
    return population
  


  def compute_statistics(self, population: pd.DataFrame) -> pd.Series:
    """
    Computes the average number of missing amino acids and the average sequence 
    identity amongst all the sequences in the given population.

    Parameters
    ----------
    population : pandas.DataFrame
        The population from which the statistics will be computed.

    Returns
    -------
    pandas.Series
        The information of the top solution in the given population, alongside
        the computed statistics for said population.
    """
    top_solution = population.iloc[0].copy()
    top_solution['pop_missing_amino_acids'] = \
      Statistics.average_missing_amino_acids(population)
    top_solution['pop_sequence_identity'] = \
      Statistics.average_sequence_identity(population)
    return top_solution
  


  def save_statistics_graph(self, stats: pd.DataFrame) -> None:
    """
    Creates a graph for the given statistics and saves it as a PNG file in the
    workspace.

    Parameters
    ----------
    stats : pandas.DataFrame
        The statistics which data will be graphed and stored.
    """
    os.makedirs(self.workspace.root_dir, exist_ok=True)
    plt.ioff() # turn off the GUI
    fig, ax = plt.subplots(ncols=2, figsize=(14, 6))
    fig.suptitle(self.workspace.root_dir)
    c = self._fitness_fn.column_name()
    series1 = ax[0].plot(stats['generation_id'], stats[c], color='C0', label=c)
    ax[0].tick_params(axis='y', labelcolor='C0')
    ax[0].set_xlabel('generation_id')
    ax1 = ax[0].twinx()
    series2 = ax1.plot(stats['generation_id'], 
                       stats['plddt'], 
                       color='C1', 
                       label='plddt')
    ax1.tick_params(axis='y', labelcolor='C1')
    series = series1 + series2
    labels = [ s.get_label() for s in series ]
    ax[0].legend(series, labels, loc='best')
    series1 = ax[1].plot(stats['generation_id'], 
                         stats['pop_sequence_identity'], 
                         color='C2',
                         label='pop_sequence_identity')
    ax[1].tick_params(axis='y', labelcolor='C2')
    sequence_length = len(stats.iloc[0]['sequence'])
    ax[1].set_ylim(bottom=0., top=sequence_length)
    ax[1].set_xlabel('generation_id')
    ax2 = ax[1].twinx()
    series2 = ax2.plot(stats['generation_id'], 
                       stats['pop_missing_amino_acids'], 
                       color='C3',
                       label='pop_missing_amino_acids')
    ax2.tick_params(axis='y', labelcolor='C3')
    ax2.set_ylim(bottom=0., top=len(Sequence.AMINO_ACIDS))
    series = series1 + series2
    labels = [ s.get_label() for s in series ]
    ax[1].legend(series, labels, loc='best')
    plt.savefig(f'{self.workspace.root_dir}/fitness_diversity.png')
    plt.close()
  


  def replacement(self, 
                  population: pd.DataFrame, 
                  children: pd.DataFrame
                  ) -> pd.DataFrame:
    # TODO documentar esta función y hacer sus pruebas unitarias
    children = self._offspring_selection(children)
    last_upper_bin = population.iloc[:self._elitism_size].copy()
    upper_bin = children.iloc[:self._elitism_size]
    last_upper_bin['generation_id'] = children.iloc[0]['generation_id']
    upper_bin = self._merge(upper_bin, last_upper_bin)
    survivors_upper = upper_bin.iloc[:self._elitism_size]
    dead_upper = upper_bin.iloc[self._elitism_size:].copy()
    dead_upper['survivor'] = False
    objs = [ 
      survivors_upper, 
      children.iloc[self._elitism_size:],
      dead_upper
    ]
    results = pd.concat(objs, axis=0, ignore_index=True)
    return results



  def _merge(self, 
             a: pd.DataFrame, 
             b: pd.DataFrame
             ) -> pd.DataFrame:
    # TODO documentar esta funcion y hacer las pruebas unitarias
    merged = pd.DataFrame(columns=a.columns)
    i, j, n = 0, 0, len(a)
    while i < n and j < n:
      fitness_a = a.iloc[i][self._fitness_fn.column_name()]
      fitness_b = b.iloc[j][self._fitness_fn.column_name()]
      if fitness_a >= fitness_b:
        merged = pd.concat([ merged, a.iloc[i].to_frame().T ],
                           ignore_index=True)
        i += 1
      else:
        merged = pd.concat([ merged, b.iloc[j].to_frame().T ],
                           ignore_index=True)
        j += 1
    while i < n:
      merged = pd.concat([ merged, a.iloc[i].to_frame().T ],
                         ignore_index=True)
      i += 1
    while j < n:
      merged = pd.concat([ merged, b.iloc[j].to_frame().T ],
                          ignore_index=True)
      j += 1
    return merged





class _OffspringSelection:
  # TODO documentar esta clase y hacer sus pruebas unitarias

  def __init__(self,
               fitnessColumn: str,
               betterFitnessBias: float = 1.0
               ) -> None:
    self._fitness_col = fitnessColumn
    self._better_bias = betterFitnessBias
  


  def __call__(self, children: pd.DataFrame) -> pd.DataFrame:
    # we assume that the recombination operation used produced two children
    # per parent pair; thus, if the number of children is odd, then disregard
    # the last element
    if len(children) % 2 != 0:
      children = children.iloc[:-1]
    survivors, dead = pd.DataFrame(), pd.DataFrame()
    for i in range(0, len(children), 2):
      better = children.iloc[i][self._fitness_col]
      worse = children.iloc[i + 1][self._fitness_col]
      b, w = (i, i + 1) if better >= worse else (i + 1, i)
      s, d = (b, w) if Random.coin_toss(self._better_bias) else (w, b)
      selected_child = children.iloc[s].copy()
      selected_child['survivor'] = True
      survivors = pd.concat([ survivors, selected_child.to_frame().T ],
                            ignore_index=True)
      dead = pd.concat([ dead, children.iloc[d].to_frame().T ],
                       ignore_index=True)
    results = pd.concat([ survivors, dead ], ignore_index=True)
    return results
