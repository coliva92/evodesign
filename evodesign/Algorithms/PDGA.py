from .Algorithm import Algorithm
from typing import List
from ..Prediction.Predictor import Predictor
from ..Fitness.FitnessFunction import FitnessFunction
from ..GA.Selection.Selection import Selection
from ..GA.Recombination.Recombination import Recombination
from ..GA.Mutation.Mutation import Mutation
from ..Random import Random
from ..Statistics import Statistics
from ..Sequence import Sequence
from Bio.PDB.Atom import Atom
import pandas as pd
import matplotlib.pyplot as plt
import os





class PDGA(Algorithm):

  @classmethod
  def _class_name(cls) -> str:
    return 'Algorithms.PGPD'
  


  def _params(self) -> dict:
    params = super()._params()
    params['elitismSize'] = self._elitism_size
    params['betterOffspringBias'] = self._better_bias
    params['fitnessFn'] = self._fitness_fn.settings()
    return params
  


  def __init__(self,
               maxGenerations: int,
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
    maxGenerations : int
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
    super().__init__(maxGenerations,
                     popSize,
                     predictor,
                     selection,
                     recombination,
                     mutation)
    self._elitism_size = elitismSize
    self._better_bias = betterOffspringBias
    self._fitness_fn = fitnessFn
    self._offspring_selection = _OffspringSelection(fitnessFn.column_name(), 
                                                    betterOffspringBias)
  


  def compute_statistics(self, population: pd.DataFrame) -> pd.Series:
    """
    Computes the average number of missing amino acids and the average sequence 
    identity amongst the sequences in the upper bin in the given population.
    Only the upper bin is considered because, if its diversity is lost, then 
    the algorithm would've converged in a local optimum.

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
    upper_bin = population.iloc[:self._elitism_size]
    top_solution['lost_amino_acids'] = \
      Statistics.average_amino_acid_loss(upper_bin)
    top_solution['sequence_identity'] = \
      Statistics.average_sequence_identity(upper_bin)
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
                         stats['sequence_identity'], 
                         color='C2',
                         label='sequence_identity')
    ax[1].tick_params(axis='y', labelcolor='C2')
    sequence_length = len(stats.iloc[0]['sequence'])
    ax[1].set_ylim(bottom=0., top=sequence_length)
    ax[1].set_xlabel('generation_id')
    ax2 = ax[1].twinx()
    series2 = ax2.plot(stats['generation_id'], 
                       stats['lost_amino_acids'], 
                       color='C3',
                       label='lost_amino_acids')
    ax2.tick_params(axis='y', labelcolor='C3')
    ax2.set_ylim(bottom=0, top=len(Sequence.AMINO_ACIDS))
    series = series1 + series2
    labels = [ s.get_label() for s in series ]
    ax[1].legend(series, labels, loc='best')
    plt.savefig(f'{self.workspace.root_dir}/fitness_diversity.png')
    plt.close()
  


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
                                   sequence=row['sequence'],
                                   plddt=row['plddt'])
        row = row.combine_first(results)
        for column, value in row.items():
          population.at[idx, column] = value
    except BaseException as e:
      self.workspace.save_population(population, **kwargs)
      raise e
    return population
  


  def initial_fitness_computation(self, 
                                  population: pd.DataFrame,
                                  reference: List[Atom]
                                  ) -> pd.DataFrame:
    """
    Computes the fitness of the population of the first generation, and sorts
    the individuals according to their fitness in descending order. If the
    given population has the fitness values already computed for some of the
    individuals, only the individuals with missing fitness values will be
    processed.

    Parameters
    ----------
    population : pandas.DataFrame
        The population for which individuals the fitness value will be computed.
    reference : List[Bio.PDB.Atom.Atom]
        The reference backbone against which the individuals in the population
        will be compared in order to compute their fitness values.

    Returns
    -------
    pandas.DataFrame
        The population updated with the fitness values for all individuals.
    """
    c = self._fitness_fn.column_name()
    if c not in population.columns or population[c].isna().sum() > 0:
      population = self.compute_fitness(population, reference)
      population.sort_values(by=self._fitness_fn.column_name(), 
                             ascending=False,
                             inplace=True,
                             ignore_index=True)
    return population
  


  def replacement(self, 
                  population: pd.DataFrame, 
                  children: pd.DataFrame,
                  reference: List[Atom]
                  ) -> pd.DataFrame:
    # TODO documentar esta función y hacer sus pruebas unitarias
    children = self.compute_fitness(children, 
                                    reference,
                                    temporary=True)
    children.sort_values(by=self._fitness_fn.column_name(), 
                         ascending=False, 
                         inplace=True,
                         ignore_index=True)
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
  


  def termination(self, population: pd.DataFrame) -> bool:
    """
    Checks if the best solution in the given population has reached the upper 
    bound for the fitness value.

    Parameters
    ----------
    population : pandas.DataFrame
        The population for which the termination condition will be checked.

    Returns
    -------
    bool
        A flag that indicates if the fitness upper bound has been reached.
    """
    f = self._fitness_fn
    return population.iloc[0][f.column_name()] >= f.upper_bound()



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
