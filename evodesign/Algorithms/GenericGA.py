from .Algorithm import Algorithm
from abc import ABC
from typing import Optional, List
from ..GA.Mutation.Mutation import Mutation
from ..GA.Recombination.Recombination import Recombination
from ..GA.Selection.Selection import Selection
from ..Prediction.Predictor import Predictor
from ..Fitness.FitnessFunction import FitnessFunction
from ..Sequence import Sequence
from Bio.PDB.Atom import Atom
import pandas as pd
import matplotlib.pyplot as plt
import os





class GenericGA(Algorithm, ABC):
  
  def _params(self) -> dict:
    params = super()._params()
    params['sortColumns'] = self._sort_cols
    params['sortAscending'] = self._sort_ascending
    params['fitnessFn'] = self._fitness_fn.settings()
    return params
  


  def __init__(self, 
               maxGenerations: int, 
               popSize: int, 
               predictor: Predictor, 
               fitnessFn: FitnessFunction,
               selection: Selection, 
               recombination: Recombination, 
               mutation: Mutation,
               sortColumns: Optional[List[str]] = None,
               sortAscending: Optional[List[bool]] = None
               ) -> None:
    super().__init__(maxGenerations,
                     popSize,
                     predictor,
                     selection,
                     recombination,
                     mutation)
    if not sortColumns:
      # if two individuals are tied in fitness, break the tie using the plddt
      sortColumns = [ fitnessFn.column_name(), 'plddt' ]
    self._sort_cols = sortColumns
    if not sortAscending:
      # by default sort in descending order; higher fitness comes at the top
      sortAscending = len(self._sort_cols) * [ False ]
    self._sort_ascending = sortAscending
    self._fitness_fn = fitnessFn



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
                      refSequence: Optional[str] = None,
                      **kwargs
                      ) -> pd.DataFrame:
    """
    Computes the fitness for each sequence in the given population.

    Parameters
    ----------
    population : pandas.DataFrame
        The population for which sequences the fitness will be computed.
    reference : List[Bio.PDB.Atom.Atom]
        The backbone of the target protein. The fitness of each individual
        in the population can be computed from this backbone.
    refSequence: str, optional
        The amino acid sequence of the target protein. The fitness of each
        individual in the population can be computed from this sequence.
        Default is `None`.
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
                                   refSequence=refSequence,
                                   sequence=row['sequence'],
                                   sequenceId=row['sequence_id'],
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
                                  reference: List[Atom],
                                  refSequence: Optional[str] = None
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
        The backbone of the target protein. The fitness of each individual
        in the population can be computed from this backbone.
    refSequence: str, optional
        The amino acid sequence of the target protein. The fitness of each
        individual in the population can be computed from this sequence.

    Returns
    -------
    pandas.DataFrame
        The population updated with the fitness values for all individuals.
    """
    c = self._fitness_fn.column_name()
    if c not in population.columns or population[c].isna().sum() > 0:
      population = self.compute_fitness(population, reference, refSequence)
      population.sort_values(by=self._sort_cols, 
                             ascending=self._sort_ascending,
                             inplace=True,
                             ignore_index=True)
    return population
  


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
  