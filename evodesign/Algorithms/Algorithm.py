from abc import ABC, abstractmethod
from ..SettingsRetrievable import SettingsRetrievable
from ..Workspace import Workspace
import evodesign.Random as Random
import evodesign.Chain as Chain
import evodesign.Population as Population
from ..Prediction.Predictor import Predictor
from ..GA.Selection.Selection import Selection
from ..GA.Recombination.Recombination import Recombination
from ..GA.Mutation.Mutation import Mutation
from Bio.PDB.Atom import Atom
import pandas as pd
from typing import Tuple, List, Optional
import math





class Algorithm(SettingsRetrievable, ABC):

  def _params(self) -> dict:
    return {
      'maxGenerations': self._max_generations,
      'popSize': self._pop_size,
      'predictor': self._predictor.settings(),
      'selection': self._selection.settings(),
      'recombination': self._recombination.settings(),
      'mutation': self._mutation.settings()
    }
  


  def __init__(self,
               maxGenerations: int,
               popSize: int,
               predictor: Predictor,
               selection: Selection,
               recombination: Recombination,
               mutation: Mutation
               ) -> None:
    """
    The generic workflow of a genetic algorithm for protein design.

    Parameters
    ----------
    maxGenerations : int
        The maximum number of generations for which the algorithm will be 
        executed.
    popSize : int
        The number of sequences in the population in each generation.
    predictor : Predictor
        The protein structure prediction algorithm that will be used.
    selection : Selection
        The parents selection operation that will be used.
    recombination : Recombination
        The parents recombination operation that will be used.
    mutation : Mutation
        The sequence mutation operation that will be used.
    """
    self._max_generations = maxGenerations
    self._pop_size = popSize
    self._predictor = predictor
    self._selection = selection
    self._recombination = recombination
    self._mutation = mutation
  


  def setup(self,
            targetPdbPath: str,
            workspaceDir: str,
            targetFastaPath: Optional[str] = None
            ) -> Tuple[List[Atom], pd.DataFrame, str | None]:
    """
    Initializes the workspace and the RNG, as well as the reference backbone
    and the population before running the evolutionary algorithm.

    Parameters
    ----------
    targetPdbPath : str
        The path to the PDB file containing the target protein backbone.
    workspaceDir : str
        The folder where all the output files generated by the evolutionary
        algorithm will be stored.
    targetFastaPath : str, optional
        The path to the FASTA file containing the amino acid sequence of the 
        target protein. Default is `None`.
    
    Returns
    -------
    Tuple[List[Bio.PDB.Atom.Atom], pandas.DataFrame, str | None]
        The reference backbone and the population from which the algorithm
        will begin its execution. If a value was provided for `targetFasta`,
        then a string containing the amino acid sequence loaded from this
        file is also returned. Otherwise, `None` is returned as a the third 
        value.
    """
    # initialize the workspace
    self.workspace = Workspace(workspaceDir, targetPdbPath, targetFastaPath)

    # initialize the RNG
    rng = Random.generator()
    state = self.workspace.load_rng_state()
    if state == None:
      state = self.workspace.load_rng_state(loadCheckpoint=False)
    if state != None:
      rng.bit_generator.state = state
    else:
      self.workspace.save_rng_state(rng.bit_generator.state, checkpoint=False)
    
    # save the algorithm's settings
    self.workspace.save_settings(self.settings())

    # load the target protein
    reference = Chain.load_structure(targetPdbPath)
    reference = Chain.backbone_atoms(reference)

    # check if we are restoring from a previous population
    population = self.workspace.load_population()

    # load the target sequence if one was provided
    ref_sequence = None
    if targetFastaPath != None:
      ref_sequence = self.workspace.load_target_fasta()
    return reference, population, ref_sequence

  

  def __call__(self, 
               reference: List[Atom], 
               population: pd.DataFrame,
               refSequence: Optional[str] = None,
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
    refSequence : str, optional
        The amino acid sequence of the target protein. Default is `None`. 
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
    # set the limit for how many generations execute
    num_generations = math.inf
    if 'numGenerations' in kwargs and kwargs['numGenerations'] > 0:
      num_generations = kwargs['numGenerations']
    
    # compute the sequence length
    sequence_length = len(reference) // len(Chain.BACKBONE_ATOMS)
    
    # set the generations counter
    t = 0

    # store the target PDB in the workspace
    self.workspace.save_target_pdb()

    # store the target FASTA in the workspace
    self.workspace.save_target_fasta()

    # load the statistics file
    stats = self.workspace.load_statistics()

    # check if we are starting from an empty population or we are continuing
    # from the last population of an earlier execution
    if population.empty:
      population = self.initial_population(sequence_length)
      self.workspace.save_population(population)
      self.workspace.save_rng_state(Random.generator().bit_generator.state)
      t = 1
    
    # if we are starting fresh, compute the fitness of all the initial 
    # individuals; if we are resuming from a previous population, continue
    # calculating the fitness of the individuals with missing fitness
    # (if there are any)
    population = self.initial_fitness_computation(population, 
                                                  reference, 
                                                  refSequence)

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
    while True:
      # check the termination conditions
      if t == num_generations:
        break
      if population.iloc[0]['generation_id'] == self._max_generations:
        break
      if stats.iloc[-1]['sequence_identity'] >= 0.95 * sequence_length:
        break
      if self.termination(population):
        break

      # check if we are resuming from an earlier execution and the children were
      # already generated; if not, apply the evolutionary operations
      children = self.workspace.load_population(temporary=True)
      if children.empty:
        children = self.next_population(population)
        self.workspace.save_rng_state(Random.generator().bit_generator.state)
        self.workspace.save_population(children, temporary=True)
      
      # compute the children's fitness and choose the surviving children
      children = self.compute_fitness(children, 
                                      reference,
                                      refSequence,
                                      temporary=True)
      population = self.replacement(population, children)
      new_stats = self.compute_statistics(population)
      stats = pd.concat([ stats, new_stats.to_frame().T ], ignore_index=True)

      # save progress
      self.workspace.save_population(population)
      self.workspace.save_statistics(stats)
      self.save_statistics_graph(stats)
      self.workspace.delete_temporary_population()
      t += 1
  


  @abstractmethod
  def compute_statistics(self, population: pd.DataFrame) -> pd.Series:
    raise NotImplementedError
  


  @abstractmethod
  def save_statistics_graph(self, stats: pd.DataFrame) -> None:
    raise NotImplementedError
  


  @abstractmethod
  def compute_fitness(self, 
                      population: pd.DataFrame,
                      reference: List[Atom],
                      refSequence: Optional[str] = None,
                      **kwargs
                      ) -> pd.DataFrame:
    raise NotImplementedError
  


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
                                          generation_id=1)
    population['survivor'] = True
    return population
  


  @abstractmethod
  def initial_fitness_computation(self, 
                                  population: pd.DataFrame,
                                  reference: List[Atom],
                                  refSequence: Optional[str] = None
                                  ) -> pd.DataFrame:
    raise NotImplementedError
  


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
  


  @abstractmethod
  def replacement(self, 
                  population: pd.DataFrame, 
                  children: pd.DataFrame
                  ) -> pd.DataFrame:
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
