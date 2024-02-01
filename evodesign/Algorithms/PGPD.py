from abc import ABC
from .Algorithm import Algorithm
from typing import List
from ..Workspace import Workspace
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





class PGPD(Algorithm, ABC):

  @classmethod
  def _class_name(cls) -> str:
    return 'Algorithms.PGPD'
  


  def _params(self) -> dict:
    return {
      'numGenerations': self._num_generations,
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
    self._num_generations = numGenerations
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
  


  def setup(self, targetPdb: str, workspacePath: str):
    # initialize the workspace and the RNG
    super().setup(targetPdb, workspacePath)
    
    # load the target protein
    reference = Chain.load_structure(targetPdb)
    sequence_length = Chain.length(reference)
    reference = Chain.backbone_atoms(reference)

    # check if we are restoring from a previous population
    workspace = Workspace.instance()
    population = workspace.load_latest_population()
    
    data = {
      'reference': reference,
      'population': population,
      'sequence_length': sequence_length
    }
    return data
  
  

  def __call__(self, **kwargs) -> None:
    """
    Starts the execution of the evolutionary algorithm.

    Parameters
    ----------
    reference : List[Bio.PDB.Atom.Atom]
        The backbone of the target protein.
    sequence_length : int
        The number of residues in each sequence in a population.
    population : pandas.DataFrame, optional
        The population from which the evolutionary algorithm will begin its 
        execution.
    """
    reference = kwargs['reference']
    population = kwargs['population']
    sequence_length = kwargs['sequence_length']
    workspace = Workspace.instance()

    # load the statistics file
    stats = workspace.load_statistics()

    # check if we are starting from an empty population or we are continuing
    # from the last population of an earlier execution
    if not population:
      population = Population.create_random(self._pop_size, 
                                            sequence_length, 
                                            generationId=1)
      population['survivor'] = True
      workspace.save_population(population)
      workspace.save_rng_state(Random.generator().bit_generator.state)
    
    # if we are starting fresh, compute the fitness of all the initial 
    # individuals; if we are resuming from a previous population, continue
    # calculating the fitness of the individuals with missing fitness
    # (if there are any)
    population = self._compute_fitness(population, reference)
    
    # sort the population by fitness in descending order
    
    population.sort_values(by=self._fitness_fn.column_name(), 
                           ascending=False, 
                           inplace=True)
    stats = pd.concat([ stats, self._compute_statistics(population) ], 
                      ignore_index=True)

    workspace.save_population(population)
    workspace.save_statistics(stats)
    self._save_statistics_plot(stats)

    # main loop
    f = self._fitness_fn
    while True:

      # check the termination conditions
      if population.iat[0, 'generation_id'] == self._num_generations:
        break
      if population.iat[0, f.column_name()] >= f.upper_bound():
        break
      if stats.iat[-1, 'pop_sequence_identity'] >= 0.9 * sequence_length:
        break
      
      children = workspace.load_population_backup()
      if not children:
        parents = self._selection(population)
        next_gen_id = population.iat[0, 'generation_id'] + 1
        children = self._recombination(parents, next_gen_id)
        children = self._mutation(children)
        workspace.save_rng_state(Random.generator().bit_generator.state)
        workspace.save_population(children, temporary=True)
      children = self._compute_fitness(children, 
                                       reference,
                                       temporary=True)
      workspace.delete_temporary_population()

      children = self._offspring_selection(children)
      last_upper_bin = population.iloc[:self._elitism_size]
      upper_bin = children.iloc[:self._elitism_size]
      last_upper_bin['generation_id'] = children.iat[0, 'generation_id']
      upper_bin = self._merge(upper_bin, last_upper_bin)
      survivors_upper = upper_bin.iloc[self._elitism_size]
      dead_upper = upper_bin.iloc[:self._elitism_size]
      dead_upper['survivor'] = False
      objs = [ survivors_upper, children.iloc[self._elitism_size:], dead_upper ]
      population = pd.concat(objs, ignore_index=True)
      workspace.save_population(population)


  
  def _compute_fitness(self, 
                       population: pd.DataFrame,
                       reference: List[Atom],
                       **kwargs
                       ) -> pd.DataFrame:
    subset = population
    if self._fitness_fn.column_name() in population.columns:
      indices = population[self._fitness_fn.column_name()].isna()
      subset = population[indices]
    
    # this loop makes me sad :( we could use DataFrame.apply, but then we
    # wouldn't be able to catch exceptions and save the progress mid-execution
    workspace = Workspace.instance()
    try:
      for idx, row in subset.iterrows():
        filename = f'{workspace.pdbs_dir}/{row["sequence_id"]}.pdb'
        model, row['plddt'] = self._predictor(row['sequence'], filename)
        results = self._fitness_fn(model=model, 
                                   reference=reference, 
                                   sequence=row['sequence'])
        row = row.combine_first(results)
        for column, value in row.items():
          population.at[idx, column] = value
    except BaseException as e:
      workspace.save_population(population, **kwargs)
      raise e
    return population
  


  def _compute_statistics(self, population: pd.DataFrame) -> pd.Series:
    top_solution = population.iloc[0]
    top_solution['pop_missing_amino_acids'] = \
      Statistics.average_missing_amino_acids(population)
    top_solution['pop_sequence_identity'] = \
      Statistics.average_sequence_identity(population)
    return top_solution
  


  def _save_statistics_plot(self, stats: pd.DataFrame) -> None:
    workspace = Workspace.instance()
    os.makedirs(workspace.root_dir, exist_ok=True)
    fig, ax = plt.subplots(ncols=2, figsize=(14, 6))
    fig.suptitle(workspace.root_dir)
    c = self._fitness_fn.column_name()
    ax[0].plot(x=stats['generation_id'], y=stats[c], label=c)
    ax[0].plot(x=stats['generation_id'], y=stats[c], label='plddt')
    ax[0].set_ylim(bottom=0., top=1.)
    ax[0].set_xlabel('generation_id')
    ax[0].legend(loc='best')
    ax[1].plot(x=stats['generation_id'], 
               y=stats['pop_sequence_identity'], 
               color='C2',
               label='pop_sequence_identity')
    ax[1].tick_params(axis='y', labelcolor='C2')
    sequence_length = len(stats.iat[0, 'sequence'])
    ax[1].set_ylim(bottom=0., top=sequence_length)
    ax[1].set_xlabel('generation_id')
    ax2 = ax[1].twinx()
    ax2.plot(x=stats['generation_id'], 
             y=stats['pop_missing_amino_acids'], 
             color='C3',
             label='pop_missing_amino_acids')
    ax2.tick_params(axis='y', labelcolor='C3')
    ax2.set_ylim(bottom=0., top=len(Sequence.AMINO_ACIDS))
    ax[1].legend(loc='best')
    fig.savefig(f'{workspace.root_dir}/fitness_diversity.png')
  


  def _merge(self, 
             a: pd.DataFrame, 
             b: pd.DataFrame
             ) -> pd.DataFrame:
    merged = pd.DataFrame(columns=a.columns)
    i, j, n = 0, 0, len(a)
    while i < n and j < n:
      fitness_a = a.iat[i, self._fitness_fn.column_name()]
      fitness_b = b.iat[j, self._fitness_fn.column_name()]
      if fitness_a >= fitness_b:
        merged = merged.append(a.iloc[i], ignore_index=True)
        i += 1
      else:
        merged = merged.append(b.iloc[j], ignore_index=True)
        j += 1
    return merged





class _OffspringSelection:

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
    survivors, dead = pd.DataFrame(), pd.DataFrame()
    if len(children) % 2 != 0:
      children = children.iloc[:-1]
    for i in range(0, len(children), 2):
      better = children.iat[i, self._fitness_col]
      worse = children.iat[i + 1, self._fitness_col]
      b, w = (i, i + 1) if better >= worse else (i + 1, i)
      s, d = (b, w) if Random.coin_toss(self._better_bias) else (w, b)
      children.iat[s, 'survivor'] = True
      survivors = pd.concat([ survivors, pd.DataFrame(children.iloc[s]) ], 
                            ignore_index=True)
      dead = pd.concat([ dead, pd.DataFrame(children.iloc[d]) ], 
                       ignore_index=True)
    return pd.concat([ survivors, dead ], ignore_index=True)
