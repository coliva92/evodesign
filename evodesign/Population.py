from typing import List, Optional
from Bio.PDB.Atom import Atom
from .Individual import Individual
from .Fitness import FitnessFunction
from .Prediction import Predictor
import os
import pandas as pd
from .Sequence import Sequence





class Population:
  
  @classmethod
  def create_random(cls, 
                    size: int,
                    sequenceLength: int,
                    generationId: int = 0
                    ) -> pd.DataFrame:
    """
    Creates a collection of randomly generated amino acid sequences. 
    Populations in EvoDesign are Pandas DataFrames in which each row corresponds 
    to a single sequence, and where additional columns typically represent
    fitness values for each sequence.

    Parameters
    ----------
    size : int
        The size of the population.
    sequenceLength : int
        The length of each amino acid sequence in the population.
    generationId : int, optional
        The unique identifier for the population being created. 
        The default value is 0.

    Returns
    -------
    pd.DataFrame
        The generated population. Column `Sequence` contains the amino acid
        sequence of each individual in the population, while column 
        `Sequence_Id` contains a unique identifier for each sequence in the
        population. The `Survived` column contains a boolean flag which 
        can be used to indicate which sequences can move on to the next 
        generation.
    """
    def pad_zeroes(n: int) -> str:
      if n < 1000:
        result = f'{0}{n}'
      if n < 100:
        result = f'{0}{result}'
      if n < 10:
        result = f'{0}{result}'
    
    data = {
      'Sequence_Id': [
        f'prot_{pad_zeroes(generationId)}_{pad_zeroes(i)}'
        for i in range(size)
      ],
      'Sequence': [ 
        Sequence.create_random(sequenceLength) 
        for _ in range(size) 
      ],
      'Survived': [ False for _ in range(size) ]
    }
    return pd.DataFrame(data)



  @classmethod
  def new_random(cls,
                 size: int, 
                 sequenceLength: int,
                 iterationId: int = 1):
    return cls([ Individual.new_random(sequenceLength) for _ in range(size) ], 
               iterationId)



  def __init__(self, 
               individuals: Optional[List[Individual]] = None,
               iterationId: int = 0, 
               ) -> None:
    if individuals is None: individuals = []
    self.iteration_id = iterationId
    self.individuals = individuals
  


  def as_dict(self) -> list:
    return [ individual.as_dict() for individual in self.individuals ]
  


  def get_filename(self, populationsFolder: Optional[str] = None) -> str:
    if populationsFolder is None: populationsFolder = ''
    return os.path.join(populationsFolder, f'pop_{self.iteration_id:04d}.csv')
  


  def __len__(self) -> int:
    return len(self.individuals)
  


  def __getitem__(self, i: int) -> Individual:
    return self.individuals[i]
  


  def __setitem__(self, i: int, value) -> None:
    if not isinstance(value, Individual):
      raise NotImplementedError
    self.individuals[i] = value
  


  def __iter__(self):
    for individual in self.individuals:
      yield individual
  


  def __iadd__(self, other):
    if isinstance(other, list):
      self.individuals += other
      return self
    if isinstance(other, Population):
      self.individuals += other.individuals
      return self
    raise NotImplementedError
  


  def sort(self, reverse: bool = False) -> None:
    self.individuals.sort(reverse=reverse)



  def update_fitness(self,
                     fitnessFn: FitnessFunction, 
                     predictor: Predictor, 
                     referenceBackbone: List[Atom],
                     pdbsFolder: Optional[str] = None) -> bool:
    missing_fitness = list(filter(lambda ind: ind.fitness is None, 
                                  self.individuals))
    for individual in missing_fitness:
      filename = individual.pdb_filename(pdbsFolder)
      individual.update_fitness(fitnessFn, 
                                predictor, 
                                referenceBackbone,
                                filename)
    return len(missing_fitness) > 0
  