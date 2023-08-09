from typing import List, Optional
from Bio.PDB.Atom import Atom
from .Individual import Individual
from .Fitness import FitnessFunction
from .Prediction import Predictor
import os





class Population:
  
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
    self.iterationId = iterationId
    self.individuals = individuals
  


  def as_json(self) -> list:
    return [ individual.as_json() for individual in self.individuals ]
  


  def get_filename(self, populationsFolder: Optional[str] = None) -> str:
    if populationsFolder is None: populationsFolder = ''
    return os.path.join(populationsFolder, f'pop_{self.iterationId}.csv')
  


  def __len__(self) -> int:
    return len(self.individuals)
  


  def __getitem__(self, i: int) -> Individual:
    return self.individuals[i]
  


  def __setitem__(self, i: int, individual: Individual) -> None:
    self.individuals[i] = individual
  


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



  def update_fitness(self,
                     fitnessFn: FitnessFunction, 
                     predictor: Predictor, 
                     referenceBackbone: List[Atom],
                     pdbsFolder: Optional[str] = None) -> bool:
    missing_fitness = list(filter(lambda ind: ind.fitness is None, 
                                  self.individuals))
    for individual in missing_fitness:
      filename = individual.get_pdb_filename(pdbsFolder)
      individual.update_fitness(fitnessFn, 
                                predictor, 
                                referenceBackbone,
                                filename)
    if missing_fitness:
      self.individuals = sorted(self.individuals)
      return True
    return False
  