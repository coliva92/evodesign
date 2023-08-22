from dataclasses import dataclass, field, asdict
from functools import total_ordering
from typing import List, Optional
from Bio.PDB.Atom import Atom
import evodesign.Sequence as Sequence
from evodesign.Fitness import FitnessFunction
from evodesign.Prediction import Predictor
import os





@total_ordering
@dataclass
class Individual:

  sequence: str
  fitness: float = field(default=None)
  metrics: dict = field(default_factory=dict)



  @classmethod
  def new_random(cls, sequenceLength: int):
    return cls(Sequence.random_sequence(sequenceLength))
  


  def __float__(self):
    return self.fitness
  


  def __len__(self):
    return len(self.sequence)
  


  def __iter__(self):
    for letter in self.sequence:
      yield letter
  


  def __getitem__(self, i: int) -> str:
    return self.sequence[i]
  


  def __setitem__(self, i: int, _) -> None:
    raise TypeError
  


  def __eq__(self, other):
    if not isinstance(other, self.__class__):
      return NotImplemented
    if self.fitness != other.fitness:
      return False
    if 'rmsd' not in self.metrics or 'rmsd' not in other.metrics:
      return True
    return self.metrics['rmsd'] == other.metrics['rmsd']
  


  def __lt__(self, other):
    if not isinstance(other, self.__class__):
      return NotImplemented
    if self.fitness == other.fitness:
      if 'rmsd' not in self.metrics or 'rmsd' not in other.metrics:
        return False
      return self.metrics['rmsd'] > other.metrics['rmsd']
    return self.fitness < other.fitness
    
  
  
  def as_json(self) -> dict:
    return asdict(self)
  


  def get_pdb_filename(self, pdbsFolder: Optional[str] = None) -> str:
    if pdbsFolder is None: pdbsFolder = ''
    return os.path.join(pdbsFolder, f'prot_{self.sequence}.pdb')
  


  def update_fitness(self, 
                     fitnessFn: FitnessFunction, 
                     predictor: Predictor, 
                     referenceBackbone: List[Atom],
                     pdbFilename: str
                     ) -> None:
    model_backbone = predictor(self.sequence, pdbFilename)
    metrics, self.fitness = fitnessFn(model_backbone, 
                                      referenceBackbone, 
                                      self.sequence)
    self.metrics = { **self.metrics, **metrics }
  