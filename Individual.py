from dataclasses import dataclass, field, asdict
from typing import List, Optional
from Bio.PDB.Atom import Atom
import evodesign.Sequence as Sequence
from evodesign.Fitness import FitnessFunction
from evodesign.Prediction import Predictor
import os





@dataclass(order=True)
class Individual:

  sequence: str = field(compare=False)
  fitness: float = field(default=None, compare=True)
  metrics: dict = field(default_factory=dict, compare=False)



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
    
  
  
  def as_json(self) -> dict:
    return asdict(self)
  


  def pdb_filename(self, pdbsFolder: Optional[str] = None) -> str:
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
  