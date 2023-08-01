from dataclasses import dataclass, field, asdict
from typing import List
from Bio.PDB.Atom import Atom
import evodesign.Sequence as Sequence
from evodesign.Fitness import FitnessFunction
from evodesign.Prediction import Predictor
import os





@dataclass(order=True)
class Individual:

  sequence: str = field(compare=False)
  fitness: float = field(default=None)
  metrics: dict = field(default_factory=dict, compare=False)



  @classmethod
  def new_random(cls, sequenceLength: int):
    return cls(Sequence.create_random_sequence(sequenceLength))
  


  def __float__(self):
    return self.fitness
  

  
  def as_dict(self) -> dict:
    return asdict(self)
  


  def get_pdb_filepath(self, parentFolder: str = str()) -> str:
    return os.path.join(parentFolder, f'prot_{self.sequence}.pdb')
  


  def update_fitness(self, 
                     fitnessFn: FitnessFunction, 
                     predictor: Predictor, 
                     referenceBackbone: List[Atom],
                     pdbFilename: str
                     ) -> None:
    model_backbone = predictor(self.sequence, pdbFilename)
    self.metrics, self.fitness = fitnessFn(model_backbone, referenceBackbone)
  