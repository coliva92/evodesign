from dataclasses import dataclass, field, asdict
import evodesign.Sequence as Sequence
import os





@dataclass(order=True)
class Individual:
  """
  La representación de un individuo (i.e., una secuencia de aminoácidos) en la 
  población del algoritmo evolutivo.
  """

  sequence: str = field(compare=False)
  fitness: float = field(default=None)
  metrics: dict = field(default_factory=dict, compare=False)



  @classmethod
  def random(cls, sequenceLength: int):
    return cls(Sequence.create_random_sequence(sequenceLength))
  


  def get_memento(self) -> dict:
    return asdict(self)
  


  def get_pdb_filepath(self, parentFolder: str = str()) -> str:
    return os.path.join(parentFolder, f'prot_{self.sequence}.pdb')
