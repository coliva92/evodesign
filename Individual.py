from dataclasses import dataclass, field, asdict
import evodesign.Sequence as Sequence
import os





@dataclass(order=True)
class Individual:

  sequence: str = field(compare=False)
  fitness: float = field(default=None)
  metrics: dict = field(default_factory=dict, compare=False)



  @classmethod
  def new_random(cls, sequenceLength: int):
    return cls(Sequence.create_random_sequence(sequenceLength))
  


  def as_dict(self) -> dict:
    return asdict(self)
  


  def get_pdb_filepath(self, parentFolder: str = str()) -> str:
    return os.path.join(parentFolder, f'prot_{self.sequence}.pdb')
