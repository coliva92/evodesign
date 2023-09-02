from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Optional
from Bio.PDB.Atom import Atom
from ..Metrics import Metric





class FitnessFunction(ABC):

  @classmethod
  @abstractmethod
  def name(cls) -> str:
    raise NotImplementedError

  

  @classmethod
  @abstractmethod
  def upper_bound(cls) -> float:
    raise NotImplementedError
  
  
  
  def __init__(self, metrics: Dict[str, Metric]) -> None:
    super().__init__()
    self._metric_calculators = metrics
  


  def params_json(self) -> dict:
    return dict()


  
  @abstractmethod
  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    raise NotImplementedError



  def __call__(self,
               modelBackbone: List[Atom], 
               referenceBackbone: List[Atom],
               sequence: Optional[str] = None
               ) -> Tuple[Dict[str, float], float]:
    metrics = { 
      key: calc(modelBackbone, referenceBackbone, sequence) \
        for key, calc in self._metric_calculators.items() 
    }
    fitness = self.compute_fitness(metrics)
    return metrics, fitness
