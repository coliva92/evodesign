from abc import ABC, abstractmethod
from typing import Dict, List
from Bio.PDB.Atom import Atom
from ..Metrics import Metric
from ..Prediction import Predictor
from ..Population import Individual
import os





class FitnessFunction(ABC):
  """
  Representación de la función de aptitud que será optimizada por el algoritmo 
  evolutivo.
  """
  
  def __init__(self, metrics: Dict[str, Metric]) -> None:
    """
    Constructor.
    - `metrics`: las métricas de calidad que serán calculadas para la 
      estructura de cada secuencia y que posteriormente serán utilizadas para 
      calcular la aptitud del individuo correspondiente.
    """
    super().__init__()
    self._metric_calculators = metrics
  


  @classmethod
  @abstractmethod
  def get_name(cls) -> str:
    pass

  

  @classmethod
  @abstractmethod
  def upper_bound(cls) -> float:
    pass


  
  @abstractmethod
  def compute_fitness(self, metrics: Dict[str, float]) -> float:
    """
    Calcula la aptitud de la secuencia especificada por `sequence`, utilizando
    los valores especificados por `metrics`.
    """
    pass



  def apply(self, 
            individual: Individual, 
            predictor: Predictor,
            reference_backbone: List[Atom],
            pdbFilename: str) -> None:
    """
    Calcula la aptitud del individuo especificado por `individual`, utilizando 
    el algoritmo de predicción de la estructura de la proteína especificado por 
    `predictor`.
    """
    backbone = predictor.get_predicted_backbone(individual.id, 
                                                individual.sequence, 
                                                pdbFilename)
    for metric, calculator in self._metric_calculators.items():
      individual.metrics[metric] = calculator.compute(backbone, 
                                                      reference_backbone)
    individual.fitness = self.compute_fitness(individual.metrics)
