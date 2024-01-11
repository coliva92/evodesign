from . import FitnessFunction
from typing import Dict, Optional
from ..Metrics import SideChainPackingEnergyScore as EnergyMetric
import math





class SideChainPackingNegativeEnergyScore(FitnessFunction):

  @classmethod
  def name(cls) -> str:
    return 'Fitness_SideChainPackingNegativeEnergyScore'
  


  @classmethod
  def upper_bound(cls) -> float:
    return math.inf
  


  def __init__(self, 
               workspaceRoot: str,
               targetPdbFilename: str,
               scwrlExecutablePath: str = './scwrl4/Scwrl4') -> None:
    metrics = {
      'sideChainPackingEnergy': EnergyMetric(workspaceRoot, 
                                             targetPdbFilename, 
                                             scwrlExecutablePath)
    }
    super().__init__(metrics)
  



  def params_as_dict(self) -> dict:
    metric = self._metric_calculators['sideChainPackingEnergy']
    return {
      'workspaceRoot': metric._root_folder,
      'targetPdbFilename': metric._scwrl_input_pdb_filename,
      'scwrlExecutablePath': metric._scwrl_executable
    }
  


  def compute_fitness(self, 
                      metrics: Dict[str, float], 
                      _: Optional[str] = None) -> float:
    return -metrics['sideChainPackingEnergy']
  