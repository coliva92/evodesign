from .FitnessFunction import FitnessFunction
from typing import List
from ..Metrics.Rmsd import Rmsd
from ..Metrics.Gdt import Gdt as GdtMetric





class Gdt(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Gdt'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_gdt'
  


  def _params(self) -> dict:
    return { 'cutoffs': self._cutoffs }
  


  @classmethod
  def upper_bound(cls) -> float:
    """
    Returns
    -------
    float
        The highest value the fitness function can achieve in order to trigger
        the termination of the evolutionary algorithm.
    """
    return 0.95
  


  def __init__(self,
               cutoffs: List[float] = [ 1.0, 2.0, 4.0, 8.0 ],
               ) -> None:
    """
    Superimposes two given protein backbones, a model, which is movable, and 
    a reference, which remains fixed, and computes the GDT for the model
    backbone.

    Parameters
    ----------
    cutoffs : List[float], optional
        The cutoff distances to use (measured in Angstroms) to obtain the GDT. 
        
        Each backbone atom in the model that deviates from the position of the corresponding atom
        in the reference by a distance that is no greater than at least one
        of the given cutoff distances contributes to a positive increase
        in the final GDT value.
        
        The default configuration is [ 1.0, 2.0, 4.0, 8.0 ].
    """
    super().__init__([ Rmsd(), GdtMetric(cutoffs) ])
    self._cutoffs = cutoffs
  


  def compute_fitness(self, **kwargs) -> float:
    """
    Uses the provided metric values to compute the final fitness value.

    Parameters
    ----------
    gdt: float
        The GDT value previously computed after superimposing the model 
        backbone into the reference backbone.

    Returns
    -------
    float
        The final fitness value.
    """
    return kwargs['gdt']
