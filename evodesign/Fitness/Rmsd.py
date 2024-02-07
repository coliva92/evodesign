from .FitnessFunction import FitnessFunction
from ..Metrics.Rmsd import Rmsd as RmsdMetric





class Rmsd(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Rmsd'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_rmsd'
  


  @classmethod
  def upper_bound(cls) -> float:
    """
    Returns
    -------
    float
        The highest value the fitness function can achieve before triggering
        the termination of the evolutionary algorithm.
    """
    return -0.5
  


  def __init__(self) -> None:
    """
    Superimposes two given protein backbones, a model, which is movable, and 
    a reference, which remains fixed, and computes the RMSD for the model
    backbone.
    """
    super().__init__([ RmsdMetric() ])
  


  def compute_fitness(self, **kwargs) -> float:
    """
    Uses the provided metric values to compute the final fitness value.

    Returns
    -------
    float
        The RMSD value previously computed after superimposing the model 
        backbone into the reference backbone.
    """
    return kwargs['rmsd']
  