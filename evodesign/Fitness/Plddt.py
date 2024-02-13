from typing import List
from evodesign.Metrics.Metric import Metric
from .FitnessFunction import FitnessFunction





class Plddt(FitnessFunction):

  @classmethod
  def _class_name(cls) -> str:
    return 'Fitness.Plddt'
  


  @classmethod
  def column_name(cls) -> str:
    return 'fitness_plddt'
  


  def __init__(self, upperBound: float = 0.9) -> None:
    """
    Computes the average plDDT value from the all-atoms plDDT values computed
    by the structure prediction algorithm.

    Parameters
    ----------
    upperBound : float, optional
        The maximum allowed value for this fitness before triggering the
        algorithm's termination. The default is 0.9.
    """
    super().__init__(upperBound, [])
  


  def compute_fitness(self, **kwargs) -> float:
    """
    Uses the provided metric values to compute the final fitness value.

    Parameters
    ----------
    plddt: float
        The plDDT value previously computed by the structure prediction 
        algorithm.

    Returns
    -------
    float
        The final fitness value.
    """
    return kwargs['plddt']
  