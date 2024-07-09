from .FitnessFunction import FitnessFunction
from ..Metrics.Metric import Metric
import numpy as np
from typing import List, Dict, Optional, Union





class WeightedMean(FitnessFunction):
  
  @classmethod
  def column_name(cls) -> str:
    return 'fitness_weighted_mean'
  


  def _params(self) -> dict:
    params = super()._params()
    params['weights'] = self._weights
    return params
  


  def __init__(self, 
               terms: List[Metric],
               upperBound: float,
               weights: Optional[List[float]] = None
               ) -> None:
    """
    Computes the fitness of a given individual as the weighted mean of the 
    specified terms.

    Parameters
    ----------
    terms : List[Metric]
        The terms for computing the fitness value.
    upperBound : float
        The maximum allowed fitness value before triggering a stop condition
        in the evolutionary algorithm.
    weights : Optional[List[float]], optional
        The weights for each individual term. If no weights are specified,
        equal weights are used for all terms, turning this fitness function 
        into a simple mean.

    Raises
    ------
    RuntimeError
        If the length for `terms` and `weights` are not the same.
    """
    super().__init__(terms, upperBound)
    if len(weights) == 0:
      weights = len(self._terms) * [ 1.0 ]
    if len(weights) != len(self._terms):
      raise RuntimeError
    self._weights = weights
  


  def compute_fitness(self, 
                      termValues: Dict[str, Union[int, float, str]] = {}
                      ) -> float:
    """
    Uses the provided term values and the configured weights to compute the 
    final fitness value as a weighted mean.

    Parameters
    ----------
    termValues : Dict[str, any], optional
        The values of all the individual terms required for the computation 
        of the fitness value. Default is an empty dictionary. 

    Returns
    -------
    float
        The final fitness value.
    """
    x = np.array([
      termValues[term.column_name()]
      for term in self._terms
    ])
    return np.average(x, weights=self._weights)
