from .FitnessFunction import FitnessFunction
from ..Metrics.Metric import Metric
from typing import Dict, List, Optional





class LinearCombination(FitnessFunction):
  
  @classmethod
  def column_name(cls) -> str:
    return 'fitness_linear'
  


  def _params(self) -> dict:
    params = super()._params()
    params['coefficients'] = self._coefficients
    return params
  


  def __init__(self, 
               terms: List[Metric],
               upperBound: float,
               coefficients: Optional[List[float]] = None
               ) -> None:
    """
    Computes the fitness of a given individual as a linear combination of the
    specified terms.

    Parameters
    ----------
    terms : List[Metric]
        The terms for computing the fitness value.
    upperBound : float
        The maximum allowed fitness value before triggering a stop condition
        in the evolutionary algorithm.
    coefficients : Optional[List[float]], optional
        The coefficients for each individual term. If no coefficients are
        specified, all terms are assigned a coefficient of 1.0, turning this
        fitness function into a simple addition of the terms.

    Raises
    ------
    RuntimeError
        If the length of `terms` and `coefficients` are not the same.
    """
    super().__init__(upperBound, terms)
    if len(coefficients) == 0:
      coefficients = len(self._terms) * [ 1.0 ]
    if len(coefficients) != len(self._terms):
      raise RuntimeError
    self._coefficients = coefficients
  


  def compute_fitness(self, 
                      termValues: Dict[str, int | float | str] = {}
                      ) -> float:
    """
    Uses the provided term values and the configured coefficients to compute the
    final fitness value as a linear combination of said term values.

    Parameters
    ----------
    termValues : Dict[str, int  |  float  |  str], optional
        The values of all the individual terms required for the computation 
        of the fitness value. Default is an empty dictionary. 

    Returns
    -------
    float
        The final fitness value.
    """
    variables = [ termValues[term.column_name()] for term in self._terms ]
    result = 0
    for c, x in zip(self._coefficients, variables):
      result += c * x
    return result
