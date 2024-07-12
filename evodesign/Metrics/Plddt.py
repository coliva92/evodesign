from .Metric import Metric
from typing import Optional, Dict, List
from ..Context import Context
from Bio.PDB.Atom import Atom
import pandas as pd





class Plddt(Metric):
  
  def __init__(self, column: Optional[str] = None) -> None:
    super().__init__(column)



  def compute_value(self, 
                    backbone: List[Atom],
                    data: pd.Series,
                    context: Context
                    ) -> Dict[str, float]:
    """
    Wrapper for retrieving the predicted lDDT value when computing the fitness
    value of an individual. The predicted lDDT values is provided by the 
    structure predictor and it represents the degree of confidence that the
    predictor has about a given prediction it provided.

    Returns
    -------
    float
        The pLDDT value.

    Raises
    ------
    RuntimeError
        If the pLDDT value was not found. This could mean the predictor did 
        not returned this value or it was not passed properly to this function.
    """
    if 'otherMetrics' in kwargs and 'plddt' in kwargs['otherMetrics']:
      return kwargs['otherMetrics']['plddt']
    if 'plddt' in data.index:
      return kwargs['plddt']
    raise RuntimeError
