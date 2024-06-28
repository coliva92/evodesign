from .Metric import Metric





class Plddt(Metric):
  
  def column_name(self) -> str:
    return 'plddt'
  


  def compute_value(self, **kwargs) -> float:
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
    if 'plddt' in kwargs:
      return kwargs['plddt']
    raise RuntimeError
