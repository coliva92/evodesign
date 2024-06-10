from .Metric import Metric





class Cyclization(Metric):
  
  @classmethod
  def _class_name(cls) -> str:
    return 'Metrics.Cyclization'
  


  def column_name(self) -> str:
    return 'cyclization'



  def compute_value(self, **kwargs) -> float:
    # distance between the N atom of the first residue and the C 
    # atom of the last
    return kwargs['model'][0] - kwargs['model'][-2]
