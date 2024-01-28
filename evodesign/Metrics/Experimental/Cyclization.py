from ..Metric import Metric





class Cyclization(Metric):

  @classmethod
  def column_name(cls) -> str:
    return 'cyclization'



  def __call__(self, kwargs) -> float:
    # distance between the N atom of the first residue and the C 
    # atom of the last
    return kwargs['model'][0] - kwargs['model'][-2]
