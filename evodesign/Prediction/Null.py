from Predictor import Predictor
from typing import List
from Bio.PDB.Atom import Atom





class Null(Predictor):
  
  @classmethod
  def _class_name(cls) -> str:
    return 'Prediction.None'
  


  def predict_structure(self, 
                        sequence: str, 
                        pdbPath: str) -> None:
    """
    Does nothing.

    Parameters
    ----------
    sequence : str
        Unused.
    pdbPath : str
        Unused.
    """
    pass



  def __call__(self,
               sequence: str, 
               pdbPath: str
               ) -> List[Atom]:
    """
    Does nothing.

    Parameters
    ----------
    sequence : str
        Unused.
    pdbPath : str
        Unused.

    Returns
    -------
    List[Atom]
        An empty list.
    """
    return []
