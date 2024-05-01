from .Descriptor import Descriptor





class BLOSUM62(Descriptor):

  @classmethod
  def column_name(cls) -> str:
    return 'blosum62'
  


  def __init__(self, ilearnDir: str) -> None:
    super().__init__(ilearnDir, 'BLOSUM62')
