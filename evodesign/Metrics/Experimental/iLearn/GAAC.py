from .Descriptor import Descriptor





class GAAC(Descriptor):

  @classmethod
  def column_name(cls) -> str:
    return 'GAAC'
  


  def __init__(self, ilearnDir: str) -> None:
    super().__init__(ilearnDir, 'GAAC')
