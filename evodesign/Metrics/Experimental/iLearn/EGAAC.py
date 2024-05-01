from .Descriptor import Descriptor





class EGAAC(Descriptor):

  @classmethod
  def column_name(cls) -> str:
    return 'EGAAC'
  


  def __init__(self, ilearnDir: str) -> None:
    super().__init__(ilearnDir, 'EGAAC')
