from .Descriptor import Descriptor





class CTDD(Descriptor):

  @classmethod
  def column_name(cls) -> str:
    return 'CTDD'
  


  def __init__(self, ilearnDir: str) -> None:
    super().__init__(ilearnDir, 'CTDD')
