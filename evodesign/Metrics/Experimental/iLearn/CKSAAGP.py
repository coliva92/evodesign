from .Descriptor import Descriptor





class CKSAAGP(Descriptor):

  @classmethod
  def column_name(cls) -> str:
    return 'CKSAAGP'
  


  def __init__(self, ilearnDir: str) -> None:
    super().__init__(ilearnDir, 'CKSAAGP')
