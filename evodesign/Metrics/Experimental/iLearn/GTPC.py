from Descriptor import Descriptor





class GTPC(Descriptor):

  @classmethod
  def column_name(cls) -> str:
    return 'GTPC'
  


  def __init__(self, ilearnDir: str) -> None:
    super().__init__(ilearnDir, 'GTPC')
