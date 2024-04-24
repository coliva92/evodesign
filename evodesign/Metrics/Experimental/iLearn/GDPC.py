from Descriptor import Descriptor





class GDPC(Descriptor):

  @classmethod
  def column_name(cls) -> str:
    return 'GDPC'
  


  def __init__(self, ilearnDir: str) -> None:
    super().__init__(ilearnDir, 'GDPC')
